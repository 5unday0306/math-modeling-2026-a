function varargout=q4_compute(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q4:ComputeAction','Invalid compute action.');
switch action
    case "parameters"
        varargout{1}=q4_parameters(varargin{:});
    case "load_inputs"
        varargout{1}=load_q4_inputs(varargin{:});
    case "select_suffix"
        varargout{1}=select_q4_stable_suffix(varargin{:});
    case "choose_suffix"
        [varargout{1:nargout}]=q4_choose_suffix(varargin{:});
    case "build_mesh"
        varargout{1}=build_q4_mesh(varargin{:});
    case "interpolate_ambient"
        varargout{1}=interpolate_q4_ambient(varargin{:});
    case "interpolate_radius"
        varargout{1}=interpolate_q4_radius(varargin{:});
    case "material_properties"
        varargout{1}=q4_material_properties(varargin{:});
    case "face_coefficients"
        varargout{1}=q4_face_coefficients(varargin{:});
    case "residual_jacobian"
        [varargout{1:nargout}]=q4_residual_jacobian(varargin{:});
    case "advance"
        [varargout{1:nargout}]=advance_q4_implicit(varargin{:});
    case "reconstruct_fields"
        varargout{1}=reconstruct_q4_fields(varargin{:});
    case "spatial_maximum"
        varargout{1}=q4_spatial_maximum(varargin{:});
    case "progress_metrics"
        varargout{1}=q4_progress_metrics(varargin{:});
    case "solve_longtime"
        varargout{1}=solve_q4_longtime(varargin{:});
    case "locate_event"
        varargout{1}=locate_q4_event(varargin{:});
    case "extract_output_grid"
        varargout{1}=extract_q4_output_grid(varargin{:});
    case "solve_case"
        varargout{1}=q4_solve_case(varargin{:});
    otherwise
        error('CUMCM:q4:ComputeAction','Unknown compute action: %s',action);
end
end

function p=q4_parameters(cfg,caseSpec)
if nargin<2,caseSpec=struct();end
p=struct('appendix',4,'geometry',"pchip",'environment',"unified_stable_suffix", ...
    'rho0',760,'rhoC',90,'cp0',1850,'cpC',2150,'k0',0.12,'kC',0.20, ...
    'D0',4.2e-4,'a',0.30,'B',3850,'hc',25,'hm',8e-7, ...
    'C0',2.55,'T0_K',301.15,'R0_m',0.02,'length_m',0.25, ...
    'threshold',0.15,'N',20,'initialStep_s',0.1,'maxStep_s',60, ...
    'localTolerance',[5e-6,5e-7],'newtonTolerance',[5e-9,5e-10], ...
    'newtonMax',30,'backtrackMin',2^-20,'armijo',1e-4,'quadratureLevel',0, ...
    'rootWidth_s',0.03,'stopLattice_s',0.36,'fieldBudget',5e-5, ...
    'wallSeconds',1800,'eventStop',false,'outputInterval_s',60);
p.selectionConfig=struct('temperatureThreshold',50,'moistureThreshold',0.05, ...
    'minimumDuration_min',30,'lambda',0.8,'driftMax',0.5,'blockSpanMax',0.5, ...
    'varianceRatioMax',3,'tieFactor',64);
names=fieldnames(caseSpec);
for k=1:numel(names),p.(names{k})=caseSpec.(names{k});end
if p.appendix==3
    p.rho0=650;p.rhoC=128;p.cp0=1450;p.cpC=2736;p.k0=.21;p.kC=.38;p.D0=.0024;p.a=.45;
end
assert(ismember(p.appendix,[3,4]));
p.projectRoot=cfg.projectRoot;
end


function in=load_q4_inputs(cfg)
paths=[fullfile(cfg.rawDataDir,'附件1.xlsx'),fullfile(cfg.rawDataDir,'附件2.xlsx')];
archive=java.util.zip.ZipFile(char(paths(1)));zipCleanup=onCleanup(@()archive.close());
stream=archive.getInputStream(archive.getEntry('xl/worksheets/sheet1.xml'));
scanner=java.util.Scanner(stream,'UTF-8');scanCleanup=onCleanup(@()scanner.close());
scanner.useDelimiter('\A');xml=char(scanner.next());
assert(isempty(regexp(xml,'<(?:[A-Za-z0-9_]+:)?f(?:\s|>)','once')), ...
    'CUMCM:q4:SelectionInput','Ambient observations must not contain Excel formulas.');
a=readmatrix(paths(1),'Sheet','Sheet1','Range','A2:C242','UseExcel',false);
r=readmatrix(paths(2),'Sheet','Sheet1','Range','A2:B146','UseExcel',false);
assert(isequal(size(a),[241,3])&&all(isfinite(a),'all')&&all(diff(a(:,1))==60));
assert(isequal(size(r),[145,2])&&all(isfinite(r),'all')&&all(diff(r(:,1))==1800));
assert(a(1,1)==0&&a(end,1)==14400&&r(1,1)==0&&r(end,1)==259200);
assert(all(r(:,2)>0)&&all(diff(r(:,2))<=0)&&r(1,2)==2);
h1=readcell(paths(1),'Sheet','Sheet1','Range','A1:C1','UseExcel',false);
h2=readcell(paths(2),'Sheet','Sheet1','Range','A1:B1','UseExcel',false);
assert(isequal(string(h1),["时间","温度","水分浓度"])&&isequal(string(h2),["时间","半径"]));
in=struct('ambient',a,'radius',r,'radius_pp',pchip(r(:,1),.01*r(:,2)));
[breaks,coefs,~,order]=unmkpp(in.radius_pp);
in.radius_derivative_pp=mkpp(breaks,coefs(:,1:end-1).*(order-1:-1:1));
p=q4_parameters(cfg);in.environment_selection=select_q4_stable_suffix(a,p.selectionConfig);
fid=fopen(paths(1),'rb');assert(fid>0);cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');
md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
in.environment_selection.selection.input=struct('path','data/raw/附件1.xlsx','sha256',hash,'bytes',numel(bytes));
end


function stats=select_q4_stable_suffix(ambient,cfg)
assert(isnumeric(ambient)&&isreal(ambient)&&isequal(size(ambient),[241,3])&& ...
    all(isfinite(ambient),'all')&&isequal(ambient(:,1),(0:60:14400)'), ...
    'CUMCM:q4:SelectionInput','Expected 241 finite one-minute observations.');
assert(cfg.minimumDuration_min==30&&cfg.tieFactor==64&& ...
    cfg.lambda>0&&cfg.lambda<1&&all(isfinite([cfg.temperatureThreshold, ...
    cfg.moistureThreshold,cfg.driftMax,cfg.blockSpanMax,cfg.varianceRatioMax]))&& ...
    all([cfg.driftMax,cfg.blockSpanMax,cfg.varianceRatioMax]>0), ...
    'CUMCM:q4:SelectionConfig','Invalid selection configuration.');
stats.temperature=oneSignal(ambient(:,2),cfg.temperatureThreshold,"temperature","degC",cfg);
stats.moisture=oneSignal(ambient(:,3),cfg.moistureThreshold,"moisture","kg/kg",cfg);
stats.windows=[stats.temperature.windows;stats.moisture.windows];
stats.selection=struct('scenario','unified_stable_suffix','configuration',cfg, ...
    'temperature',stats.temperature.selection,'moisture',stats.moisture.selection, ...
    'temperature_C',stats.temperature.selection.mean_value, ...
    'temperature_K',stats.temperature.selection.mean_value+273.15, ...
    'moisture_kgkg',stats.moisture.selection.mean_value, ...
    'temperature_jump_C',stats.temperature.selection.mean_value-ambient(end,2), ...
    'moisture_jump_kgkg',stats.moisture.selection.mean_value-ambient(end,3), ...
    'temperature_sensitivity',stats.temperature.sensitivity, ...
    'moisture_sensitivity',stats.moisture.sensitivity, ...
    'history_through_4h_unchanged',true,'physical_equilibrium_onset_verified',false, ...
    'prediction_improvement_verified',false);
end

function result=oneSignal(x,threshold,label,unit,cfg)
first=find(x>threshold,1)-1;
assert(~isempty(first),'CUMCM:q4:NoExceedance','No strictly exceeding observation.');
last=240-cfg.minimumDuration_min;
assert(first<=last,'CUMCM:q4:NoCandidates','Threshold crossing leaves no 30-minute suffix.');
starts=(first:last)';count=numel(starts);rows=cell(count,1);
for ii=1:count
    k=starts(ii);v=x(k+1:end);n=numel(v);mu=mean(v);variance=sum((v-mu).^2)/(n-1);
    assert(isfinite(mu)&&isfinite(variance)&&variance>0&&any(v~=v(1)), ...
        'CUMCM:q4:DegenerateVariance','Nonpositive or nonfinite candidate variance.');
    tm=(k:240)';center=tm-mean(tm);beta=sum(center.*(v-mu))/sum(center.^2);
    sizes=repmat(floor(n/3),1,3);sizes(1:mod(n,3))=sizes(1:mod(n,3))+1;
    bm=zeros(1,3);bv=bm;lo=1;
    for b=1:3
        z=v(lo:lo+sizes(b)-1);bm(b)=mean(z);bv(b)=sum((z-bm(b)).^2)/(sizes(b)-1);
        if all(z==z(1)),bv(b)=0;end
        lo=lo+sizes(b);
    end
    drift=abs(beta)*(240-k)/sqrt(variance);span=(max(bm)-min(bm))/sqrt(variance);
    assert(all(isfinite([beta,drift,span])),'CUMCM:q4:SelectionOverflow','Nonfinite overall statistics.');
    reason="";ratio=NaN;
    if any(sizes<2)||any(~isfinite([bm,bv]))||any(bv<=0)
        reason="DEGENERATE_BLOCK";
    else
        ratio=max(bv)/min(bv);
        if ~isfinite(ratio),reason="NONFINITE_BLOCK_RATIO";end
    end
    if drift>cfg.driftMax,reason=reason+";DRIFT";end
    if span>cfg.blockSpanMax,reason=reason+";BLOCK_MEAN";end
    if isfinite(ratio)&&ratio>cfg.varianceRatioMax,reason=reason+";BLOCK_VARIANCE";end
    reason=strip(reason,'left',';');eligible=strlength(reason)==0;
    rows{ii}=struct('signal',label,'unit',unit,'threshold',threshold,'first_exceed_min',first, ...
        'start_min',k,'end_min',240,'start_s',k*60,'end_s',14400,'duration_min',240-k, ...
        'sample_count',n,'mean_value',mu,'sample_variance',variance,'slope_per_min',beta, ...
        'block_sizes',sizes,'block_means',bm,'block_variances',bv, ...
        'standardized_drift',drift,'three_block_mean_span',span,'three_block_variance_ratio',ratio, ...
        'screen_pass',eligible,'rejection_reason',reason);
end
windows=struct2table(vertcat(rows{:}));vmin=min(windows.sample_variance);nmax=max(windows.sample_count);
windows.variance_penalty=cfg.lambda*log(windows.sample_variance/vmin);
windows.sample_penalty=(1-cfg.lambda)*log(nmax./windows.sample_count);
windows.total_loss=windows.variance_penalty+windows.sample_penalty;
assert(all(isfinite(windows.total_loss)),'CUMCM:q4:SelectionOverflow','Nonfinite loss.');
[chosen,tol]=q4_choose_suffix(windows.total_loss,windows.screen_pass,windows.sample_count,cfg.tieFactor);
assert(~isempty(chosen),'CUMCM:q4:NoStableWindow','No suffix meets the frozen stability conditions.');
windows.selected=false(count,1);windows.selected(chosen)=true;
windows.loss_rank=zeros(count,1);[~,order]=sortrows([windows.total_loss,-windows.sample_count],[1,2]);
windows.loss_rank(order)=(1:count)';
selection=table2struct(windows(chosen,:));selection.candidate_count=count;
selection.screen_pass_count=sum(windows.screen_pass);selection.variance_reference=vmin;
selection.maximum_sample_count=nmax;selection.tie_tolerance=tol;
selection.effective_loss_for_ineligible='excluded; not a third penalty';
weights=[.60,.70,.75,.80,.85,.90];weightStarts=NaN(size(weights));
for j=1:numel(weights)
    losses=weights(j)*log(windows.sample_variance/vmin)+(1-weights(j))*log(nmax./windows.sample_count);
    index=q4_choose_suffix(losses,windows.screen_pass,windows.sample_count,cfg.tieFactor);
    if ~isempty(index),weightStarts(j)=windows.start_min(index);end
end
screen=zeros(6,3);jj=0;
for d=[.25,.5,1]
    for b=[.5,1]
        jj=jj+1;ok=windows.standardized_drift<=d&windows.three_block_mean_span<=b& ...
            isfinite(windows.three_block_variance_ratio)&windows.three_block_variance_ratio<=3;
        index=q4_choose_suffix(windows.total_loss,ok,windows.sample_count,cfg.tieFactor);start=NaN;
        if ~isempty(index),start=windows.start_min(index);end
        screen(jj,:)=[d,b,start];
    end
end
result=struct('windows',windows,'selection',selection,'sensitivity', ...
    struct('weights',weights,'weight_selected_start_min',weightStarts, ...
    'screen_columns',{{'drift_max','block_span_max','selected_start_min'}}, ...
    'screen_variance_ratio_max',3,'screen_results',screen));
end


function [index,tolerance]=q4_choose_suffix(loss,eligible,n,tieFactor)
indices=find(eligible);index=[];tolerance=NaN;
if isempty(indices),return;end
assert(all(isfinite(loss(indices)))&&all(isfinite(n(indices)))&&all(n(indices)>0));
tolerance=tieFactor*eps(max(1,max(abs(loss(indices)))));
candidates=indices(loss(indices)<=min(loss(indices))+tolerance);
[~,j]=max(n(candidates));index=candidates(j);
end


function mesh=build_q4_mesh(N)
assert(N>=20&&N==fix(N)&&abs(log2(N/20)-round(log2(N/20)))<eps);
xi=(0:N)'/N; bounds=[0;(xi(1:end-1)+xi(2:end))/2;1];
mesh=struct('N',N,'xi',xi,'dr',1/N,'faces',bounds(2:end-1), ...
    'w',diff(bounds.^2)/2,'radius',.02);
assert(abs(sum(mesh.w)-.5)<32*eps);
end


function bc=interpolate_q4_ambient(t,in,scenario)
a=in.ambient; tt=min(t,14400);
assert(isscalar(t)&&isfinite(t)&&t>=0);
scenario=string(scenario);
assert(ismember(scenario,["unified_stable_suffix","temperature_last", ...
    "temperature_mean30","temperature_mean58","temperature_mean59","temperature_mean60"]), ...
    'CUMCM:q4:Ambient','Unknown or obsolete scenario; old complete boundary is not reusable.');
v=interp1(a(:,1),a(:,2:3),tt,'linear');
if t>14400
    v=[a(end,2),in.environment_selection.selection.moisture_kgkg];
    if scenario=="unified_stable_suffix"
        v(1)=in.environment_selection.selection.temperature_C;
    elseif scenario~="temperature_last"
        minutes=str2double(extractAfter(string(scenario),'temperature_mean'));
        values=a(a(:,1)>=14400-60*minutes,2);
        assert(numel(values)==minutes+1);v(1)=mean(values);
    end
end
bc=struct('temperature_K',v(1)+273.15,'moisture_kgkg',v(2),'extrapolated',t>14400);
end


function geo=interpolate_q4_radius(t,in,mode)
assert(all(t>=0&isfinite(t),'all'));
if mode=="fixed"
    R=.02+zeros(size(t)); d=zeros(size(t));
else
    clamped=min(t,259200);
    if mode=="pchip"
        R=ppval(in.radius_pp,clamped);d=ppval(in.radius_derivative_pp,clamped);
    elseif mode=="linear"
        R=.01*interp1(in.radius(:,1),in.radius(:,2),clamped,'linear');
        idx=min(floor(clamped/1800)+1,144);
        slopes=.01*diff(in.radius(:,2))/1800;d=reshape(slopes(idx),size(t));
    else
        error('CUMCM:q4:Geometry','Unknown geometry.');
    end
    d(t>=259200)=0;
end
assert(all(R>0&isfinite(R),'all'));
geo=struct('R_m',R,'Rdot_m_s',d,'radius_extrapolated',t>259200);
end


function material=q4_material_properties(moisture,temperature_K,p)
validState=isreal(moisture)&&isreal(temperature_K) && ...
    isequal(size(moisture),size(temperature_K)) && ...
    all(isfinite(moisture)&isfinite(temperature_K)&moisture>0&temperature_K>0,'all');
assert(validState, ...
    'CUMCM:q4:MaterialDomain','Material states must be positive and finite.');
material.rho=p.rho0+p.rhoC*moisture;
material.cp=p.cp0+p.cpC*moisture./(1+moisture);
material.k=p.k0+p.kC*moisture./(1+moisture);
material.b=material.rho.*material.cp;
material.D=p.D0*exp(-p.a./moisture-p.B./temperature_K);
material.alpha=material.k./material.b;
material.rho_C=p.rhoC+zeros(size(moisture));
material.cp_C=p.cpC./(1+moisture).^2;
material.k_C=p.kC./(1+moisture).^2;
material.b_C=material.rho_C.*material.cp+material.rho.*material.cp_C;
material.D_C=material.D.*p.a./moisture.^2;
material.D_T=material.D.*p.B./temperature_K.^2;
assert(all(material.D>0&isfinite(material.b)&material.b>0,'all'), ...
    'CUMCM:q4:MaterialDomain','Nonpositive or nonfinite coefficients.');
end


function face=q4_face_coefficients(CL,CR,TL,TR,p,level)
persistent cachedN cachedX cachedW
n=8*2^level;
if isempty(cachedN)||cachedN~=n
    j=(1:n-1)';b=j./sqrt(4*j.^2-1);
    [V,L]=eig(diag(b,1)+diag(b,-1));[x,ind]=sort(diag(L));
    cachedX=(x'+1)/2;cachedW=V(1,ind).^2;cachedN=n;
end
x=cachedX;w=cachedW;
C=CL+(CR-CL).*x;T=TL+(TR-TL).*x;
assert(all(C>0&T>0&isfinite(C)&isfinite(T),'all'));
k=p.k0+p.kC*C./(1+C);D=p.D0*exp(-p.a./C-p.B./T);
kC=p.kC./(1+C).^2;DC=D.*p.a./C.^2;DT=D.*p.B./T.^2;
face=struct('k',k*w','D',D*w','k_CL',(kC.*(1-x))*w', ...
    'k_CR',(kC.*x)*w','D_CL',(DC.*(1-x))*w','D_CR',(DC.*x)*w', ...
    'D_TL',(DT.*(1-x))*w','D_TR',(DT.*x)*w');
assert(all(face.k>0&face.D>0));
end


function [residual,jacobian,audit]=q4_residual_jacobian(candidate,previous,mesh,dt_s,boundary,p,quadratureLevel)
nodeCount=mesh.N+1;
systemSize=2*nodeCount;
temperatureIndex=(1:2:systemSize)';
moistureIndex=temperatureIndex+1;
assert(isequal(size(candidate),[systemSize,1]) && ...
    isequal(size(previous),[systemSize,1])&&dt_s>0);

temperature_K=candidate(temperatureIndex);
moisture=candidate(moistureIndex);
previousTemperature_K=previous(temperatureIndex);
previousMoisture=previous(moistureIndex);
material=q4_material_properties(moisture,temperature_K,p);
face=q4_face_coefficients(moisture(1:end-1),moisture(2:end), ...
    temperature_K(1:end-1),temperature_K(2:end),p,quadratureLevel);

temperatureJump_K=diff(temperature_K);
moistureJump=diff(moisture);
radialGeometry=mesh.faces/(mesh.dr*mesh.radius^2);
heatFaceFlux=-radialGeometry.*face.k.*temperatureJump_K;
moistureFaceFlux=-radialGeometry.*face.D.*moistureJump;
surfaceHeatFlux=p.hc/mesh.radius*(temperature_K(end)-boundary.temperature_K);
surfaceMoistureFlux=p.hm/mesh.radius*(moisture(end)-boundary.moisture_kgkg);
heatStorage=mesh.w.*material.b.*(temperature_K-previousTemperature_K)/dt_s;
moistureStorage=mesh.w.*(moisture-previousMoisture)/dt_s;
heatResidual=heatStorage+diff([0;heatFaceFlux;surfaceHeatFlux]);
moistureResidual=moistureStorage+diff([0;moistureFaceFlux;surfaceMoistureFlux]);
residual=zeros(systemSize,1);
residual(temperatureIndex)=heatResidual;
residual(moistureIndex)=moistureResidual;

heatFaceJacobian=[radialGeometry.*face.k, ...
    -radialGeometry.*face.k_CL.*temperatureJump_K, ...
    -radialGeometry.*face.k, ...
    -radialGeometry.*face.k_CR.*temperatureJump_K];
moistureFaceJacobian=[-radialGeometry.*face.D_TL.*moistureJump, ...
    radialGeometry.*(face.D-face.D_CL.*moistureJump), ...
    -radialGeometry.*face.D_TR.*moistureJump, ...
    -radialGeometry.*(face.D+face.D_CR.*moistureJump)];
faceColumns=[temperatureIndex(1:end-1),moistureIndex(1:end-1), ...
    temperatureIndex(2:end),moistureIndex(2:end)];
faceRows=[repmat(temperatureIndex(1:end-1),1,4); ...
    repmat(temperatureIndex(2:end),1,4); ...
    repmat(moistureIndex(1:end-1),1,4); ...
    repmat(moistureIndex(2:end),1,4)];
faceColumnPattern=repmat(faceColumns,4,1);
faceValues=[heatFaceJacobian;-heatFaceJacobian; ...
    moistureFaceJacobian;-moistureFaceJacobian];
jacobian=sparse(faceRows(:),faceColumnPattern(:),faceValues(:),systemSize,systemSize);

jacobian=jacobian+sparse( ...
    [temperatureIndex;temperatureIndex;moistureIndex; ...
    temperatureIndex(end);moistureIndex(end)], ...
    [temperatureIndex;moistureIndex;moistureIndex; ...
    temperatureIndex(end);moistureIndex(end)], ...
    [mesh.w.*material.b/dt_s; ...
    mesh.w.*material.b_C.*(temperature_K-previousTemperature_K)/dt_s; ...
    mesh.w/dt_s;p.hc/mesh.radius;p.hm/mesh.radius], ...
    systemSize,systemSize);

rowScale=zeros(systemSize,1);
rowScale(temperatureIndex)=dt_s./(mesh.w.*material.b);
rowScale(moistureIndex)=dt_s./mesh.w;
scaledResidual=residual.*rowScale;
audit=struct('row_scale',rowScale, ...
    'scaled_residual',[max(abs(scaledResidual(temperatureIndex))), ...
    max(abs(scaledResidual(moistureIndex)))], ...
    'balance',[sum(heatStorage)+surfaceHeatFlux, ...
    sum(moistureStorage)+surfaceMoistureFlux], ...
    'summed_residual',[sum(heatResidual),sum(moistureResidual)], ...
    'material',material,'heat_face_flux',heatFaceFlux, ...
    'moisture_face_flux',moistureFaceFlux);
end


function [state,audit]=advance_q4_implicit(previous,t_s,dt_s,mesh,p,in)
assert(isscalar(dt_s)&&isfinite(dt_s)&&dt_s>0&&t_s+dt_s>t_s);
assert(isequal(size(previous),[2*(mesh.N+1),1])&&all(previous>0&isfinite(previous)));
boundary=interpolate_q4_ambient(t_s+dt_s,in,p.environment);
geo=interpolate_q4_radius(t_s+dt_s,in,p.geometry); mesh.radius=geo.R_m;
state=previous; dimension=numel(state); ti=1:2:dimension; ci=2:2:dimension;
audit=struct('converged',false,'iterations',0,'backtracks',0,'correction',[Inf,Inf], ...
    'scaled_residual',[Inf,Inf],'balance',[NaN,NaN],'balance_defect',[NaN,NaN]);
cached=false;
for iteration=1:p.newtonMax
    if ~cached
        [res,J,raw]=q4_residual_jacobian(state,previous,mesh,dt_s,boundary,p,p.quadratureLevel);
    end
    scale=raw.row_scale;
    [oldMessage,oldId]=lastwarn; lastwarn('');
    update=-(spdiags(scale,0,dimension,dimension)*J)\(scale.*res);
    [message,~]=lastwarn; lastwarn(oldMessage,oldId);
    audit.iterations=iteration;
    if ~isempty(message)||any(~isfinite(update)),return;end
    audit.correction=[max(abs(update(ti))),max(abs(update(ci)))];
    audit.scaled_residual=raw.scaled_residual;
    if all(audit.correction<=p.newtonTolerance)&&all(raw.scaled_residual<=p.newtonTolerance)
        audit.converged=true; audit.balance=raw.balance;
        audit.balance_defect=abs(raw.balance-raw.summed_residual);
        return
    end
    eta=1; accepted=false;
    while eta>=p.backtrackMin
        trial=state+eta*update;
        if all(trial>0&isfinite(trial))
            [trialRes,trialJ,trialRaw]=q4_residual_jacobian(trial,previous,mesh,dt_s,boundary,p,p.quadratureLevel);
            if max(trialRaw.scaled_residual)<=(1-p.armijo*eta)*max(raw.scaled_residual)
                state=trial;res=trialRes;J=trialJ;raw=trialRaw;cached=true;accepted=true;break
            end
        end
        eta=eta/2; audit.backtracks=audit.backtracks+1;
    end
    if ~accepted,return;end
end
end


function f=reconstruct_q4_fields(y,mesh,geo,bc,p,query)
T=y(1:2:end);C=y(2:2:end);m=q4_material_properties(C,T,p);
sC=zeros(size(C));sT=zeros(size(T));h=mesh.dr;
sC(2:end-1)=(C(3:end)-C(1:end-2))/(2*h);
sT(2:end-1)=(T(3:end)-T(1:end-2))/(2*h);
sC(end)=-geo.R_m*p.hm*(C(end)-bc.moisture_kgkg)/m.D(end);
sT(end)=-geo.R_m*p.hc*(T(end)-bc.temperature_K)/m.k(end);
cC=coeff(C,sC,h);cT=coeff(T,sT,h);
f=struct('Ccoef',cC,'Tcoef',cT);
if nargin>=6&&~isempty(query)
    assert(all(query>=0&query<=1));
    ix=min(floor(query(:)/h)+1,mesh.N);x=(query(:)-mesh.xi(ix))/h;
    f.C=((cC(ix,1).*x+cC(ix,2)).*x+cC(ix,3)).*x+cC(ix,4);
    f.T=((cT(ix,1).*x+cT(ix,2)).*x+cT(ix,3)).*x+cT(ix,4);
end
end
function a=coeff(u,s,h)
a=[2*u(1:end-1)-2*u(2:end)+h*(s(1:end-1)+s(2:end)), ...
    -3*u(1:end-1)+3*u(2:end)-h*(2*s(1:end-1)+s(2:end)), ...
    h*s(1:end-1),u(1:end-1)];
end


function m=q4_spatial_maximum(y,mesh,geo,bc,p)
f=reconstruct_q4_fields(y,mesh,geo,bc,p,[]);
a=f.Ccoef;values=y(2:2:end);locations=mesh.xi;
for i=1:mesh.N
    b=[3*a(i,1),2*a(i,2),a(i,3)];start=find(b~=0,1);
    if isempty(start),continue;end
    rootsHere=roots(b(start:end)/max(abs(b)));
    x=real(rootsHere(abs(imag(rootsHere))<=64*eps(max(1,abs(real(rootsHere))))));
    x=x(x>0&x<1);
    values=[values;polyval(a(i,:),x)];locations=[locations;mesh.xi(i)+mesh.dr*x]; %#ok<AGROW>
end
peak=max(values);ix=abs(values-peak)<=64*eps(max(1,abs(peak)));
m=struct('maximum_C',peak,'minimum_C',min(values),'node_max',max(y(2:2:end)), ...
    'argmax_xi',unique(locations(ix)),'coefficients',a);
end


function d=q4_progress_metrics(y,mesh,geo,bc,p,withFraction)
if nargin<6,withFraction=true;end
f=reconstruct_q4_fields(y,mesh,geo,bc,p,[]);h=mesh.dr;a=f.Ccoef;
meanC=2*h*sum(mesh.xi(1:end-1).*(a*[1/4;1/3;1/2;1])+h*(a*[1/5;1/4;1/3;1/2]));
fraction=0;
if withFraction
    for i=1:mesh.N
        poly=a(i,:);poly(end)=poly(end)-p.threshold;
        first=find(poly~=0,1);
        x=[];
        if ~isempty(first)
            rr=roots(poly(first:end)/max(abs(poly)));
            x=real(rr(abs(imag(rr))<=64*eps(max(1,abs(real(rr))))));x=x(x>0&x<1);
        end
        edges=unique([0;x(:);1]);mid=(edges(1:end-1)+edges(2:end))/2;
        dry=polyval(a(i,:),mid)<p.threshold;
        left=mesh.xi(i)+h*edges(1:end-1);right=mesh.xi(i)+h*edges(2:end);
        fraction=fraction+sum((right(dry).^2-left(dry).^2));
    end
else
    fraction=NaN;
end
f=reconstruct_q4_fields(y,mesh,geo,bc,p,[0;.5;1]);
c=p.a*(1/p.C0-1./f.C);v=p.B*(1/p.T0_K-1./f.T);
g=2*log(p.R0_m/geo.R_m)+zeros(3,1);
logNow=log(p.D0)-p.a./f.C-p.B./f.T-2*log(geo.R_m);
logInitial=log(p.D0)-p.a/p.C0-p.B/p.T0_K-2*log(p.R0_m);
direct=logNow-logInitial;
assert(max(abs(c+v+g-direct))<1e-10);
d=struct('mean_C',meanC,'dry_fraction',fraction,'moisture_log',c, ...
    'temperature_log',v,'geometry_log',g,'total_log',direct, ...
    'C_at_material_points',f.C,'T_at_material_points',f.T, ...
    'Lambda',exp(logNow),'closure_residual',c+v+g-direct);
end


function sol=solve_q4_longtime(initial,mesh,in,p,context)
t=initial.time_s;y=initial.state(:);finish=context.stop_time_s;timer=tic;
if isfield(context,'total_timer'),timer=context.total_timer;end
stream=isfield(context,'trace_path')&&strlength(string(context.trace_path))>0;
trace=[];dim=numel(y);ti=1:2:dim;ci=2:2:dim;
geo=interpolate_q4_radius(t,in,p.geometry);bc=interpolate_q4_ambient(t,in,p.environment);
if t==0
    maximum=max(y(ci));meanC=2*sum(mesh.w.*y(ci));
else
    m=q4_spatial_maximum(y,mesh,geo,bc,p);maximum=m.maximum_C;
    prog=q4_progress_metrics(y,mesh,geo,bc,p,false);meanC=prog.mean_C;
end
if stream
    trace=struct('path',context.trace_path,'maxBytes',16*1024^3);
    trace=q4_trace_store('init',trace,[t;y;0;geo.R_m;bc.temperature_K;bc.moisture_kgkg;maximum;meanC]);
end
times=((floor(t/60)+1)*60:60:finish)';
if isempty(times)||times(end)~=finish,times(end+1)=finish;end
states=zeros(dim,numel(times));out=1;accepted=0;rejected=0;maxNewton=0;
proposed=min(p.initialStep_s,p.maxStep_s);bracket=[];avgBracket=[];
errMax=[0,0];overshoot=0;lastPrint=tic;
try
    while t<finish
        assert(toc(timer)<p.wallSeconds,'CUMCM:q4:ResourceLimit','Wall budget exceeded; checkpoint retained.');
        dt=min([proposed,p.maxStep_s,times(out)-t,finish-t]);
        if t<14400,dt=min(dt,60*(floor(t/60)+1)-t);end
        if t<259200,dt=min(dt,1800*(floor(t/1800)+1)-t);end
        assert(dt>0&&t+dt/2>t,'CUMCM:q4:StepResolution','Cannot represent smaller time step.');
        ht=[t+dt/2,t+dt];hd=[ht(1)-t,ht(2)-ht(1)];
        [full,af]=advance_q4_implicit(y,t,dt,mesh,p,in);
        [half,ah]=advance_q4_implicit(y,t,hd(1),mesh,p,in);
        if ~af.converged||~ah.converged,rejected=rejected+1;proposed=dt/2;continue;end
        [two,at]=advance_q4_implicit(half,ht(1),hd(2),mesh,p,in);
        if ~at.converged,rejected=rejected+1;proposed=dt/2;continue;end
        err=[max(abs(two(ti)-full(ti))),max(abs(two(ci)-full(ci)))];
        ratio=max(err./p.localTolerance);
        if ratio>1,rejected=rejected+1;proposed=dt*max(.2,.9/sqrt(ratio));continue;end
        for k=1:2
            oldt=t;oldy=y;oldmax=maximum;oldmean=meanC;
            if k==1,y=half;aud=ah;else,y=two;aud=at;end
            t=ht(k);geo=interpolate_q4_radius(t,in,p.geometry);bc=interpolate_q4_ambient(t,in,p.environment);
            assert(all(y(ti)>=301.15-1e-8&y(ti)<=323.396+1e-8)&& ...
                all(y(ci)>=.01963-1e-9&y(ci)<=p.C0+1e-9),'CUMCM:q4:Range','State outside input range.');
            m=q4_spatial_maximum(y,mesh,geo,bc,p);maximum=m.maximum_C;
            prog=q4_progress_metrics(y,mesh,geo,bc,p,false);meanC=prog.mean_C;
            overshoot=max(overshoot,maximum-m.node_max);
            if isempty(avgBracket)&&oldmean>=p.threshold&&meanC<p.threshold
                avgBracket=struct('left_s',oldt,'right_s',t,'left_state',oldy,'right_state',y);
            end
            if isempty(bracket)&&oldmax>=p.threshold&&maximum<p.threshold
                bracket=struct('left_s',oldt,'right_s',t,'left_state',oldy,'right_state',y);
            end
            accepted=accepted+1;maxNewton=max(maxNewton,aud.iterations);
            if stream,trace=q4_trace_store('append',trace,[t;y;hd(k);geo.R_m;bc.temperature_K;bc.moisture_kgkg;maximum;meanC]);end
            if p.eventStop&&~isempty(bracket),break;end
        end
        errMax=max(errMax,err);
        if t==times(out),states(:,out)=y;out=out+1;end
        if p.eventStop&&~isempty(bracket),break;end
        proposed=dt*min(2,max(.2,.9/sqrt(max(ratio,1e-12))));
        if toc(lastPrint)>30,fprintf('q4 N%d: %.3f h, max C=%.8g, steps=%d, elapsed=%.1fs\n',mesh.N,t/3600,maximum,accepted,toc(timer));lastPrint=tic;end
    end
catch exception
    if stream
        q4_trace_store('checkpoint',trace,[]);
        checkpoint=struct('time_s',t,'state',y,'parameters',p,'mesh',mesh, ...
            'accepted',accepted,'rejected',rejected,'elapsed_s',toc(timer),'error',exception.identifier);
        save(fullfile(fileparts(context.trace_path),'checkpoint.mat'),'checkpoint','-v7');
    end
    rethrow(exception)
end
if stream,trace=q4_trace_store('finalize',trace,[]);trace=rmfield(trace,'buffer');end
sol=struct('time_s',times(1:out-1),'states',states(:,1:out-1),'final_time_s',t, ...
    'final_state',y,'event_bracket',bracket,'average_bracket',avgBracket,'trace',trace, ...
    'accepted_half_steps',accepted,'rejected_macro_steps',rejected,'max_newton_iterations',maxNewton, ...
    'local_error_max',errMax,'reconstruction_overshoot_max',overshoot,'runtime_seconds',toc(timer));
end


function event=locate_q4_event(bracket,mesh,in,p,context,kind)
if nargin<6,kind="maximum";end
initial=struct('time_s',bracket.left_s,'state',bracket.left_state);
left=bracket.left_s;right=bracket.right_s;p.eventStop=false;
vl=evaluate(left);vr=evaluate(right);
assert(vl.value>=p.threshold,'CUMCM:q4:Bracket','Left state is already below threshold.');
extension=right-left;extensions=0;
while vr.value>=p.threshold
    right=right+extension;extension=2*extension;extensions=extensions+1;vr=evaluate(right);
    assert(isfinite(right),'CUMCM:q4:Bracket','Bracket extension overflow.');
end
while right-left>p.rootWidth_s
    mid=left+(right-left)/2;v=evaluate(mid);
    if v.value<p.threshold,right=mid;else,left=mid;end
end
stop=round(p.stopLattice_s*ceil(right/p.stopLattice_s)*1000)/1000;
if stop<right,stop=stop+p.stopLattice_s;end
v=evaluate(stop);
while v.value>=p.threshold,stop=stop+p.stopLattice_s;v=evaluate(stop);end
rootTime=(left+right)/2;rootValue=evaluate(rootTime);
event=struct('kind',kind,'root_s',rootTime,'left_s',left,'right_s',right, ...
    'stop_time_s',stop,'stop_time_h',stop/3600,'value',v.value,'stop_state',v.state, ...
    'root_state',rootValue.state,'strict_nominal_stop',v.value<p.threshold,'numerical_precision_verified',false);
event.reintegration_bracket_extensions=extensions;
    function v=evaluate(t)
        if t==initial.time_s,state=initial.state;
        else
            local=context;local.stop_time_s=t;local.trace_path="";
            s=solve_q4_longtime(initial,mesh,in,p,local);state=s.final_state;
        end
        geo=interpolate_q4_radius(t,in,p.geometry);bc=interpolate_q4_ambient(t,in,p.environment);
        if kind=="maximum",m=q4_spatial_maximum(state,mesh,geo,bc,p);val=m.maximum_C;
        else,m=q4_progress_metrics(state,mesh,geo,bc,p,false);val=m.mean_C;end
        v=struct('state',state,'value',val);
    end
end


function fields=extract_q4_output_grid(sol,event,in,p,mesh)
times=sol.time_s;states=sol.states;
if ~isempty(event)
    keep=times<event.stop_time_s;times=[times(keep);event.stop_time_s];states=[states(:,keep),event.stop_state];
end
assert(~isempty(times)&&all(diff(times)>0)&&all(times>0));
g=interpolate_q4_radius(times,in,p.geometry);radii=100*g.R_m;
distance=(0:floor(max(radii)/.1))*.1;
mask=distance<=radii+64*eps(max(1,radii));
C=NaN(numel(times),numel(distance)+1);T=C;progress=zeros(numel(times),3);
mechanism=zeros(numel(times)*3,11);
for j=1:numel(times)
    geo=interpolate_q4_radius(times(j),in,p.geometry);bc=interpolate_q4_ambient(times(j),in,p.environment);
    query=[distance(mask(j,:))/radii(j),1]';
    query(abs(query-1)<=64*eps)=1;
    f=reconstruct_q4_fields(states(:,j),mesh,geo,bc,p,query);
    C(j,[find(mask(j,:)),numel(distance)+1])=f.C;
    T(j,[find(mask(j,:)),numel(distance)+1])=f.T-273.15;
    m=q4_spatial_maximum(states(:,j),mesh,geo,bc,p);
    d=q4_progress_metrics(states(:,j),mesh,geo,bc,p,true);
    progress(j,:)=[m.maximum_C,d.mean_C,d.dry_fraction];
    ix=(j-1)*3+(1:3);
    mechanism(ix,:)=[repmat(times(j),3,1),[0;.5;1],repmat(radii(j),3,1), ...
        d.C_at_material_points,d.T_at_material_points,d.Lambda,d.total_log, ...
        d.moisture_log,d.temperature_log,d.geometry_log,d.closure_residual];
end
fields=struct('time_s',times,'radius_cm',radii,'distance_cm',distance, ...
    'mask',mask,'C',C,'temperature_C',T,'progress',progress,'mechanism',mechanism);
end


function candidate=q4_solve_case(cfg,in,p,caseRoot,totalTimer)
mesh=build_q4_mesh(p.N);p.eventStop=true;
initial=struct('time_s',0,'state',repmat([p.T0_K;p.C0],p.N+1,1));
times=[];states=zeros(2*(p.N+1),0);traces=strings(0,1);average=[];window=0;
accepted=0;rejected=0;meanEvent=[];
while true
    if window==0,finish=259200;else,finish=259200+21600*2^(window-1);end
    relative=string(caseRoot)+"/window_"+window+"/accepted_states.h5";
    context=struct('stop_time_s',finish,'trace_path',fullfile(cfg.projectRoot,relative),'total_timer',totalTimer);
    s=solve_q4_longtime(initial,mesh,in,p,context);
    times=[times;s.time_s];states=[states,s.states];traces(end+1)=relative; %#ok<AGROW>
    accepted=accepted+s.accepted_half_steps;rejected=rejected+s.rejected_macro_steps;
    if isempty(average)&&~isempty(s.average_bracket)
        average=s.average_bracket;
        meanEvent=locate_q4_event(average,mesh,in,p,context,"mean");
    end
    if ~isempty(s.event_bracket),break;end
    initial=struct('time_s',s.final_time_s,'state',s.final_state);window=window+1;
    fprintf('q4 %s: not dry by %.1f h; declared radius continuation applied.\n',caseRoot,initial.time_s/3600);
end
event=locate_q4_event(s.event_bracket,mesh,in,p,context,"maximum");
s.time_s=times;s.states=states;
fields=extract_q4_output_grid(s,event,in,p,mesh);
if ~isempty(meanEvent)
    meanEvent.root_gap_s=event.root_s-meanEvent.root_s;
    meanEvent.gap_interval_s=[event.left_s-meanEvent.right_s,event.right_s-meanEvent.left_s];
    meanEvent.first_crossing_numerically_certified=false;
    meanEvent.first_crossing_note="Accepted half-step scan; global time refinement and near-threshold history still required.";
end
candidate=struct('parameters',rmfield(p,'projectRoot'),'mesh',mesh,'fields',fields, ...
    'event',event,'mean_event',meanEvent,'trace_paths',traces,'states',states, ...
    'time_s',times,'accepted_half_steps',accepted,'rejected_macro_steps',rejected, ...
    'runtime_seconds',toc(totalTimer),'numerically_validated',false);
candidate.environment_selection=in.environment_selection.selection;
end
