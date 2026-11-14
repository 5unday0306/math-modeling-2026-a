function varargout=q3_compute(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q3:ComputeAction','Invalid compute action.');
switch action
    case "config"
        varargout{1}=q3_config(varargin{:});
    case "parameters"
        varargout{1}=q3_parameters(varargin{:});
    case "q2_parameters"
        varargout{1}=q2_parameters();
    case "build_mesh"
        varargout{1}=build_q2_mesh(varargin{:});
    case "q2_material_properties"
        varargout{1}=q2_material_properties(varargin{:});
    case "q2_face_coefficients"
        varargout{1}=q2_face_coefficients(varargin{:});
    case "q2_residual_jacobian"
        [varargout{1:nargout}]=q2_residual_jacobian(varargin{:});
    case "interpolate_ambient"
        varargout{1}=interpolate_q3_ambient(varargin{:});
    case "advance_step"
        [varargout{1:nargout}]=advance_q3_implicit(varargin{:});
    case "solve_longtime"
        varargout{1}=solve_q3_longtime(varargin{:});
    case "solve_event_case"
        varargout{1}=q3_solve_event_case(varargin{:});
    case "event_evaluate"
        varargout{1}=q3_event_evaluate(varargin{:});
    case "locate_event"
        varargout{1}=locate_q3_event(varargin{:});
    case "spatial_maximum"
        varargout{1}=q3_spatial_maximum(varargin{:});
    case "output_times"
        varargout{1}=q3_output_times(varargin{:});
    case "extract_fields"
        varargout{1}=extract_q3_fields(varargin{:});
    case "select_window"
        varargout{1}=select_equilibrium_window(varargin{:});
    case "tail_maximum_history"
        varargout{1}=tail_maximum_history(varargin{:});
    otherwise
        error('CUMCM:q3:ComputeAction','Unknown compute action: %s',action);
end
end

function cfg=q3_config(root)
cfg=struct('projectRoot',string(root),'rawDataDir',fullfile(root,'data','raw'), ...
    'resultsDir',fullfile(root,'附件'),'processedDataDir',fullfile(root,'附件'), ...
    'figuresDir',fullfile(root,'图片源文件'),'figureQADir',fullfile(root,'计算缓存','图片预览'));
end

function p=q3_parameters(cfg,upstream)
p=q2_parameters();
p.projectRoot=cfg.projectRoot; p.threshold=0.15; p.boundaryBreak_s=14400;
p.restartTime_s=10800; p.mvpEnd_s=14460; p.outputInterval_s=60;
p.paperInterval_s=21600; p.stopLattice_s=0.36; p.rootWidth_s=0.03;
p.jointEventBudget_s=0.18; p.fieldBudget=5e-5;
p.newtonTolerance=5e-9; p.newtonTightTolerance=5e-11;
p.localTolerance=[5e-6,5e-7]; p.quadratureLevel=0;
p.wallSeconds=1200; p.validationWallSeconds=3600;
p.initialStep_s=1;
if ~isempty(upstream)
    assert(isfield(upstream,'last_step_s')&&upstream.last_step_s>0);
    p.initialStep_s=upstream.last_step_s;
end
p.maxStep_s=60; p.eventStop=false;
end

function p=q2_parameters()
p.radius_m=0.02; p.length_m=0.25; p.T0_K=301.15; p.C0=2.55;
p.rho0=650; p.rhoC=128; p.cp0=1450; p.cpC=2736;
p.k0=0.21; p.kC=0.38; p.D0=0.0024; p.a=0.45; p.B=3850;
p.hc=25; p.hm=8e-7; p.tFinal=10800;
p.newtonMax=30; p.backtrackMin=2^-20; p.armijo=1e-4;
p.algebraicTolerance=[5e-7,5e-9];
p.quadRelTol=[1e-10,1e-12]; p.quadAbsTol=[1e-12,1e-14];
p.minStep=2^-20; p.maxWorkingBytes=2*1024^3;
p.maxTraceBytes=16*1024^3; p.wallSeconds=600;
end

function mesh=build_q2_mesh(radius_m,radialSegments)
assert(isscalar(radius_m)&&isfinite(radius_m)&&radius_m>0 && ...
    isscalar(radialSegments)&&radialSegments>=20&&radialSegments<=flintmax && ...
    mod(log2(radialSegments/20),1)==0,'CUMCM:q2:Mesh','Invalid grid.');
dr=radius_m/radialSegments; r=(0:radialSegments)'*dr;
faces=((0:radialSegments-1)'+0.5)*dr; edges=[0;faces;radius_m];
w=diff(edges.^2)/2;
assert(all(w>0)&&abs(sum(w)-radius_m^2/2)<=64*eps(radius_m^2) && ...
    abs(w(1)-dr^2/8)<=64*eps(dr^2),'CUMCM:q2:Mesh','Volume invariant failed.');
mesh=struct('r',r,'faces',faces,'boundaries',edges,'w',w,'dr',dr, ...
    'radius',radius_m,'N',radialSegments,'output_columns',1+(0:20)*(radialSegments/20));
end

function material=q2_material_properties(C,T_K,p)
assert(isreal(C)&&isreal(T_K)&&isequal(size(C),size(T_K))&&all(isfinite(C)&isfinite(T_K)&C>0&T_K>0,'all'), ...
    'CUMCM:q2:MaterialDomain','Material states must be positive and finite.');
material.rho=p.rho0+p.rhoC*C;
material.cp=p.cp0+p.cpC*C./(1+C); material.k=p.k0+p.kC*C./(1+C);
material.b=material.rho.*material.cp;
material.D=p.D0*exp(-p.a./C-p.B./T_K); material.alpha=material.k./material.b;
material.rho_C=p.rhoC+zeros(size(C)); material.cp_C=p.cpC./(1+C).^2;
material.k_C=p.kC./(1+C).^2;
material.b_C=material.rho_C.*material.cp+material.rho.*material.cp_C;
material.D_C=material.D.*p.a./C.^2; material.D_T=material.D.*p.B./T_K.^2;
assert(all(material.D>0&isfinite(material.b)&material.b>0,'all'), ...
    'CUMCM:q2:MaterialDomain','Nonpositive or nonfinite coefficients.');
end

function face=q2_face_coefficients(Cleft,Cright,Tleft_K,Tright_K,p,quadratureLevel)
assert(ismember(quadratureLevel,[0,1])&&isequal(size(Cleft),size(Cright),size(Tleft_K),size(Tright_K)));
assert(all(isfinite([Cleft;Cright;Tleft_K;Tright_K])&[Cleft;Cright;Tleft_K;Tright_K]>0,'all'), ...
    'CUMCM:q2:MaterialDomain','Face state outside positive domain.');
cmax=max(Cleft,Cright); cmin=min(Cleft,Cright);
tmax=max(Tleft_K,Tright_K); tmin=min(Tleft_K,Tright_K);
bound=q2_material_properties(cmax,tmax,p);
ks=bound.k; ds=bound.D; kcs=p.kC./(1+cmin).^2;
dcs=ds.*p.a./cmin.^2; dts=ds.*p.B./tmin.^2;
dc=Cright-Cleft; dt=Tright_K-Tleft_K;
[oldMessage,oldId]=lastwarn; cleanup=onCleanup(@()lastwarn(oldMessage,oldId));
lastwarn('');
integrals=integral(@integrand,0,1,'ArrayValued',true, ...
    'RelTol',p.quadRelTol(quadratureLevel+1),'AbsTol',p.quadAbsTol(quadratureLevel+1));
[message,warningId]=lastwarn;
assert(isempty(message),'CUMCM:q2:QuadratureFailure','Face quadrature warning %s: %s',warningId,message);
face=struct('k',integrals(:,1).*ks,'D',integrals(:,2).*ds, ...
    'k_CL',integrals(:,3).*kcs,'k_CR',integrals(:,4).*kcs, ...
    'D_CL',integrals(:,5).*dcs,'D_CR',integrals(:,6).*dcs, ...
    'D_TL',integrals(:,7).*dts,'D_TR',integrals(:,8).*dts);
assert(all(face.k>0&face.D>0),'CUMCM:q2:QuadratureFailure','Nonpositive face mean.');
    function values=integrand(s)
        C=Cleft+s*dc; T=Tleft_K+s*dt;
        conductivity=p.k0+p.kC*C./(1+C);
        diffusion=p.D0*exp(-p.a./C-p.B./T);
        derivK=p.kC./(1+C).^2; derivC=diffusion.*p.a./C.^2;
        derivT=diffusion.*p.B./T.^2;
        values=[conductivity./ks,diffusion./ds, ...
            derivK.*(1-s)./kcs,derivK.*s./kcs, ...
            derivC.*(1-s)./dcs,derivC.*s./dcs, ...
            derivT.*(1-s)./dts,derivT.*s./dts];
    end
end

function [residual,jacobian,audit]=q2_residual_jacobian(candidate,previous,mesh,dt_s,boundary,p,quadratureLevel)
m=mesh.N+1; d=2*m; ti=(1:2:d)'; ci=ti+1;
assert(isequal(size(candidate),[d,1])&&isequal(size(previous),[d,1])&&dt_s>0);
T=candidate(ti); C=candidate(ci); oldT=previous(ti); oldC=previous(ci);
material=q2_material_properties(C,T,p);
face=q2_face_coefficients(C(1:end-1),C(2:end),T(1:end-1),T(2:end),p,quadratureLevel);
deltaT=diff(T); deltaC=diff(C); geometry=mesh.faces/mesh.dr;
FT=-geometry.*face.k.*deltaT; FC=-geometry.*face.D.*deltaC;
surfaceT=mesh.radius*p.hc*(T(end)-boundary.temperature_K);
surfaceC=mesh.radius*p.hm*(C(end)-boundary.moisture_kgkg);
storageT=mesh.w.*material.b.*(T-oldT)/dt_s;
storageC=mesh.w.*(C-oldC)/dt_s;
RT=storageT+diff([0;FT;surfaceT]); RC=storageC+diff([0;FC;surfaceC]);
residual=zeros(d,1); residual(ti)=RT; residual(ci)=RC;
gT=[geometry.*face.k,-geometry.*face.k_CL.*deltaT, ...
    -geometry.*face.k,-geometry.*face.k_CR.*deltaT];
gC=[-geometry.*face.D_TL.*deltaC,geometry.*(face.D-face.D_CL.*deltaC), ...
    -geometry.*face.D_TR.*deltaC,-geometry.*(face.D+face.D_CR.*deltaC)];
columns=[ti(1:end-1),ci(1:end-1),ti(2:end),ci(2:end)];
rows=[repmat(ti(1:end-1),1,4);repmat(ti(2:end),1,4); ...
    repmat(ci(1:end-1),1,4);repmat(ci(2:end),1,4)];
cols=repmat(columns,4,1); values=[gT;-gT;gC;-gC];
jacobian=sparse(rows(:),cols(:),values(:),d,d);
jacobian=jacobian+sparse([ti;ti;ci;ti(end);ci(end)], ...
    [ti;ci;ci;ti(end);ci(end)], ...
    [mesh.w.*material.b/dt_s;mesh.w.*material.b_C.*(T-oldT)/dt_s; ...
    mesh.w/dt_s;mesh.radius*p.hc;mesh.radius*p.hm],d,d);
scale=zeros(d,1); scale(ti)=dt_s./(mesh.w.*material.b); scale(ci)=dt_s./mesh.w;
scaled=residual.*scale;
audit=struct('row_scale',scale,'scaled_residual',[max(abs(scaled(ti))),max(abs(scaled(ci)))], ...
    'balance',[sum(storageT)+surfaceT,sum(storageC)+surfaceC], ...
    'summed_residual',[sum(RT),sum(RC)],'material',material, ...
    'heat_face_flux',FT,'moisture_face_flux',FC);
end

function boundary=interpolate_q3_ambient(t_s,ambient,scenario)
assert(isreal(t_s)&&all(isfinite(t_s)&t_s>=0,'all'));
T=interp1(ambient.time_s,ambient.temperature_C,min(t_s,14400),'linear');
C=interp1(ambient.time_s,ambient.moisture_kgkg,min(t_s,14400),'linear');
T(t_s>14400)=scenario.temperature_C; C(t_s>14400)=scenario.moisture_kgkg;
boundary=struct('temperature_K',T+273.15,'moisture_kgkg',C);
end

function [state,audit]=advance_q3_implicit(previous,t_s,dt_s,mesh,p,ambient,scenario)
assert(isscalar(dt_s)&&isfinite(dt_s)&&dt_s>0&&t_s+dt_s>t_s);
assert(isequal(size(previous),[2*(mesh.N+1),1])&&all(previous>0&isfinite(previous)));
boundary=interpolate_q3_ambient(t_s+dt_s,ambient,scenario);
state=previous; dimension=numel(state); ti=1:2:dimension; ci=2:2:dimension;
audit=struct('converged',false,'iterations',0,'backtracks',0,'correction',[Inf,Inf], ...
    'scaled_residual',[Inf,Inf],'balance',[NaN,NaN],'balance_defect',[NaN,NaN]);
cached=false;
for iteration=1:p.newtonMax
    if ~cached
        [res,J,raw]=q2_residual_jacobian(state,previous,mesh,dt_s,boundary,p,p.quadratureLevel);
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
            [trialRes,trialJ,trialRaw]=q2_residual_jacobian(trial,previous,mesh,dt_s,boundary,p,p.quadratureLevel);
            if max(trialRaw.scaled_residual)<=(1-p.armijo*eta)*max(raw.scaled_residual)
                state=trial;res=trialRes;J=trialJ;raw=trialRaw;cached=true;accepted=true;break
            end
        end
        eta=eta/2; audit.backtracks=audit.backtracks+1;
    end
    if ~accepted,return;end
end
end

function solution=solve_q3_longtime(initial,history,mesh,p,ambient,scenario,runContext)
t=initial.time_s;y=initial.state(:);finish=runContext.stop_time_s;
assert(isfinite(finish)&&finish>t&&numel(y)==2*(mesh.N+1));
timer=tic;if isfield(runContext,'total_timer'),timer=runContext.total_timer;end
streaming=isfield(runContext,'trace_path')&&strlength(string(runContext.trace_path))>0;
trace=[];memoryTrace=[]; dimension=numel(y);ti=1:2:dimension;ci=2:2:dimension;
boundary=interpolate_q3_ambient(t,ambient,scenario);
maximum=q3_spatial_maximum(y(ci),y(ti),mesh,boundary,p);
payload=[t;y;0;0;0;0;maximum.reconstructed_max];
if streaming
    trace=struct('path',runContext.trace_path,'maxBytes',p.maxTraceBytes);
    trace=q3_trace_store('init',trace,payload);
else
    memoryTrace=payload;
end
outputTimes=((floor(t/60)+1)*60:60:finish)';
if isempty(outputTimes)||outputTimes(end)~=finish,outputTimes(end+1,1)=finish;end
assert(2*8*(numel(outputTimes)*42+128*(dimension+6)+100*dimension)<=p.maxWorkingBytes, ...
    'CUMCM:q3:ResourceLimit','Output/work arrays exceed budget.');
outT=zeros(numel(outputTimes),21);outC=outT;outMax=zeros(numel(outputTimes),1);outIndex=1;
proposed=min(p.initialStep_s,p.maxStep_s);accepted=0;rejected=0;maxNewton=0;
localErrorMax=[0,0];dtMin=Inf;dtMax=0;bracket=struct([]);
Tbound=[min([p.T0_K;initial.state(ti);ambient.temperature_C+273.15;scenario.temperature_C+273.15]), ...
    max([p.T0_K;initial.state(ti);ambient.temperature_C+273.15;scenario.temperature_C+273.15])];
Cbound=[min([p.C0;initial.state(ci);ambient.moisture_kgkg;scenario.moisture_kgkg]), ...
    max([p.C0;initial.state(ci);ambient.moisture_kgkg;scenario.moisture_kgkg])];
try
    while t<finish
        assert(toc(timer)<=p.wallSeconds,'CUMCM:q3:ResourceLimit','Run wall-clock budget exceeded.');
        dt=min([proposed,p.maxStep_s,outputTimes(outIndex)-t,finish-t]);
        if t<14400,dt=min(dt,14400-t);end
        assert(dt/2>=p.minStep&&t+dt/2>t,'CUMCM:q3:StepResolution','Cannot halve step further.');
        halfTimes=[t+dt/2,t+dt];halfDurations=[halfTimes(1)-t,halfTimes(2)-halfTimes(1)];
        [full,af]=advance_q3_implicit(y,t,dt,mesh,p,ambient,scenario);
        [half,ah]=advance_q3_implicit(y,t,halfDurations(1),mesh,p,ambient,scenario);
        if ~af.converged||~ah.converged
            rejected=rejected+1;proposed=dt/2;continue
        end
        [two,at]=advance_q3_implicit(half,halfTimes(1),halfDurations(2),mesh,p,ambient,scenario);
        if ~at.converged,rejected=rejected+1;proposed=dt/2;continue;end
        err=[max(abs(two(ti)-full(ti))),max(abs(two(ci)-full(ci)))];
        if any(err>p.localTolerance),rejected=rejected+1;proposed=dt/2;continue;end
        for halfIndex=1:2
            if halfIndex==1,next=half;a=ah;else,next=two;a=at;end
            previousTime=t;previousState=y;previousMaximum=maximum.reconstructed_max;
            t=halfTimes(halfIndex);y=next;
            assert(all(y(ti)>=Tbound(1)-64*eps(Tbound(1))&y(ti)<=Tbound(2)+64*eps(Tbound(2)))&& ...
                all(y(ci)>=Cbound(1)-64*eps(Cbound(1))&y(ci)<=Cbound(2)+64*eps(Cbound(2))), ...
                'CUMCM:q3:Range','Accepted state outside physical input bounds; not clipped.');
            boundary=interpolate_q3_ambient(t,ambient,scenario);
            maximum=q3_spatial_maximum(y(ci),y(ti),mesh,boundary,p);
            if isempty(bracket)&&previousMaximum>=p.threshold&&maximum.reconstructed_max<p.threshold
                bracket=struct('left_s',previousTime,'right_s',t,'left_state',previousState, ...
                    'right_state',y,'status','NOMINAL_NOT_VALIDATED');
            end
            payload=[t;y;halfDurations(halfIndex);a.iterations;a.balance(:);maximum.reconstructed_max];
            if streaming
                trace=q3_trace_store('append',trace,payload);
            else
                assert(2*8*numel(payload)*(size(memoryTrace,2)+1)<p.maxWorkingBytes, ...
                    'CUMCM:q3:ResourceLimit','In-memory trace exceeds budget.');
                memoryTrace(:,end+1)=payload; %#ok<AGROW>
            end
            accepted=accepted+1;maxNewton=max(maxNewton,a.iterations);
        end
        localErrorMax=max(localErrorMax,err);dtMin=min([dtMin,halfDurations]);dtMax=max([dtMax,halfDurations]);
        if abs(t-outputTimes(outIndex))<=8*eps(max(1,t))
            assert(t==outputTimes(outIndex),'CUMCM:q3:OutputTime','Integrator did not land on output time.');
            outT(outIndex,:)=y(ti(mesh.output_columns));outC(outIndex,:)=y(ci(mesh.output_columns));
            outMax(outIndex)=maximum.reconstructed_max;outIndex=outIndex+1;
        end
        if p.eventStop&&~isempty(bracket),break;end
        if all(err<=p.localTolerance/4),proposed=2*dt;else,proposed=dt;end
    end
catch exception
    if streaming,trace=q3_trace_store('checkpoint',trace,[]);end %#ok<NASGU>
    rethrow(exception)
end
if streaming,trace=q3_trace_store('finalize',trace,[]);else,trace=struct('memory',memoryTrace);end
used=outIndex-1;
solution=struct('time_s',outputTimes(1:used),'temperature_K',outT(1:used,:), ...
    'moisture_kgkg',outC(1:used,:),'maximum_C',outMax(1:used), ...
    'final_state',y,'final_time_s',t,'radius_m',mesh.r,'trace',trace,'history',history, ...
    'event_bracket',bracket,'accepted_half_steps',accepted,'rejected_macro_steps',rejected, ...
    'newton_max_iterations',maxNewton,'local_error_max',localErrorMax, ...
    'dt_min_s',dtMin,'dt_max_s',dtMax,'runtime_seconds',toc(timer));
end

function value=q3_event_evaluate(t,initial,mesh,p,ambient,scenario,totalTimer)
assert(t>=initial.time_s);p.eventStop=false;
if t==initial.time_s
    state=initial.state;
else
    context=struct('stop_time_s',t,'trace_path','');
    if nargin>=7,context.total_timer=totalTimer;end
    solution=solve_q3_longtime(initial,[],mesh,p,ambient,scenario,context);state=solution.final_state;
end
b=interpolate_q3_ambient(t,ambient,scenario);
m=q3_spatial_maximum(state(2:2:end),state(1:2:end),mesh,b,p);
value=struct('state',state,'maximum_C',m.reconstructed_max,'argmax_radius_m',m.argmax_radius_m);
end

function event=locate_q3_event(solution,evaluator,errorIndicators,p)
assert(~isempty(solution.event_bracket),'CUMCM:q3:NoBracket','No nominal first-event bracket.');
left=solution.event_bracket.left_s;right=solution.event_bracket.right_s;
assert(right>left);leftValue=evaluator(left);rightValue=evaluator(right);
assert(leftValue.maximum_C>=p.threshold&&rightValue.maximum_C<p.threshold);
timer=tic;ambiguous=false;
while right-left>p.rootWidth_s
    assert(toc(timer)<=p.wallSeconds,'CUMCM:q3:ResourceLimit','Event evaluation budget exceeded.');
    mid=left+(right-left)/2;
    assert(mid>left&&mid<right,'CUMCM:q3:EventResolution','Event bracket below floating resolution.');
    value=evaluator(mid);g=value.maximum_C-p.threshold;
    if isempty(errorIndicators)
        if g<0,right=mid;rightValue=value;else,left=mid;leftValue=value;end
    else
        epsilon=errorIndicators(mid);assert(isscalar(epsilon)&&isfinite(epsilon)&&epsilon>=0);
        if g-epsilon>=0,left=mid;leftValue=value;
        elseif g+epsilon<0,right=mid;rightValue=value;
        else,ambiguous=true;break
        end
    end
end
stop=p.stopLattice_s*ceil(right/p.stopLattice_s);
stop=round(stop*1000)/1000;
if stop<right,stop=stop+p.stopLattice_s;end
value=evaluator(stop);
epsilon=[];gplus=[];gminus=[];status='CANDIDATE_UNVERIFIED';
if ~isempty(errorIndicators)
    epsilon=errorIndicators(stop);gplus=value.maximum_C-p.threshold+epsilon;
    gminus=leftValue.maximum_C-p.threshold-errorIndicators(left);
    rightPlus=rightValue.maximum_C-p.threshold+errorIndicators(right);
    if ~ambiguous&&gminus>=0&&rightPlus<0&&gplus<0&&right-left<=p.rootWidth_s
        status='NUMERICAL_SCREEN_PASSED_NOT_FULL_VALIDATION';
    else,status='PRECISION_UNRESOLVED';end
end
event=struct('nominal_root_s',left+(right-left)/2,'bracket_left_s',left,'bracket_right_s',right, ...
    'stop_time_s',stop,'stop_time_h',stop/3600,'maximum_C',value.maximum_C, ...
    'epsilon_C',epsilon,'G_minus',gminus,'G_plus',gplus,'status',status, ...
    'strict_nominal_stop',value.maximum_C<p.threshold,'stop_state',value.state);
end

function maximum=q3_spatial_maximum(C,T_K,mesh,boundary,p)
C=C(:);T_K=T_K(:); assert(numel(C)==mesh.N+1&&numel(T_K)==numel(C));
assert(all(isfinite(C)&C>0&T_K>0&isfinite(T_K)));
material=q2_material_properties(C,T_K,p);
d=zeros(size(C)); d(2:end-1)=(C(3:end)-C(1:end-2))/(2*mesh.dr);
d(end)=-p.hm*(C(end)-boundary.moisture_kgkg)/material.D(end);
left=C(1:end-1);right=C(2:end);sl=mesh.dr*d(1:end-1);sr=mesh.dr*d(2:end);
coeff=[2*left-2*right+sl+sr,-3*left+3*right-2*sl-sr,sl,left];
values=C; positions=mesh.r;
for k=1:mesh.N
    polynomial=[3*coeff(k,1),2*coeff(k,2),coeff(k,3)];
    first=find(polynomial~=0,1);
    if isempty(first),continue;end
    stationary=roots(polynomial(first:end)/max(abs(polynomial)));
    realMask=abs(imag(stationary))<=64*eps(max(1,abs(real(stationary))));
    x=real(stationary(realMask)); x=x(x>0&x<1);
    if ~isempty(x)
        values=[values;polyval(coeff(k,:),x)]; %#ok<AGROW>
        positions=[positions;mesh.r(k)+mesh.dr*x]; %#ok<AGROW>
    end
end
peak=max(values); equal=abs(values-peak)<=64*eps(max(1,abs(peak)));
maximum=struct('node_max',max(C),'reconstructed_max',peak,'linear_max',max(C), ...
    'argmax_radius_m',unique(positions(equal)),'overshoot',peak-max(C), ...
    'minimum_reconstructed_candidate',min(values),'coefficients',coeff);
end

function times=q3_output_times(stop_time_s)
assert(isscalar(stop_time_s)&&isfinite(stop_time_s)&&stop_time_s>0);
endMs=round(stop_time_s*1000);
assert(endMs<=flintmax&&abs(endMs/1000-stop_time_s)<=32*eps(max(1,stop_time_s)), ...
    'CUMCM:q3:OutputTime','Stop time must be on the declared millisecond lattice.');
minuteMs=unique([(60000:60000:endMs)';endMs]);
paperMs=unique([(21600000:21600000:endMs)';endMs]);
assert(numel(minuteMs)==ceil(endMs/60000)&&numel(paperMs)==ceil(endMs/21600000));
times=struct('time_s',minuteMs/1000,'paper_time_s',paperMs/1000, ...
    'paper_time_h',paperMs/3600000,'end_ms',endMs);
end

function fields=extract_q3_fields(solution,upstream,event,mesh)
finish=solution.final_time_s;if ~isempty(event),finish=event.stop_time_s;end
times=q3_output_times(finish);
history=upstream.history;
allTime=[history.time_s(:);solution.time_s(:)];
allC=[history.moisture_kgkg;solution.moisture_kgkg];
if ~isempty(event)&&~any(allTime==finish)
    assert(numel(event.stop_state)==2*(mesh.N+1),'CUMCM:q3:MissingEndpoint','Missing true stop state.');
    allTime(end+1)=finish;allC(end+1,:)=event.stop_state(2*mesh.output_columns);
end
assert(numel(unique(allTime))==numel(allTime),'CUMCM:q3:DuplicateHistory','Overlapping upstream and continuation.');
[found,index]=ismember(times.time_s,allTime);assert(all(found),'CUMCM:q3:MissingOutput','Required times were not solved.');
[paperFound,paperIndex]=ismember(times.paper_time_s,times.time_s);assert(all(paperFound));
fields=struct('time_s',times.time_s,'radius_cm',0:.1:2,'moisture_kgkg',allC(index,:), ...
    'paper_time_h',times.paper_time_h,'paper_moisture_kgkg',allC(index(paperIndex),[1,6,11,16,21]), ...
    'end_time_s',finish);
assert(all(isfinite(fields.moisture_kgkg)&fields.moisture_kgkg>0,'all'));
end

function candidate=q3_solve_event_case(cfg,upstream,p,ambient,scenario,caseRoot,totalTimer,prefix)
if nargin<8,prefix=[];end
p.wallSeconds=p.validationWallSeconds-120;p.eventStop=true;
mesh=upstream.mesh;initial=upstream.initial;allTimes=zeros(0,1);allT=zeros(0,21);allC=allT;
allMax=zeros(0,1);traces=strings(0,1);
cases=struct('end_time_s',{},'accepted_half_steps',{},'rejected_macro_steps',{}, ...
    'local_error_max',{},'newton_max_iterations',{},'trace_path',{});accepted=0;rejected=0;
if ~isempty(prefix)
    initial=prefix.initial;allTimes=prefix.time_s;allT=prefix.temperature_K;
    allC=prefix.moisture_kgkg;allMax=prefix.maximum_C;traces=prefix.trace_paths;
end
window=0;first=initial.time_s<14400;
while true
    if first,finish=14400;first=false;else,finish=14400+21600*2^window;window=window+1;end
    if finish<=initial.time_s,continue;end
    assert(isfinite(finish),'CUMCM:q3:TimeOverflow','Search window cannot be represented.');
    relative=string(caseRoot)+"/window_"+numel(cases)+"/accepted_states.h5";
    context=struct('stop_time_s',finish,'trace_path',fullfile(cfg.projectRoot,relative),'total_timer',totalTimer);
    sol=solve_q3_longtime(initial,[],mesh,p,ambient,scenario,context);
    allTimes=[allTimes;sol.time_s];allT=[allT;sol.temperature_K];allC=[allC;sol.moisture_kgkg]; %#ok<AGROW>
    allMax=[allMax;sol.maximum_C];traces(end+1,1)=relative; %#ok<AGROW>
    accepted=accepted+sol.accepted_half_steps;rejected=rejected+sol.rejected_macro_steps;
    cases(end+1)=struct('end_time_s',sol.final_time_s,'accepted_half_steps',sol.accepted_half_steps, ...
        'rejected_macro_steps',sol.rejected_macro_steps,'local_error_max',sol.local_error_max, ...
        'newton_max_iterations',sol.newton_max_iterations,'trace_path',relative); %#ok<AGROW>
    fprintf('q3 %s: reached %.3f h; reconstructed maximum C %.9g; elapsed %.1f s.\n', ...
        scenario.id,sol.final_time_s/3600,maxAt(sol.final_state,sol.final_time_s),toc(totalTimer));
    initial=struct('time_s',sol.final_time_s,'state',sol.final_state);
    if sol.final_time_s==14400&&isempty(prefix)
        prefix=struct('initial',initial,'time_s',allTimes,'temperature_K',allT, ...
            'moisture_kgkg',allC,'maximum_C',allMax,'trace_paths',traces);
    end
    if ~isempty(sol.event_bracket),break;end
end
left=struct('time_s',sol.event_bracket.left_s,'state',sol.event_bracket.left_state);
evaluator=@(t)q3_event_evaluate(t,left,mesh,p,ambient,scenario,totalTimer);
event=locate_q3_event(sol,evaluator,[],p);
assert(event.strict_nominal_stop,'CUMCM:q3:StrictStop','Candidate stop does not satisfy strict threshold.');
stopMaximum=q3_spatial_maximum(event.stop_state(2:2:end),event.stop_state(1:2:end),mesh, ...
    interpolate_q3_ambient(event.stop_time_s,ambient,scenario),p);
event.argmax_radius_m=stopMaximum.argmax_radius_m;
event.nominal_moisture_margin=p.threshold-event.maximum_C;
event.conditional_on_manual_upstream=upstream.conditional_on_manual_upstream;
whole=sol;whole.time_s=allTimes;whole.temperature_K=allT;whole.moisture_kgkg=allC;whole.maximum_C=allMax;
fields=extract_q3_fields(whole,upstream,event,mesh);
candidate=struct('scenario',scenario,'fields',fields,'event',event,'trace_paths',traces, ...
    'windows',cases,'prefix',prefix,'runtime_seconds',toc(totalTimer), ...
    'accepted_half_steps',accepted,'rejected_macro_steps',rejected,'settings',rmfield(p,'projectRoot'));
    function result=maxAt(y,t)
        m=q3_spatial_maximum(y(2:2:end),y(1:2:end),mesh,interpolate_q3_ambient(t,ambient,scenario),p);
        result=m.reconstructed_max;
    end
end

function result = select_equilibrium_window(time_s, values, threshold, alpha, minimum_minutes)
if nargin < 4, alpha = 0.8; end
if nargin < 5, minimum_minutes = 30; end
time_s = time_s(:); values = values(:);
assert(numel(time_s)==numel(values) && numel(values)>=2);
assert(all(isfinite(time_s)) && all(isfinite(values)) && all(diff(time_s)>0));
assert(0<alpha && alpha<1 && minimum_minutes>0);
first = find(values>threshold,1,'first');
assert(~isempty(first),'观测中没有严格超过阈值的记录。');
N = numel(values); nmax = N-first+1;
reference_variance = var(values(first:end),0);
starts = (first:N-1)';
starts = starts(time_s(end)-time_s(starts)>=minimum_minutes*60);
assert(~isempty(starts),'首次越阈后没有足够长的观测窗口。');
stats = zeros(numel(starts),8);
for k=1:numel(starts)
    j=starts(k); x=values(j:end); t=time_s(j:end); n=numel(x);
    s2=var(x,0);
    if reference_variance==0, variance_ratio=0; else, variance_ratio=s2/reference_variance; end
    variance_penalty=alpha*variance_ratio;
    sample_penalty=(1-alpha)*(nmax/n-1);
    stats(k,:)=[time_s(j)/60,n,s2,mean(x),trapz(t,x)/(t(end)-t(1)), ...
        variance_penalty,sample_penalty,variance_penalty+sample_penalty];
end
[~,best]=min(stats(:,8));
result=struct('first_crossing_min',time_s(first)/60, ...
    'first_crossing_value',values(first),'maximum_samples',nmax, ...
    'reference_variance',reference_variance,'alpha',alpha, ...
    'minimum_duration_min',minimum_minutes, ...
    'start_min',stats(best,1),'end_min',time_s(end)/60, ...
    'sample_count',stats(best,2),'sample_variance',stats(best,3), ...
    'sample_mean',stats(best,4),'time_mean',stats(best,5), ...
    'variance_penalty',stats(best,6),'sample_penalty',stats(best,7), ...
    'loss',stats(best,8),'candidate_statistics',stats);
end

function history=tail_maximum_history(root,earlyMaximum,candidate)
history=earlyMaximum;
for k=1:numel(candidate.trace_paths)
    file=fullfile(root,candidate.trace_paths(k));info=h5info(file,'/accepted');
    count=info.Dataspace.Size(2);width=info.Dataspace.Size(1);
    t=h5read(file,'/accepted',[1,2],[1,count-1])';
    maximum=h5read(file,'/accepted',[width,2],[1,count-1])';
    history=[history;t,maximum]; %#ok<AGROW>
end
history=history(history(:,1)<candidate.event.stop_time_s,:);
history(end+1,:)=[candidate.event.stop_time_s,candidate.event.maximum_C];
assert(all(diff(history(:,1))>0)&&all(isfinite(history),'all'));
end
