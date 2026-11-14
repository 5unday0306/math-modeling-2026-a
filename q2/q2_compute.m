function varargout=q2_compute(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q2:ComputeAction','Invalid compute action.');
switch action
    case "parameters"
        varargout{1}=q2_parameters();
    case "load_inputs"
        [varargout{1:nargout}]=load_q2_inputs(varargin{:});
    case "build_mesh"
        varargout{1}=build_q2_mesh(varargin{:});
    case "interpolate_ambient"
        [varargout{1:nargout}]=interpolate_q2_ambient(varargin{:});
    case "material_properties"
        varargout{1}=q2_material_properties(varargin{:});
    case "face_coefficients"
        varargout{1}=q2_face_coefficients(varargin{:});
    case "residual_jacobian"
        [varargout{1:nargout}]=q2_residual_jacobian(varargin{:});
    case "solve"
        varargout{1}=solve_q2_coupled(varargin{:});
    case "extract_output_grid"
        varargout{1}=extract_q2_output_grid(varargin{:});
    case "verify_trace"
        varargout{1}=q2_verify_trace(varargin{:});
    case "balance_summary"
        varargout{1}=q2_balance_summary(varargin{:});
    case "compare_fields"
        varargout{1}=q2_compare_fields(varargin{:});
    otherwise
        error('CUMCM:q2:ComputeAction','Unknown compute action: %s',action);
end
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

function [ambient,audit]=load_q2_inputs(cfg)
file=fullfile(cfg.rawDataDir,'附件1.xlsx');
assert(isequal(string(sheetnames(file)),"Sheet1"),'CUMCM:q2:InputStructure','Unexpected sheets.');
cells=readcell(file,'Sheet','Sheet1');
assert(isequal(size(cells),[242,3]) && ...
    isequal(string(cells(1,:)),["时间","温度","水分浓度"]),'CUMCM:q2:InputStructure','Wrong attachment structure.');
assert(all(cellfun(@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x),cells(2:end,:)),'all'), ...
    'CUMCM:q2:InputStructure','Missing or nonnumeric input.');
data=cell2mat(cells(2:end,:));
assert(isequal(data(:,1),(0:60:14400)') && isequal(data(1,:),[0,28,0.01963]) && ...
    isequal(data(181,:),[10800,50.195,0.04977]) && isequal(data(end,:),[14400,50.165,0.04986]), ...
    'CUMCM:q2:InputStructure','Input times/endpoints differ.');
template=fullfile(cfg.rawDataDir,'附件3','result2.xlsx');
names=string(sheetnames(template)); assert(isequal(names(:),["温度";"水分浓度"]));
for k=1:2
    first=readcell(template,'Sheet',names(k),'Range','A1','UseExcel',false);
    assert(string(first{1})=="时间\到药材中心的距离",'CUMCM:q2:InputStructure','Template header differs.');
end
ambient=struct('time_s',data(:,1),'temperature_C',data(:,2), ...
    'moisture_kgkg',data(:,3),'q2_indices',(1:181)', ...
    'boundary_version','LINEAR_TO_4H_THEN_LAST_VALUE');
audit=struct('total_rows',241,'q2_rows',181,'missing_cells',0, ...
    'time_interval_s',60,'template_read_only',true,'post_3h_used_in_3h_interpolation',false);
end

function mesh=build_q2_mesh(radius_m,radialSegments)
assert(isscalar(radius_m)&&isfinite(radius_m)&&radius_m>0 && ...
    isscalar(radialSegments)&&radialSegments>=20&&radialSegments<=flintmax && ...
    mod(log2(radialSegments/20),1)==0,'CUMCM:q2:Mesh','Invalid grid.');
dr=radius_m/radialSegments;
radialNodes=(0:radialSegments)'*dr;
radialFaces=((0:radialSegments-1)'+0.5)*dr;
controlVolumeBounds=[0;radialFaces;radius_m];
volumeWeights=diff(controlVolumeBounds.^2)/2;
assert(all(volumeWeights>0) && ...
    abs(sum(volumeWeights)-radius_m^2/2)<=64*eps(radius_m^2) && ...
    abs(volumeWeights(1)-dr^2/8)<=64*eps(dr^2), ...
    'CUMCM:q2:Mesh','Volume invariant failed.');
outputColumns=1+(0:20)*(radialSegments/20);
mesh=struct('r',radialNodes,'faces',radialFaces, ...
    'boundaries',controlVolumeBounds,'w',volumeWeights,'dr',dr, ...
    'radius',radius_m,'N',radialSegments,'output_columns',outputColumns);
end

function [temperature_C,moisture_kgkg]=interpolate_q2_ambient(time_s,ambient)
assert(isreal(time_s)&&all(isfinite(time_s)&time_s>=0,'all'), ...
    'CUMCM:q2:BoundaryTime','Invalid time.');
requestedShape=size(time_s);
requestedTime_s=time_s(:);
observationTime_s=ambient.time_s;
assert(observationTime_s(1)==0&&all(diff(observationTime_s)==60), ...
    'CUMCM:q2:BoundaryInput','Boundary knots differ.');
if isfield(ambient,'validation_interpolation')
    assert(string(ambient.validation_interpolation)=="PCHIP_SENSITIVITY_ONLY", ...
        'CUMCM:q2:BoundaryInput','Unsupported interpolation experiment.');
    boundedTime_s=min(requestedTime_s,observationTime_s(end));
    temperature_C=reshape(ppval(ambient.temperature_pp,boundedTime_s),requestedShape);
    moisture_kgkg=reshape(ppval(ambient.moisture_pp,boundedTime_s),requestedShape);
    return
end
boundedTime_s=min(requestedTime_s,observationTime_s(end));
leftIndex=floor(boundedTime_s/60)+1;
leftIndex=min(leftIndex,numel(observationTime_s)-1);
fraction=(boundedTime_s-observationTime_s(leftIndex))/60;
temperature_C=reshape((1-fraction).*ambient.temperature_C(leftIndex)+ ...
    fraction.*ambient.temperature_C(leftIndex+1),requestedShape);
moisture_kgkg=reshape((1-fraction).*ambient.moisture_kgkg(leftIndex)+ ...
    fraction.*ambient.moisture_kgkg(leftIndex+1),requestedShape);
end

function material=q2_material_properties(moisture,temperature_K,p)
validState=isreal(moisture)&&isreal(temperature_K) && ...
    isequal(size(moisture),size(temperature_K)) && ...
    all(isfinite(moisture)&isfinite(temperature_K)&moisture>0&temperature_K>0,'all');
assert(validState, ...
    'CUMCM:q2:MaterialDomain','Material states must be positive and finite.');
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
    'CUMCM:q2:MaterialDomain','Nonpositive or nonfinite coefficients.');
end

function face=q2_face_coefficients(moistureLeft,moistureRight,temperatureLeft_K,temperatureRight_K,p,quadratureLevel)
sameShape=isequal(size(moistureLeft),size(moistureRight), ...
    size(temperatureLeft_K),size(temperatureRight_K));
assert(ismember(quadratureLevel,[0,1])&&sameShape);
faceStates=[moistureLeft;moistureRight;temperatureLeft_K;temperatureRight_K];
assert(all(isfinite(faceStates)&faceStates>0,'all'), ...
    'CUMCM:q2:MaterialDomain','Face state outside positive domain.');
maximumMoisture=max(moistureLeft,moistureRight);
minimumMoisture=min(moistureLeft,moistureRight);
maximumTemperature_K=max(temperatureLeft_K,temperatureRight_K);
minimumTemperature_K=min(temperatureLeft_K,temperatureRight_K);
scalingMaterial=q2_material_properties(maximumMoisture,maximumTemperature_K,p);
conductivityScale=scalingMaterial.k;
diffusivityScale=scalingMaterial.D;
conductivityDerivativeScale=p.kC./(1+minimumMoisture).^2;
moistureDerivativeScale=diffusivityScale.*p.a./minimumMoisture.^2;
temperatureDerivativeScale=diffusivityScale.*p.B./minimumTemperature_K.^2;
moistureJump=moistureRight-moistureLeft;
temperatureJump_K=temperatureRight_K-temperatureLeft_K;
[oldMessage,oldId]=lastwarn; cleanup=onCleanup(@()lastwarn(oldMessage,oldId));
lastwarn('');
integrals=integral(@integrand,0,1,'ArrayValued',true, ...
    'RelTol',p.quadRelTol(quadratureLevel+1),'AbsTol',p.quadAbsTol(quadratureLevel+1));
[message,warningId]=lastwarn;
assert(isempty(message),'CUMCM:q2:QuadratureFailure','Face quadrature warning %s: %s',warningId,message);
face=struct('k',integrals(:,1).*conductivityScale, ...
    'D',integrals(:,2).*diffusivityScale, ...
    'k_CL',integrals(:,3).*conductivityDerivativeScale, ...
    'k_CR',integrals(:,4).*conductivityDerivativeScale, ...
    'D_CL',integrals(:,5).*moistureDerivativeScale, ...
    'D_CR',integrals(:,6).*moistureDerivativeScale, ...
    'D_TL',integrals(:,7).*temperatureDerivativeScale, ...
    'D_TR',integrals(:,8).*temperatureDerivativeScale);
assert(all(face.k>0&face.D>0),'CUMCM:q2:QuadratureFailure','Nonpositive face mean.');
    function values=integrand(s)
        moisture=moistureLeft+s*moistureJump;
        temperature_K=temperatureLeft_K+s*temperatureJump_K;
        conductivity=p.k0+p.kC*moisture./(1+moisture);
        diffusivity=p.D0*exp(-p.a./moisture-p.B./temperature_K);
        conductivityDerivative=p.kC./(1+moisture).^2;
        moistureDerivative=diffusivity.*p.a./moisture.^2;
        temperatureDerivative=diffusivity.*p.B./temperature_K.^2;
        values=[conductivity./conductivityScale,diffusivity./diffusivityScale, ...
            conductivityDerivative.*(1-s)./conductivityDerivativeScale, ...
            conductivityDerivative.*s./conductivityDerivativeScale, ...
            moistureDerivative.*(1-s)./moistureDerivativeScale, ...
            moistureDerivative.*s./moistureDerivativeScale, ...
            temperatureDerivative.*(1-s)./temperatureDerivativeScale, ...
            temperatureDerivative.*s./temperatureDerivativeScale];
    end
end

function [residual,jacobian,audit]=q2_residual_jacobian(candidate,previous,mesh,dt_s,boundary,p,quadratureLevel)
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
material=q2_material_properties(moisture,temperature_K,p);
face=q2_face_coefficients(moisture(1:end-1),moisture(2:end), ...
    temperature_K(1:end-1),temperature_K(2:end),p,quadratureLevel);

temperatureJump_K=diff(temperature_K);
moistureJump=diff(moisture);
radialGeometry=mesh.faces/mesh.dr;
heatFaceFlux=-radialGeometry.*face.k.*temperatureJump_K;
moistureFaceFlux=-radialGeometry.*face.D.*moistureJump;
surfaceHeatFlux=mesh.radius*p.hc*(temperature_K(end)-boundary.temperature_K);
surfaceMoistureFlux=mesh.radius*p.hm*(moisture(end)-boundary.moisture_kgkg);
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
    mesh.w/dt_s;mesh.radius*p.hc;mesh.radius*p.hm], ...
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

function solution=solve_q2_coupled(p,ambient,mesh,settings)
maximumStep_s=settings.maxStep_s;
algebraicLevel=settings.algebraicLevel;
quadratureLevel=settings.quadratureLevel;
assert(ismember(algebraicLevel,[0,1])&&ismember(quadratureLevel,[0,1]) && ...
    maximumStep_s>=p.minStep&&maximumStep_s<=1 && ...
    mod(log2(1/maximumStep_s),1)==0 && p.tFinal==floor(p.tFinal), ...
    'CUMCM:q2:Settings','Invalid time grid/settings.');
nodeCount=mesh.N+1;
systemSize=2*nodeCount;
temperatureIndex=(1:2:systemSize)';
moistureIndex=temperatureIndex+1;
assert(2*8*((p.tFinal+1)*42+256*(systemSize+12)+100*systemSize)<=p.maxWorkingBytes, ...
    'CUMCM:q2:ResourceLimit','Planned working memory exceeds budget.');
timer=tic;
if isfield(settings,'totalTimer'), timer=settings.totalTimer; end
algebraicTolerance=p.algebraicTolerance(algebraicLevel+1);
t=0;
stateVector=zeros(systemSize,1);
stateVector(temperatureIndex)=p.T0_K;
stateVector(moistureIndex)=p.C0;
temperatureHistory_K=zeros(p.tFinal+1,21);
moistureHistory=temperatureHistory_K;
temperatureHistory_K(1,:)=p.T0_K;
moistureHistory(1,:)=p.C0;
streaming=isfield(settings,'trace_path')&&strlength(string(settings.trace_path))>0;
if streaming
    trace=struct('path',settings.trace_path,'maxBytes',p.maxTraceBytes);
    trace=q2_trace_store('init',trace,[0;stateVector;zeros(12,1)]);
    memoryTrace=[];
else
    memoryTrace=[0;stateVector]';
end
iterationsMax=0; rejected=0; accepted=0; rejectedIntervals=zeros(0,2);
correctionMax=[0,0]; residualMax=[0,0];
balanceMax=[0,0]; balanceSumError=[0,0];
minimumStep=Inf; maximumStep=0;
proposedStep_s=maximumStep_s; backtracksMax=0;
used=ambient.time_s<=p.tFinal;
minimumTemperature_K=min([p.T0_K;ambient.temperature_C(used)+273.15]);
maximumTemperature_K=max([p.T0_K;ambient.temperature_C(used)+273.15]);
minimumMoisture=min([p.C0;ambient.moisture_kgkg(used)]);
maximumMoisture=max([p.C0;ambient.moisture_kgkg(used)]);
while t<p.tFinal
    assert(toc(timer)<=p.wallSeconds,'CUMCM:q2:ResourceLimit','MVP wall-clock budget exceeded.');
    step_s=min([proposedStep_s,p.tFinal-t,ceil(t+eps(max(t,1)))-t]);
    assert(step_s>=p.minStep&&t+step_s>t, ...
        'CUMCM:q2:SolverNonConvergence','Cannot refine step further.');
    [ambientTemperature_C,ambientMoisture]=interpolate_q2_ambient(t+step_s,ambient);
    boundary=struct('temperature_K',ambientTemperature_C+273.15, ...
        'moisture_kgkg',ambientMoisture);
    previousState=stateVector;
    candidate=previousState;
    converged=false; cached=false; backtracks=0;
    for newtonIteration=1:p.newtonMax
        if ~cached
            [residual,jacobian,audit]=q2_residual_jacobian( ...
                candidate,previousState,mesh,step_s,boundary,p,quadratureLevel);
        end
        rowScale=audit.row_scale;
        scaledJacobian=spdiags(rowScale,0,systemSize,systemSize)*jacobian;
        [oldWarning,oldWarningId]=lastwarn;
        lastwarn('');
        newtonUpdate=-(scaledJacobian\(rowScale.*residual));
        [linearWarning,~]=lastwarn;
        lastwarn(oldWarning,oldWarningId);
        if ~isempty(linearWarning)||any(~isfinite(newtonUpdate))
            break
        end
        correction=[max(abs(newtonUpdate(temperatureIndex))), ...
            max(abs(newtonUpdate(moistureIndex)))];
        if all(correction<=algebraicTolerance) && ...
                all(audit.scaled_residual<=algebraicTolerance)
            converged=true;
            break
        end
        residualNormBefore=max(audit.scaled_residual);
        damping=1;
        acceptedTrial=false;
        while damping>=p.backtrackMin
            trialState=candidate+damping*newtonUpdate;
            if all(isfinite(trialState)&trialState>0)
                [trialResidual,trialJacobian,trialAudit]=q2_residual_jacobian( ...
                    trialState,previousState,mesh,step_s,boundary,p,quadratureLevel);
                armijoBound=(1-p.armijo*damping)*residualNormBefore;
                if max(trialAudit.scaled_residual)<=armijoBound
                    candidate=trialState;
                    residual=trialResidual;
                    jacobian=trialJacobian;
                    audit=trialAudit;
                    cached=true;
                    acceptedTrial=true;
                    break
                end
            end
            damping=damping/2;
            backtracks=backtracks+1;
        end
        if ~acceptedTrial, break; end
    end
    if ~converged
        rejected=rejected+1;
        rejectedIntervals(end+1,:)=[t,t+step_s];
        proposedStep_s=step_s/2;
        continue
    end
    temperature_K=candidate(temperatureIndex);
    moisture=candidate(moistureIndex);
    temperatureRoundoff=64*eps(max(abs([minimumTemperature_K,maximumTemperature_K])));
    moistureRoundoff=64*eps(max(abs([minimumMoisture,maximumMoisture])));
    assert(all(temperature_K>=minimumTemperature_K-temperatureRoundoff & ...
        temperature_K<=maximumTemperature_K+temperatureRoundoff) && ...
        all(moisture>=minimumMoisture-moistureRoundoff & ...
        moisture<=maximumMoisture+moistureRoundoff), ...
        'CUMCM:q2:Range','State outside input-history bounds. No clipping is performed.');
    balanceCancellationError=abs(audit.balance-audit.summed_residual);
    assert(balanceCancellationError(1)<=128*eps(max(1,sum(abs(audit.heat_face_flux)))) && ...
        balanceCancellationError(2)<=1e-18, ...
        'CUMCM:q2:Balance','Internal face cancellation failed.');
    assert(abs(audit.balance(1))<= ...
        sum(mesh.w.*audit.material.b)*algebraicTolerance/step_s+1e-10 && ...
        abs(audit.balance(2))<=sum(mesh.w)*algebraicTolerance/step_s+1e-18, ...
        'CUMCM:q2:Balance','Balance exceeds algebraic bound.');
    t=t+step_s;
    stateVector=candidate;
    accepted=accepted+1;
    proposedStep_s=maximumStep_s;
    iterationsMax=max(iterationsMax,newtonIteration);
    backtracksMax=max(backtracksMax,backtracks);
    correctionMax=max(correctionMax,correction);
    residualMax=max(residualMax,audit.scaled_residual);
    balanceMax=max(balanceMax,abs(audit.balance));
    balanceSumError=max(balanceSumError,balanceCancellationError);
    minimumStep=min(minimumStep,step_s);
    maximumStep=max(maximumStep,step_s);
    if streaming
        trace=q2_trace_store('append',trace,[t;stateVector;step_s;newtonIteration; ...
            correction';audit.scaled_residual';audit.balance'; ...
            min(temperature_K);max(temperature_K);min(moisture);max(moisture)]);
    else
        assert(2*8*(size(memoryTrace,1)+1)*(systemSize+1)<=p.maxWorkingBytes, ...
            'CUMCM:q2:ResourceLimit','Memory trace exceeds budget.');
        memoryTrace(end+1,:)=[t;stateVector]';
    end
    if t==floor(t)
        temperatureHistory_K(t+1,:)=temperature_K(mesh.output_columns)';
        moistureHistory(t+1,:)=moisture(mesh.output_columns)';
    end
end
assert(t==p.tFinal&&all(temperatureHistory_K>0,'all') && ...
    all(moistureHistory>0,'all'),'CUMCM:q2:MissingOutput','Missing integer output.');
if streaming
    trace=q2_trace_store('finalize',trace,[]);
else
    trace=struct('memory',memoryTrace);
end
solution=struct('time_s',(0:p.tFinal)','temperature_K',temperatureHistory_K, ...
    'moisture_kgkg',moistureHistory,'final_state',stateVector, ...
    'final_time_s',t,'radius_m',mesh.r,'trace',trace,'streaming',streaming, ...
    'accepted_steps',accepted,'rejected_steps',rejected,'rejected_intervals',rejectedIntervals, ...
    'newton_max_iterations',iterationsMax,'backtracks_max',backtracksMax, ...
    'max_correction',correctionMax,'max_scaled_residual',residualMax, ...
    'max_balance',balanceMax,'max_balance_sum_difference',balanceSumError, ...
    'dt_min_s',minimumStep,'dt_max_s',maximumStep,'runtime_seconds',toc(timer));
end

function fields=extract_q2_output_grid(solution,mesh)
assert(solution.final_time_s==10800&&isequal(solution.time_s,(0:10800)'), ...
    'CUMCM:q2:OutputMapping','Not the required full time horizon.');
assert(isequal(size(solution.temperature_K),[10801,21])&&isequal(size(solution.moisture_kgkg),[10801,21]));
outputRadius_cm=mesh.r(mesh.output_columns)'*100;
fields=struct('time_s',(1:10800)','radius_cm',outputRadius_cm, ...
    'temperature_C',solution.temperature_K(2:end,:)-273.15, ...
    'moisture_kgkg',solution.moisture_kgkg(2:end,:), ...
    'paper_time_s',(1800:1800:10800)','paper_columns',[1,6,11,16,21]);
assert(all(isfinite(fields.temperature_C),'all')&&all(isfinite(fields.moisture_kgkg),'all'));
end

function audit=q2_verify_trace(trace,p,mesh,ambient)
assert(h5readatt(trace.path,'/accepted','complete')==1,'CUMCM:q2:TraceIncomplete','Incomplete trace.');
info=h5info(trace.path,'/accepted'); count=info.Dataspace.Size(2);d=2*(mesh.N+1);
ti=2:2:d+1;ci=ti+1; maximumBalance=[0,0]; recordedDifference=[0,0];
range=[Inf,-Inf,Inf,-Inf]; steps=0; previous=[];oldTime=0;
for first=1:256:count
    n=min(256,count-first+1); block=h5read(trace.path,'/accepted',[1,first],[d+13,n]);
    for j=1:n
        t=block(1,j); T=block(ti,j);C=block(ci,j);
        assert(all(isfinite([T;C])&[T;C]>0),'CUMCM:q2:TraceRange','Bad recorded state.');
        props=q2_material_properties(C,T,p);
        range=[min(range(1),min(T)),max(range(2),max(T)),min(range(3),min(C)),max(range(4),max(C))];
        if isempty(previous)
            assert(t==0&&all(T==p.T0_K)&&all(C==p.C0),'CUMCM:q2:TraceInitial','Initial state differs.');
        else
            dt=t-oldTime;assert(dt>0&&dt==block(d+2,j),'CUMCM:q2:TraceTime','Invalid accepted time.');
            [Ta,Ce]=interpolate_q2_ambient(t,ambient);
            BT=sum(mesh.w.*props.b.*(T-previous(1:2:end))/dt)+mesh.radius*p.hc*(T(end)-(Ta+273.15));
            BC=sum(mesh.w.*(C-previous(2:2:end))/dt)+mesh.radius*p.hm*(C(end)-Ce);
            maximumBalance=max(maximumBalance,abs([BT,BC]));
            recordedDifference=max(recordedDifference,abs([BT,BC]-block(d+8:d+9,j)'));
            steps=steps+1;
        end
        previous=block(2:d+1,j);oldTime=t;
    end
end
assert(oldTime==p.tFinal&&steps==count-1 && recordedDifference(1)<1e-10 && recordedDifference(2)<1e-18, ...
    'CUMCM:q2:TraceBalance','Re-read evidence differs.');
audit=struct('accepted_states_reread',count,'accepted_steps',steps,'last_time_s',oldTime, ...
    'maximum_balance_rate',maximumBalance,'recorded_balance_difference',recordedDifference, ...
    'temperature_range_K',range(1:2),'moisture_range_kgkg',range(3:4), ...
    'heat_balance_unit','W/m','moisture_balance_unit','(kg/kg)*m^2/s', ...
    'complete_trace',true,'all_states_finite_positive',true);
end

function report=q2_balance_summary(trace,p,mesh,ambient)
info=h5info(trace.path,'/accepted');count=info.Dataspace.Size(2);d=2*(mesh.N+1);
heatStorage=0;wetChange=0;heatBE=0;wetBE=0;heatTrap=0;wetTrap=0;previous=[];
for first=1:256:count
    block=h5read(trace.path,'/accepted',[1,first],[d+13,min(256,count-first+1)]);
    for j=1:size(block,2)
        t=block(1,j);T=block(2:2:d+1,j);C=block(3:2:d+1,j);
        material=q2_material_properties(C,T,p);[Ta,Ce]=interpolate_q2_ambient(t,ambient);
        flux=[mesh.radius*p.hc*(Ta+273.15-T(end)),mesh.radius*p.hm*(C(end)-Ce)];
        if ~isempty(previous)
            dt=t-oldTime;
            heatStorage=heatStorage+sum(mesh.w.*material.b.*(T-previous(1:2:end)));
            wetChange=wetChange+sum(mesh.w.*(previous(2:2:end)-C));
            heatBE=heatBE+dt*flux(1);wetBE=wetBE+dt*flux(2);
            heatTrap=heatTrap+dt*(flux(1)+oldFlux(1))/2;
            wetTrap=wetTrap+dt*(flux(2)+oldFlux(2))/2;
        end
        previous=block(2:d+1,j);oldTime=t;oldFlux=flux;
    end
end
denom=[abs(heatStorage),abs(wetChange)];be=abs([heatStorage-heatBE,wetChange-wetBE]);
trap=abs([heatStorage-heatTrap,wetChange-wetTrap]);relative=NaN(1,2);valid=denom>1e-14;
relative(valid)=trap(valid)./denom(valid);
report=struct('heat_PDE_storage_J_per_m',heatStorage,'effective_moisture_loss_kgkg_m2',wetChange, ...
    'backward_euler_flux_integrals',[heatBE,wetBE],'trapezoidal_flux_integrals',[heatTrap,wetTrap], ...
    'discrete_absolute_defects',be,'trapezoidal_absolute_defects',trap, ...
    'trapezoidal_relative_defects',relative,'relative_denominator_valid',valid, ...
    'not_total_multiphase_energy_or_true_water_mass',true,'states_reread',count);
end

function report=q2_compare_fields(a,b)
names={'temperature_C','moisture_kgkg'};maximum=zeros(1,2);count=maximum;locations=cell(1,2);
for j=1:2
    first=a.(names{j});second=b.(names{j});assert(isequal(size(first),size(second)));
    [maximum(j),ix]=max(abs(first-second),[],'all','linear');[row,col]=ind2sub(size(first),ix);
    unequal=round(first,4)~=round(second,4);count(j)=nnz(unequal);
    selected=find(unequal,10);[rr,cc]=ind2sub(size(first),selected);
    locations{j}=struct('maximum_at_time_s',row,'maximum_at_radius_cm',(col-1)/10, ...
        'first_value_at_maximum',first(ix),'second_value_at_maximum',second(ix), ...
        'first_10_rounding_difference_time_s',rr,'first_10_rounding_difference_radius_cm',(cc-1)/10, ...
        'first_10_values_a',first(selected),'first_10_values_b',second(selected));
end
report=struct('maximum_difference',maximum,'rounded_different_count',count, ...
    'locations',{locations},'all_453600_compared',sum(numel(a.temperature_C)+numel(a.moisture_kgkg))==453600, ...
    'rounded_equal',all(count==0));
end
