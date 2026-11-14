function varargout = q1_compute(action,varargin)
action=string(action);
assert(isscalar(action),"CUMCM:q1:InvalidComputeAction","计算动作必须是标量字符串。");
switch action
    case "parameters"
        varargout{1}=q1_parameters();
    case "load_inputs"
        [varargout{1:nargout}]=load_q1_inputs(varargin{:});
    case "build_mesh"
        varargout{1}=build_q1_mesh(varargin{:});
    case "interpolate_ambient"
        [varargout{1:nargout}]=interpolate_q1_ambient(varargin{:});
    case "solve_heat"
        varargout{1}=solve_q1_heat(varargin{:});
    case "solve_moisture"
        varargout{1}=solve_q1_moisture(varargin{:});
    case "moisture_face_flux"
        [varargout{1:nargout}]=q1_moisture_face_flux(varargin{:});
    case "verify_heat_trace"
        varargout{1}=q1_verify_heat_trace(varargin{:});
    case "verify_moisture_trace"
        varargout{1}=q1_verify_moisture_trace(varargin{:});
    case "extract_output_grid"
        varargout{1}=extract_q1_output_grid(varargin{:});
    otherwise
        error("CUMCM:q1:InvalidComputeAction","未知计算动作：%s",action);
end
end

function [ambient, audit] = load_q1_inputs(cfg)
file=fullfile(cfg.rawDataDir,"附件1.xlsx");
assert(isequal(string(sheetnames(file)),"Sheet1"), ...
    "CUMCM:q1:InputStructureMismatch","Unexpected attachment sheets.");
cells=readcell(file,"Sheet","Sheet1");
assert(isequal(size(cells),[242,3]), ...
    "CUMCM:q1:InputStructureMismatch","Expected A1:C242.");
assert(isequal(string(cells(1,:)),["时间","温度","水分浓度"]), ...
    "CUMCM:q1:InputStructureMismatch","Unexpected headers.");
assert(all(cellfun(@(x) isnumeric(x)&&isscalar(x)&&isreal(x)&&isfinite(x),cells(2:end,:)),"all"), ...
    "CUMCM:q1:InputStructureMismatch","Nonfinite or nonnumeric cell.");
data=cell2mat(cells(2:end,:));
assert(isequal(data(:,1),(0:60:14400)'), ...
    "CUMCM:q1:InputStructureMismatch","Unexpected time grid.");
q=data(data(:,1)<=1800,:);
assert(isequal(size(q),[31,3]) && isequal(q(1,:),[0,28,0.01963]) && ...
    isequal(q(end,:),[1800,41.513,0.03307]), ...
    "CUMCM:q1:InputStructureMismatch","q1 endpoints differ.");
ambient=struct("time_s",q(:,1),"temperature_C",q(:,2),"moisture_kgkg",q(:,3));
audit=struct("total_rows",241,"q1_rows",31,"headers_checked",true, ...
    "missing_cells",0,"post_1800_used_in_solver",false);
end

function p = q1_parameters()
p.radius_m=0.02; p.length_m=0.25; p.theta0=28; p.C0=2.55;
p.rho=820; p.cp=2600; p.k=0.36; p.hc=25; p.hm=8e-7;
p.D0=7e-9; p.a=0.89; p.tFinal=1800;
p.displayDelta=5e-5; p.newtonMax=30; p.backtrackMin=2^-20;
p.algebraicTolerance=[5e-7,5e-9];
p.quadRelTol=[1e-10,1e-12]; p.quadAbsTol=[1e-12,1e-14];
p.minStep=2^-20; p.maxSegments=flintmax; p.wallSeconds=3600;
p.maxWorkingBytes=2*1024^3; p.maxTraceBytes=8*1024^3;
end

function mesh = build_q1_mesh(radius_m, radialSegments)
assert(isscalar(radius_m)&&isfinite(radius_m)&&radius_m>0 && ...
    isscalar(radialSegments)&&radialSegments>=20&& ...
    mod(log2(radialSegments/20),1)==0, ...
    "CUMCM:q1:MeshInvariantFailure","Invalid radial mesh.");
dr=radius_m/radialSegments;
radialNodes=(0:radialSegments)'*dr;
radialFaces=((0:radialSegments-1)'+0.5)*dr;
controlVolumeBounds=[0;radialFaces;radius_m];

volumeWeights=diff(controlVolumeBounds.^2)/2;
mesh=struct("r",radialNodes,"dr",dr,"w",volumeWeights,"faces",radialFaces, ...
    "radius",radius_m,"N",radialSegments);
assert(abs(sum(volumeWeights)-radius_m^2/2)<=64*eps(radius_m^2) && ...
    abs(volumeWeights(1)-dr^2/8)<=64*eps(dr^2) && all(volumeWeights>0), ...
    "CUMCM:q1:MeshInvariantFailure","Control volume weights inconsistent.");
end

function [ambientTemperature_C, equivalentMoisture_kgkg] = interpolate_q1_ambient(t_s, ambient)
isInsideObservationWindow=all(t_s>=0 & t_s<=1800,"all");
assert(isreal(t_s)&&all(isfinite(t_s),"all")&&isInsideObservationWindow, ...
    "CUMCM:q1:InputStructureMismatch","Time outside q1 interval.");
isNondecreasing=all(diff(t_s(:))>=0);
assert(isNondecreasing,"CUMCM:q1:InputStructureMismatch","Time not sorted.");
ambientTemperature_C=interp1(ambient.time_s,ambient.temperature_C,t_s,"linear");
equivalentMoisture_kgkg=interp1(ambient.time_s,ambient.moisture_kgkg,t_s,"linear");
end

function solution = solve_q1_heat(p,ambient,radialSegments,maxStep_s,algebraicLevel)
assert(ismember(algebraicLevel,[0,1])&&maxStep_s>=p.minStep&&maxStep_s<=1 && ...
    mod(log2(1/maxStep_s),1)==0 && radialSegments<=p.maxSegments, ...
    "CUMCM:q1:MeshInvariantFailure","Invalid heat discretization.");
mesh=build_q1_mesh(p.radius_m,radialSegments);
nodeCount=radialSegments+1;
stepCount=round(p.tFinal/maxStep_s);
streaming=isfield(p,"trace_path")&&strlength(string(p.trace_path))>0;
if isfield(p,"totalTimer"), started=p.totalTimer; else, started=tic; end
if streaming
    estimatedBytes=2*8*((p.tFinal+1)*nodeCount+256*(nodeCount+2)+20*nodeCount);
else
    estimatedBytes=2*8*((stepCount+1)*(nodeCount+5)+20*nodeCount);
end
assert(estimatedBytes<=p.maxWorkingBytes, ...
    "CUMCM:q1:ResourceLimit","Heat working arrays exceed 2 GiB.");
if streaming
    time=(0:p.tFinal)'; state=zeros(numel(time),nodeCount);
    trace=q1_trace_store(p.trace_path,nodeCount+2,stepCount+1,p);
else
    time=(0:stepCount)'*maxStep_s; state=zeros(stepCount+1,nodeCount);
end
state(1,:)=p.theta0;
previousTemperature_C=p.theta0*ones(nodeCount,1);

faceConductance=p.k*mesh.faces/mesh.dr;
leftNode=(1:nodeCount-1)';
transportMatrix=sparse( ...
    [leftNode;leftNode;leftNode+1;leftNode+1], ...
    [leftNode;leftNode+1;leftNode;leftNode+1], ...
    [faceConductance;-faceConductance;-faceConductance;faceConductance], ...
    nodeCount,nodeCount);
transportMatrix(end,end)=transportMatrix(end,end)+p.radius_m*p.hc;
timeCoefficient=p.rho*p.cp*mesh.w/maxStep_s;
systemMatrix=spdiags(timeCoefficient,0,nodeCount,nodeCount)+transportMatrix;
systemFactor=decomposition(systemMatrix,"lu");
algebraicTolerance=p.algebraicTolerance(algebraicLevel+1);
maxBalance=0; maxCorrection=0; maxResidual=0;
if streaming
    [initialAmbient_C,~]=interpolate_q1_ambient(0,ambient);
    trace.append([0;previousTemperature_C;initialAmbient_C]);
    balance=[];
else
    balance=zeros(stepCount,1);
end
minimum=min([p.theta0;ambient.temperature_C]); maximum=max([p.theta0;ambient.temperature_C]);
for stepIdx=1:stepCount
    currentTime_s=stepIdx*maxStep_s;
    [ambientTemperature_C,~]=interpolate_q1_ambient(currentTime_s,ambient);
    rightHandSide=timeCoefficient.*(previousTemperature_C-p.theta0);
    rightHandSide(end)=rightHandSide(end)+ ...
        p.radius_m*p.hc*(ambientTemperature_C-p.theta0);
    temperatureRise_C=systemFactor\rightHandSide;
    for refineIdx=1:3
        linearResidual=systemMatrix*temperatureRise_C-rightHandSide;
        linearUpdate=-(systemFactor\linearResidual);
        if max(abs(linearUpdate))<=algebraicTolerance && ...
                max(abs(linearResidual)./timeCoefficient)<=algebraicTolerance
            break
        end
        temperatureRise_C=temperatureRise_C+linearUpdate;
    end
    nextTemperature_C=p.theta0+temperatureRise_C;
    correction=max(abs(linearUpdate));
    residual=max(abs(linearResidual)./timeCoefficient);
    maxCorrection=max(maxCorrection,correction); maxResidual=max(maxResidual,residual);
    assert(correction<=algebraicTolerance && residual<=algebraicTolerance, ...
        "CUMCM:q1:SolverNonConvergence","Heat linear correction failed.");
    assert(all(isfinite(nextTemperature_C))&&all(nextTemperature_C>=minimum) && ...
        all(nextTemperature_C<=maximum), ...
        "CUMCM:q1:RangeOrFiniteFailure","Heat range failure.");

    radialBalance=p.rho*p.cp*sum(mesh.w.*(nextTemperature_C-previousTemperature_C))+ ...
        maxStep_s*p.radius_m*p.hc*(nextTemperature_C(end)-ambientTemperature_C);
    energyDefect_J=2*pi*p.length_m*radialBalance;
    maxBalance=max(maxBalance,abs(energyDefect_J));
    assert(abs(radialBalance)<=10*maxStep_s*sum(timeCoefficient)*algebraicTolerance, ...
        "CUMCM:q1:BalanceFailure","Heat balance inconsistent with algebraic tolerance.");
    if streaming
        trace.append([currentTime_s;nextTemperature_C;ambientTemperature_C]);
        if currentTime_s==floor(currentTime_s)
            state(currentTime_s+1,:)=nextTemperature_C';
        end
    else
        state(stepIdx+1,:)=nextTemperature_C';
        balance(stepIdx)=energyDefect_J;
    end
    previousTemperature_C=nextTemperature_C;
    assert(toc(started)<p.wallSeconds,"CUMCM:q1:ResourceLimit","Total heat/validation time limit.");
end
if streaming
    trace.finish();
    balance=maxBalance;
end
solution=struct("time_s",time,"r_m",mesh.r,"state",state,"balance",balance, ...
    "max_correction",maxCorrection,"max_scaled_residual",maxResidual, ...
    "dt_s",maxStep_s,"mesh",mesh,"streaming",streaming);
end

function solution = solve_q1_moisture(p,ambient,radialSegments,maxStep_s,algebraicLevel,quadratureLevel)
assert(ismember(algebraicLevel,[0,1])&&ismember(quadratureLevel,[0,1]) && ...
    maxStep_s>=p.minStep&&maxStep_s<=1&&mod(1/maxStep_s,1)==0 && ...
    radialSegments<=p.maxSegments,"CUMCM:q1:MeshInvariantFailure","Invalid moisture discretization.");
mesh=build_q1_mesh(p.radius_m,radialSegments);
plannedRows=round(p.tFinal/maxStep_s)+1;
streaming=isfield(p,"trace_path")&&strlength(string(p.trace_path))>0;
if streaming
    plannedBytes=2*8*((p.tFinal+1)*(radialSegments+1)+256*(radialSegments+3)+40*(radialSegments+1));
else
    plannedBytes=2*8*plannedRows*(radialSegments+8);
end
assert(plannedBytes<=p.maxWorkingBytes,"CUMCM:q1:ResourceLimit","Moisture work arrays exceed budget.");
time=0; state=p.C0*ones(1,radialSegments+1); balance=zeros(0,1);
iterations=zeros(0,1); corrections=zeros(0,1); residuals=zeros(0,1); rejected=0;
previousMoisture=state'; minimumStep=Inf; maximumStep=0; maximumIterations=0;
maximumBalance=0; maximumCorrection=0; maximumResidual=0; acceptedSteps=0;
rejectedIntervals=zeros(0,2);
if streaming
    trace=q1_trace_store(p.trace_path,radialSegments+3,plannedRows,p);
    [~,initialAmbientMoisture]=interpolate_q1_ambient(0,ambient);
    trace.append([0;previousMoisture;initialAmbientMoisture]);
    time=(0:p.tFinal)'; state=zeros(numel(time),radialSegments+1);
    state(1,:)=previousMoisture';
end
t=0; proposedStep_s=maxStep_s; started=tic;
algebraicTolerance=p.algebraicTolerance(algebraicLevel+1);
if isfield(p,"totalTimer"), started=p.totalTimer; end
minimum=min([p.C0;ambient.moisture_kgkg]); maximum=max([p.C0;ambient.moisture_kgkg]);
while t<p.tFinal
    step_s=min([proposedStep_s,p.tFinal-t,ceil(t+eps(max(1,t)))-t]);
    assert(step_s>=p.minStep,"CUMCM:q1:SolverNonConvergence", ...
        "Rejected step below resource guard.");
    [~,ambientMoisture]=interpolate_q1_ambient(t+step_s,ambient);
    candidate=previousMoisture;
    converged=false; lastCorrection=Inf; cached=false;
    for newtonIdx=1:p.newtonMax
        if ~cached
            [residual,jacobian,audit]=q1_moisture_residual_jacobian( ...
                candidate,previousMoisture,mesh,step_s,ambientMoisture,p,quadratureLevel);
        end
        newtonUpdate=-(jacobian\residual);
        lastCorrection=max(abs(newtonUpdate));
        if lastCorrection<=algebraicTolerance && ...
                audit.scaled_residual<=algebraicTolerance
            converged=true; break
        end
        damping=1; accepted=false;
        while damping>=p.backtrackMin
            trialState=candidate+damping*newtonUpdate;
            if all(isfinite(trialState)&trialState>0)
                [trialResidual,trialJacobian,trialAudit]=q1_moisture_residual_jacobian( ...
                    trialState,previousMoisture,mesh,step_s,ambientMoisture,p,quadratureLevel);
                if trialAudit.scaled_residual<audit.scaled_residual
                    candidate=trialState; residual=trialResidual;
                    jacobian=trialJacobian; audit=trialAudit;
                    cached=true; accepted=true; break
                end
            end
            damping=damping/2;
        end
        if ~accepted, break; end
    end
    if ~converged
        proposedStep_s=step_s/2; rejected=rejected+1;
        rejectedIntervals(end+1,:)=[t,t+step_s];
        assert(toc(started)<p.wallSeconds,"CUMCM:q1:ResourceLimit","Moisture time limit.");
        continue
    end
    assert(all(isfinite(candidate)&candidate>0) && all(candidate>=minimum) && ...
        all(candidate<=maximum), ...
        "CUMCM:q1:RangeOrFiniteFailure","Moisture range failure.");
    assert(abs(audit.balance)<=10*sum(mesh.w)*algebraicTolerance && ...
        abs(audit.balance-audit.weighted_residual)<=1e-18, ...
        "CUMCM:q1:BalanceFailure","Moisture balance mismatch.");
    t=t+step_s; acceptedSteps=acceptedSteps+1;
    minimumStep=min(minimumStep,step_s); maximumStep=max(maximumStep,step_s);
    maximumIterations=max(maximumIterations,newtonIdx);
    maximumBalance=max(maximumBalance,abs(audit.balance));
    maximumCorrection=max(maximumCorrection,lastCorrection);
    maximumResidual=max(maximumResidual,audit.scaled_residual);
    if streaming
        if trace.written+trace.used>=trace.maxRows
            trace.maxRows=trace.maxRows+1;
            assert(trace.otherBytes+8*trace.width*trace.maxRows+1024^2<=p.maxTraceBytes, ...
                "CUMCM:q1:ResourceLimit","Rejected-step trace exceeds budget.");
        end
        trace.append([t;candidate;ambientMoisture]);
        if t==floor(t), state(t+1,:)=candidate'; end
    else
        assert(2*8*(numel(time)+1)*(radialSegments+8)<=p.maxWorkingBytes, ...
            "CUMCM:q1:ResourceLimit","Rejected steps exceed memory budget.");
        time(end+1,1)=t; state(end+1,:)=candidate';
        balance(end+1,1)=audit.balance; iterations(end+1,1)=newtonIdx;
        corrections(end+1,1)=lastCorrection; residuals(end+1,1)=audit.scaled_residual;
    end
    previousMoisture=candidate;
    proposedStep_s=maxStep_s;
    assert(toc(started)<p.wallSeconds,"CUMCM:q1:ResourceLimit","Moisture time limit.");
end
if streaming
    trace.finish(); balance=maximumBalance; iterations=maximumIterations;
    corrections=maximumCorrection; residuals=maximumResidual;
end
solution=struct("time_s",time,"r_m",mesh.r,"state",state,"balance",balance, ...
    "iterations",iterations,"rejected_steps",rejected, ...
    "max_correction",max(corrections),"max_scaled_residual",max(residuals), ...
    "dt_min_s",minimumStep,"dt_max_s",maximumStep,"mesh",mesh, ...
    "streaming",streaming,"accepted_steps",acceptedSteps,"rejected_intervals",rejectedIntervals);
end

function [residual,jacobian,audit] = q1_moisture_residual_jacobian(candidate,previous,mesh,dt_s,ambientAtEnd,p,quadratureLevel)
nodeCount=numel(mesh.w);
assert(isequal(size(candidate),[nodeCount,1])&&isequal(size(previous),[nodeCount,1]), ...
    "CUMCM:q1:MeshInvariantFailure","Wrong state shape.");
[faceFlux,dFlux_dLeft,dFlux_dRight]=q1_moisture_face_flux( ...
    candidate(1:end-1),candidate(2:end),mesh.dr,p,quadratureLevel);
radiallyWeightedFlux=[0;mesh.faces.*faceFlux; ...
    mesh.radius*p.hm*(candidate(end)-ambientAtEnd)];
residual=mesh.w.*(candidate-previous)/dt_s+diff(radiallyWeightedFlux);
timeDiagonal=mesh.w/dt_s;
jacobian=spdiags(timeDiagonal,0,nodeCount,nodeCount);
leftNode=(1:nodeCount-1)';
leftDerivative=mesh.faces.*dFlux_dLeft;
rightDerivative=mesh.faces.*dFlux_dRight;
jacobian=jacobian+sparse( ...
    [leftNode;leftNode;leftNode+1;leftNode+1], ...
    [leftNode;leftNode+1;leftNode;leftNode+1], ...
    [leftDerivative;rightDerivative;-leftDerivative;-rightDerivative], ...
    nodeCount,nodeCount);
jacobian(end,end)=jacobian(end,end)+mesh.radius*p.hm;
audit=struct( ...
    "balance",sum(mesh.w.*(candidate-previous))+dt_s*radiallyWeightedFlux(end), ...
    "weighted_residual",dt_s*sum(residual), ...
    "scaled_residual",max(abs(residual)./timeDiagonal));
end

function [moistureFlux, dFlux_dLeft, dFlux_dRight, audit] = q1_moisture_face_flux(cLeft,cRight,spacing_m,p,quadratureLevel)
assert(isequal(size(cLeft),size(cRight))&&isreal(cLeft)&&isreal(cRight) && ...
    all(isfinite(cLeft)&isfinite(cRight)&cLeft>0&cRight>0,"all"), ...
    "CUMCM:q1:NonPositiveTrial","Flux requires finite positive C.");
assert(isscalar(spacing_m)&&isfinite(spacing_m)&&spacing_m>0 && ...
    ismember(quadratureLevel,[0,1]),"CUMCM:q1:MeshInvariantFailure","Invalid flux options.");
diffusivityLeft=p.D0*exp(-p.a./cLeft);
diffusivityRight=p.D0*exp(-p.a./cRight);
assert(all(diffusivityLeft>0 & diffusivityRight>0,"all"), ...
    "CUMCM:q1:SolverNonConvergence","Diffusivity underflow; no artificial floor.");
concentrationJump=cRight-cLeft;
normalizedIntegral=ones(size(concentrationJump));
nonUniformFace=concentrationJump~=0;
if any(nonUniformFace,"all")
    leftEndpoint=cLeft(nonUniformFace);
    intervalSpan=concentrationJump(nonUniformFace);
    scalingConcentration=max(leftEndpoint,cRight(nonUniformFace));
    scaledIntegrand=@(s) exp(-p.a./(leftEndpoint+s.*intervalSpan)+p.a./scalingConcentration);
    normalizedIntegral(nonUniformFace)=integral(scaledIntegrand,0,1,"ArrayValued",true, ...
        "RelTol",p.quadRelTol(quadratureLevel+1),"AbsTol",p.quadAbsTol(quadratureLevel+1));
end
assert(all(isfinite(normalizedIntegral)&normalizedIntegral>0& ...
    normalizedIntegral<=1+64*eps,"all"), ...
    "CUMCM:q1:SolverNonConvergence","Scaled integral out of bounds.");
moistureFlux=-concentrationJump./spacing_m.* ...
    max(diffusivityLeft,diffusivityRight).*normalizedIntegral;
dFlux_dLeft=diffusivityLeft/spacing_m;
dFlux_dRight=-diffusivityRight/spacing_m;
audit=struct("integral_min",min(normalizedIntegral,[],"all"), ...
    "integral_max",max(normalizedIntegral,[],"all"));
end

function audit = q1_verify_heat_trace(file,p,mesh)
info=h5info(file,"/trajectory"); count=info.Dataspace.Size(2); width=info.Dataspace.Size(1);
last=[]; maximum=0; checked=0;
for first=1:256:count
    block=h5read(file,"/trajectory",[1,first],[width,min(256,count-first+1)]);
    for k=1:size(block,2)
        row=block(:,k);
        if ~isempty(last)
            dt=row(1)-last(1);
            value=2*pi*p.length_m*(p.rho*p.cp*sum(mesh.w.*(row(2:end-1)-last(2:end-1))) ...
                +dt*p.radius_m*p.hc*(row(end-1)-row(end)));
            maximum=max(maximum,abs(value)); checked=checked+1;
        end
        last=row;
    end
end
audit=struct("maximum_balance_J",maximum,"steps",checked);
end

function audit=q1_verify_moisture_trace(file,p,mesh)
info=h5info(file,"/trajectory"); count=info.Dataspace.Size(2); width=info.Dataspace.Size(1);
last=[]; maximum=0; checked=0; maximumTime=0;
for first=1:256:count
    block=h5read(file,"/trajectory",[1,first],[width,min(256,count-first+1)]);
    for k=1:size(block,2)
        row=block(:,k);
        if ~isempty(last)
            dt=row(1)-last(1); assert(dt>0);
            value=sum(mesh.w.*(row(2:end-1)-last(2:end-1))) ...
                +dt*p.radius_m*p.hm*(row(end-1)-row(end));
            if abs(value)>maximum, maximum=abs(value); maximumTime=row(1); end
            checked=checked+1;
        end
        last=row;
    end
end
assert(last(1)==p.tFinal,"CUMCM:q1:OutputMappingFailure","Incomplete moisture trace.");
audit=struct("maximum_balance_kgkg_m2",maximum,"time_s",maximumTime,"steps",checked);
end

function fields = extract_q1_output_grid(heatSolution,moistureSolution,radialSegments)
targetTime_s=(1:1800)';
radialColumnIndex=1+(0:20)*(radialSegments/20);
[hasTemperatureRow,temperatureRowIndex]=ismember(targetTime_s,heatSolution.time_s);
[hasMoistureRow,moistureRowIndex]=ismember(targetTime_s,moistureSolution.time_s);
assert(all(hasTemperatureRow)&all(hasMoistureRow), ...
    "CUMCM:q1:OutputMappingFailure","Missing integer-second node.");
temperature_C=heatSolution.state(temperatureRowIndex,radialColumnIndex);
moisture_kgkg=moistureSolution.state(moistureRowIndex,radialColumnIndex);
assert(isequal(size(temperature_C),[1800,21])&&isequal(size(moisture_kgkg),[1800,21]), ...
    "CUMCM:q1:OutputMappingFailure","Full field shape mismatch.");
fields=struct("time_s",targetTime_s,"radius_cm",(0:20)/10, ...
    "temperature_C",temperature_C,"moisture_kgkg",moisture_kgkg, ...
    "paper_time_s",[100;300;600;900;1200;1500;1800], ...
    "paper_columns",[1,6,11,16,21]);
end
