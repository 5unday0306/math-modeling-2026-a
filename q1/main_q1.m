function output = main_q1(mode)
arguments
    mode (1,1) string = "reproduce"
end
assert(ismember(mode,["mvp","reproduce","validation"]), ...
    "CUMCM:q1:InvalidMode","mode 必须为 mvp、reproduce 或 validation。");

projectRoot=string(fileparts(mfilename("fullpath")));
addpath(projectRoot);
cfg=q1_config(projectRoot);
output=run_q1(cfg,mode);
end

function cfg = q1_config(projectRoot)
arguments
    projectRoot (1,1) string
end

cfg.problem = "A";
cfg.projectRoot = projectRoot;
cfg.rawDataDir = fullfile(projectRoot, "data", "raw");
cfg.processedDataDir = fullfile(projectRoot, "data", "processed");
cfg.resultsDir = fullfile(projectRoot, "results");
cfg.figuresDir = fullfile(projectRoot, "figures");
cfg.figureQADir = fullfile(projectRoot, "artifacts", "validation", "figure_previews");
cfg.scriptsDir = fullfile(projectRoot, "scripts");
cfg.testsDir = fullfile(projectRoot, "tests");
cfg.logDir = fullfile(projectRoot, "artifacts", "logs");
cfg.validationDir = fullfile(projectRoot, "artifacts", "validation");
cfg.questions = strings(1,0);
cfg.randomSeed = 42;
cfg.requiredFeatures = ["data", "visualization"];
cfg.requiredFunctions = ["exportgraphics", "jsonencode", "rng", "readcell", ...
    "readmatrix", "sheetnames", "interp1", "spdiags", "integral", "expint", ...
    "h5create", "h5write", "h5read", "h5info", "decomposition", ...
    "pdepe", "pdeval", "ode15s", "odeset", "besselj", "fzero", "pchip"];
cfg.requiredPythonModules = strings(1,0);
cfg.publicReproductionCommand = "main_q1(""reproduce"")";
cfg.questionCommandTemplate = "main_q1(""<mode>"")";
end

function output = run_q1(cfg,mode)
assert(ismember(string(mode),["mvp","reproduce","validation"]),"CUMCM:q1:ModelContractConflict", ...
    "Final mode requires completed validation.");
p=q1_compute("parameters");
if string(mode)=="validation"
    [ambient,~]=q1_compute("load_inputs",cfg);
    output=q1_validation("resources",cfg,p,ambient);
    return
end
for name=["readcell","readmatrix","interp1","spdiags","integral","expint"]
    assert(~isempty(which(char(name))),"CUMCM:q1:InputStructureMismatch","Missing required MATLAB function.");
end
assert(license('test','MATLAB'),"CUMCM:q1:InputStructureMismatch","Missing MATLAB license.");
[ambient,inputAudit]=q1_compute("load_inputs",cfg);
if string(mode)=="reproduce"
    output=reproduce_q1(cfg,p,ambient,inputAudit);
    return
end
fprintf('q1 MVP: 31 environmental rows; N=20; dt=1 s; horizon=1800 s.\n');
h=q1_compute("solve_heat",p,ambient,20,1,0);
c=q1_compute("solve_moisture",p,ambient,20,1,0,0);
f=q1_compute("extract_output_grid",h,c,20);
candidate=struct("heat",h,"moisture",c,"fields",f);
checks=q1_validation("solution",cfg,p,ambient,candidate,"mvp");
checks.input_audit=inputAudit;
checks.code_snapshot=codeRecords(cfg);
sampleTimes=[100,900,1800]; sampleColumns=[1,11,21]; values=zeros(9,4); row=0;
for t=sampleTimes
    for col=sampleColumns
        row=row+1; values(row,:)=[t,f.radius_cm(col),f.temperature_C(t,col),f.moisture_kgkg(t,col)];
    end
end
folder=fullfile(cfg.projectRoot,"artifacts","q1","mvp");
if ~isfolder(folder), mkdir(folder); end
csvPath=fullfile(folder,"q1_mvp_slice.csv");
fid=fopen(csvPath,"w","n","UTF-8");
assert(fid>=0,"CUMCM:q1:OutputMappingFailure","Cannot write MVP CSV.");
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,"时间（秒）,半径（厘米）,温度（℃）,干基含水率（kg/kg）\n");
fprintf(fid,"%.17g,%.17g,%.17g,%.17g\n",values');
clear cleanup
actual=readmatrix(csvPath,"NumHeaderLines",1);
headers=readcell(csvPath,"Range","A1:D1","Encoding","UTF-8");
assert(isequal(string(headers),["时间（秒）","半径（厘米）","温度（℃）","干基含水率（kg/kg）"]), ...
    "CUMCM:q1:OutputMappingFailure","MVP Chinese headers mismatch.");
assert(isequal(actual,values),"CUMCM:q1:OutputMappingFailure","MVP CSV round-trip mismatch.");
checks.slice_rows=9; checks.slice_values=values;
jsonPath=fullfile(folder,"q1_mvp_checks.json");
fid=fopen(jsonPath,"w","n","UTF-8");
assert(fid>=0,"CUMCM:q1:OutputMappingFailure","Cannot write MVP checks.");
cleanup=onCleanup(@() fclose(fid)); fprintf(fid,"%s\n",jsonencode(checks,"PrettyPrint",true)); clear cleanup
output=struct("artifacts",["artifacts/q1/mvp/q1_mvp_slice.csv", ...
    "artifacts/q1/mvp/q1_mvp_checks.json"],"metrics",checks);
disp(values); fprintf('q1 MVP checks complete; four-decimal stability remains unverified.\n');
end

function output=reproduce_q1(cfg,p,ambient,inputAudit)
N=10240; dt=0.25; algebraicLevel=1; quadratureLevel=0;
traceParent=fullfile(cfg.projectRoot,"artifacts","q1","reproduce_traces");
if ~isfolder(traceParent), mkdir(traceParent); end
uuid=char(java.util.UUID.randomUUID());
traceDir=fullfile(traceParent,"reproduce_"+string(uuid(1:8)));
mkdir(traceDir);
cleanup=onCleanup(@()removeTraceDir(traceDir,traceParent));
p.totalTimer=tic;
fprintf('q1 reproduction: N=%d; dt=%.2f s; horizon=%d s.\n',N,dt,p.tFinal);
trial=p; trial.trace_path=fullfile(traceDir,"heat.h5");
h=q1_compute("solve_heat",trial,ambient,N,dt,algebraicLevel);
heatAudit=q1_compute("verify_heat_trace",trial.trace_path,trial,h.mesh);
trial.trace_path=fullfile(traceDir,"moisture.h5");
c=q1_compute("solve_moisture",trial,ambient,N,dt,algebraicLevel,quadratureLevel);
moistureAudit=q1_compute("verify_moisture_trace",trial.trace_path,trial,c.mesh);
assert(heatAudit.steps==round(p.tFinal/dt)&&moistureAudit.steps==c.accepted_steps, ...
    "CUMCM:q1:BalanceFailure","Trace does not contain all accepted steps.");
assert(heatAudit.maximum_balance_J<=10*2*pi*p.length_m*p.rho*p.cp*sum(h.mesh.w)*p.algebraicTolerance(algebraicLevel+1) && ...
    moistureAudit.maximum_balance_kgkg_m2<=10*sum(c.mesh.w)*p.algebraicTolerance(algebraicLevel+1), ...
    "CUMCM:q1:BalanceFailure","Independent full-trace balance failed.");
fields=q1_compute("extract_output_grid",h,c,N);
report=struct("scope","q1_fixed_reproduction","status","REPRODUCED", ...
    "radial_segments",N,"max_step_s",dt,"algebraic_level",algebraicLevel, ...
    "quadrature_level",quadratureLevel,"time_steps",round(p.tFinal/dt), ...
    "input_audit",inputAudit,"heat_balance",heatAudit,"moisture_balance",moistureAudit, ...
    "temperature_range_C",[min(h.state,[],"all"),max(h.state,[],"all")], ...
    "moisture_range_kgkg",[min(c.state,[],"all"),max(c.state,[],"all")], ...
    "moisture_rejected_steps",c.rejected_steps,"runtime_seconds",toc(p.totalTimer), ...
    "public_command",cfg.publicReproductionCommand);
paths=q1_outputs("write_csv",cfg,fields,report);
paths=[paths,q1_outputs("write_result1",cfg,fields),q1_outputs("plot_results",cfg,fields,report)];
report.runtime_seconds=toc(p.totalTimer);
output=struct("artifacts",paths,"metrics",report);
fprintf('q1 reproduction complete in %.1f s.\n',report.runtime_seconds);
clear cleanup
end

function removeTraceDir(traceDir,traceParent)
actual=string(java.io.File(char(traceDir)).getCanonicalPath());
parent=string(java.io.File(char(traceParent)).getCanonicalPath());
[~,name]=fileparts(actual);
if startsWith(lower(actual),lower(parent+filesep))&&startsWith(name,"reproduce_")&&isfolder(actual)
    rmdir(actual,"s");
end
end

function records=codeRecords(cfg)
names=["main_q1.m","q1_compute.m","q1_validation.m","q1_outputs.m", ...
    "q1_trace_store.m"];
records=repmat(struct("path","","sha256","","bytes",0),numel(names),1);
for k=1:numel(names)
    pathValue=fullfile(cfg.projectRoot,names(k));
    fid=fopen(pathValue,"rb");
    assert(fid>=0,"CUMCM:q1:InputStructureMismatch","Cannot read code.");
    cleanup=onCleanup(@() fclose(fid)); bytes=fread(fid,Inf,"*uint8"); clear cleanup
    digest=java.security.MessageDigest.getInstance("SHA-256"); digest.update(typecast(bytes,"int8"));
    hash=lower(string(reshape(dec2hex(typecast(digest.digest(),"uint8"),2)',1,[])));
    records(k)=struct("path",names(k),"sha256",hash,"bytes",numel(bytes));
end
end
