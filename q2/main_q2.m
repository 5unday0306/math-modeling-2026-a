function output=main_q2(mode)
arguments
    mode (1,1) string = "reproduce"
end
assert(ismember(mode,["mvp","reproduce","validation"]),'CUMCM:q2:InvalidMode','mode must be mvp, reproduce or validation.');
projectRoot=string(fileparts(mfilename("fullpath")));
addpath(projectRoot);
cfg=q2_config(projectRoot);
output=run_q2(cfg,mode);
end

function cfg=q2_config(projectRoot)
arguments
    projectRoot (1,1) string
end
cfg.problem="A";
cfg.projectRoot=projectRoot;
cfg.rawDataDir=fullfile(projectRoot,"data","raw");
cfg.processedDataDir=fullfile(projectRoot,"data","processed");
cfg.resultsDir=fullfile(projectRoot,"results");
cfg.figuresDir=fullfile(projectRoot,"figures");
cfg.figureQADir=fullfile(projectRoot,"artifacts","validation","figure_previews");
cfg.validationDir=fullfile(projectRoot,"artifacts","validation");
cfg.randomSeed=42;
cfg.publicReproductionCommand="main_q2(""reproduce"")";
cfg.questionCommandTemplate="main_q2(""<mode>"")";
end

function output=run_q2(cfg,mode)
if string(mode)=="validation"
    p=q2_compute("parameters");[ambient,~]=q2_compute("load_inputs",cfg);
    output=q2_validation("resources",cfg,p,ambient);return
end
if string(mode)=="reproduce"
    output=reproduce_q2(cfg);return
end
assert(string(mode)=="mvp",'CUMCM:q2:NotImplemented','Final requires the stage31/freeze interface.');
started=tic;p=q2_compute("parameters");[ambient,inputAudit]=q2_compute("load_inputs",cfg);
mesh=q2_compute("build_mesh",p.radius_m,20);
stamp=datetime('now','TimeZone','UTC','Format',"yyyyMMdd'T'HHmmssSSS'Z'");
nonce=char(java.util.UUID.randomUUID());caseId="q2_mvp_trace_"+string(stamp)+"_"+string(nonce(1:8));
traceRelative="artifacts/q2/validation/"+caseId+"/accepted_states.h5";
tracePath=fullfile(cfg.projectRoot,traceRelative);
disk=java.io.File(char(cfg.projectRoot)).getUsableSpace();
assert(disk>8*(2*(mesh.N+1)+13)*(p.tFinal+1)+1024^3, ...
    'CUMCM:q2:ResourceLimit','Insufficient trace space plus 1GiB reserve.');
settings=struct('maxStep_s',1,'algebraicLevel',0,'quadratureLevel',0, ...
    'totalTimer',started,'trace_path',tracePath);
fprintf('q2 MVP: coupled Appendix3; N=20; dt=1 s; horizon=10800 s.\n');
solution=q2_compute("solve",p,ambient,mesh,settings);
fields=q2_compute("extract_output_grid",solution,mesh);
candidate=struct('mesh',mesh,'solution',solution,'fields',fields);
checks=q2_validation("solution",cfg,p,ambient,candidate,mode);checks.input_audit=inputAudit;
checks.code_snapshot=codeRecords(cfg);checks.trace_path=traceRelative;
values=zeros(30,4);row=0;
for t=fields.paper_time_s'
    for col=fields.paper_columns
        row=row+1;values(row,:)=[t,fields.radius_cm(col),fields.temperature_C(t,col),fields.moisture_kgkg(t,col)];
    end
end
folder=fullfile(cfg.projectRoot,'artifacts','q2','mvp');if ~isfolder(folder),mkdir(folder);end
csvPath=fullfile(folder,'q2_mvp_slice.csv');fid=fopen(csvPath,'w','n','UTF-8');assert(fid>=0);
c=onCleanup(@()fclose(fid));fprintf(fid,'时间（秒）,半径（厘米）,温度（℃）,干基含水率（kg/kg）\n');
fprintf(fid,'%.17g,%.17g,%.17g,%.17g\n',values');clear c
assert(isequal(readmatrix(csvPath,'NumHeaderLines',1),values),'CUMCM:q2:OutputMapping','CSV values differ on reread.');
header=readcell(csvPath,'Range','A1:D1','Encoding','UTF-8');
assert(isequal(string(header),["时间（秒）","半径（厘米）","温度（℃）","干基含水率（kg/kg）"]));
checks.slice_rows=30;checks.slice_values=values;
final_state=struct('time_s',solution.final_time_s,'radius_m',mesh.r, ...
    'temperature_K',solution.final_state(1:2:end),'moisture_kgkg',solution.final_state(2:2:end), ...
    'purpose','MVP_ONLY_NOT_VALIDATED_Q3_RESTART');
save(fullfile(folder,'q2_mvp_final_state.mat'),'final_state','-v7');
checks.runtime_seconds=toc(started);assert(checks.runtime_seconds<=p.wallSeconds,'CUMCM:q2:ResourceLimit','MVP total time exceeded.');
fid=fopen(fullfile(folder,'q2_mvp_checks.json'),'w','n','UTF-8');assert(fid>=0);
c=onCleanup(@()fclose(fid));fprintf(fid,'%s\n',jsonencode(checks,'PrettyPrint',true));clear c
output=struct('artifacts',["artifacts/q2/mvp/q2_mvp_slice.csv", ...
    "artifacts/q2/mvp/q2_mvp_checks.json","artifacts/q2/mvp/q2_mvp_final_state.mat",traceRelative], ...
    'metrics',checks);
disp(values);fprintf('q2 MVP complete: %d accepted steps, %d rejected; no final precision claim.\n',solution.accepted_steps,solution.rejected_steps);
end

function output=reproduce_q2(cfg)
started=tic;p=q2_compute("parameters");[ambient,inputAudit]=q2_compute("load_inputs",cfg);
N=160;dt=0.25;algebraicLevel=0;quadratureLevel=0;
mesh=q2_compute("build_mesh",p.radius_m,N);
traceParent=fullfile(cfg.projectRoot,"artifacts","q2","reproduce_traces");
if ~isfolder(traceParent),mkdir(traceParent);end
nonce=char(java.util.UUID.randomUUID());runId="q2_reproduce_"+string(nonce(1:8));
traceDir=fullfile(traceParent,runId);mkdir(traceDir);
cleanup=onCleanup(@()removeTraceDir(traceDir,traceParent));
settings=struct('maxStep_s',dt,'algebraicLevel',algebraicLevel, ...
    'quadratureLevel',quadratureLevel,'totalTimer',started, ...
    'trace_path',fullfile(traceDir,'accepted_states.h5'));
fprintf('q2 reproduction: N=%d; dt=%.2f s; horizon=%d s.\n',N,dt,p.tFinal);
solution=q2_compute("solve",p,ambient,mesh,settings);
fields=q2_compute("extract_output_grid",solution,mesh);
traceAudit=q2_compute("verify_trace",solution.trace,p,mesh,ambient);
balance=q2_compute("balance_summary",solution.trace,p,mesh,ambient);
report=struct('scope','q2_fixed_reproduction','run_id',runId,'status','REPRODUCED', ...
    'radial_segments',N,'max_step_s',dt,'algebraic_level',algebraicLevel, ...
    'quadrature_level',quadratureLevel,'time_steps',round(p.tFinal/dt), ...
    'accepted_steps',solution.accepted_steps,'rejected_steps',solution.rejected_steps, ...
    'input_audit',inputAudit,'trace_audit',traceAudit,'balance',balance, ...
    'temperature_range_C',[min(fields.temperature_C,[],"all"),max(fields.temperature_C,[],"all")], ...
    'moisture_range_kgkg',[min(fields.moisture_kgkg,[],"all"),max(fields.moisture_kgkg,[],"all")], ...
    'files_written',false,'public_command',cfg.publicReproductionCommand, ...
    'runtime_seconds',toc(started));
output=struct('artifacts',strings(1,0),'fields',fields,'metrics',report);
fprintf('q2 reproduction complete in %.1f s.\n',report.runtime_seconds);
clear cleanup
end

function removeTraceDir(traceDir,traceParent)
actual=string(java.io.File(char(traceDir)).getCanonicalPath());
parent=string(java.io.File(char(traceParent)).getCanonicalPath());
[~,name]=fileparts(actual);
if startsWith(lower(actual),lower(parent+filesep))&&startsWith(name,"q2_reproduce_")&&isfolder(actual)
    rmdir(actual,'s');
end
end

function records=codeRecords(cfg)
names=["main_q2.m","q2_compute.m","q2_validation.m","q2_outputs.m","q2_trace_store.m"];
records=repmat(struct('path','','sha256','','bytes',0),numel(names),1);
for k=1:numel(names)
    pathValue=fullfile(cfg.projectRoot,names(k));
    bytes=java.nio.file.Files.readAllBytes(java.io.File(char(pathValue)).toPath());
    md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);
    hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
    info=dir(pathValue);
    records(k)=struct('path',names(k),'sha256',hash,'bytes',info.bytes);
end
end
