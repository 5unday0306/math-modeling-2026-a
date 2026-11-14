function output=main_q4(mode)
arguments
    mode (1,1) string = "result"
end
assert(ismember(mode,["result","check","export"]),'CUMCM:q4:Mode','mode must be result, check or export.');
projectRoot=string(fileparts(mfilename("fullpath")));
addpath(projectRoot);
cfg=q4_config(projectRoot);
rng(cfg.randomSeed,'twister');
in=q4_compute("load_inputs",cfg);
if mode=="check",caseSpec=struct('N',20);else,caseSpec=struct('N',160);end
p=q4_compute("parameters",cfg,caseSpec);
p.wallSeconds=3600;
environment=q4_validation("environment",in,p);
traceParent=fullfile(projectRoot,'artifacts','q4','runtime');
if ~isfolder(traceParent),mkdir(traceParent);end
nonce=char(java.util.UUID.randomUUID());runId="q4_"+mode+"_"+string(nonce(1:8));
traceDir=fullfile(traceParent,runId);
relative="artifacts/q4/runtime/"+runId;
cleanup=onCleanup(@()removeRuntime(traceDir,traceParent));
timer=tic;
candidate=q4_compute("solve_case",cfg,in,p,relative,timer);
answerMode="result";if mode=="check",answerMode="check";end
answer=q4_validation("answer",candidate,in,p,answerMode);
paths=strings(1,0);filesWritten=false;
if mode=="export"
    summary=struct('primary_only',true,'status','NUMERICAL_CANDIDATE', ...
        'precision_verified',false,'case_name','N160','candidate_N',160, ...
        'candidate_stop_h',answer.stop_time_h);
    paths=q4_outputs("write_all",cfg,candidate,summary,in);
    filesWritten=true;
end
candidate=rmfield(candidate,'trace_paths');
output=struct('answer',answer,'candidate',candidate,'environment_checks',environment, ...
    'artifacts',paths,'files_written',filesWritten,'runtime_seconds',toc(timer));
fprintf('q4 %s: %.4f h; max C %.12g; accepted %d; rejected %d.\n', ...
    mode,answer.stop_time_h,answer.maximum_moisture_kgkg, ...
    answer.accepted_half_steps,answer.rejected_macro_steps);
clear cleanup
end

function cfg=q4_config(projectRoot)
cfg=struct('projectRoot',projectRoot,'rawDataDir',fullfile(projectRoot,'data','raw'), ...
    'resultsDir',fullfile(projectRoot,'results'),'processedDataDir',fullfile(projectRoot,'data','processed'), ...
    'figuresDir',fullfile(projectRoot,'figures'), ...
    'figureQADir',fullfile(projectRoot,'artifacts','validation','figure_previews'), ...
    'randomSeed',42);
end

function removeRuntime(pathValue,parentValue)
actual=string(java.io.File(char(pathValue)).getCanonicalPath());
parent=string(java.io.File(char(parentValue)).getCanonicalPath());
[~,name]=fileparts(actual);
if startsWith(lower(actual),lower(parent+filesep))&&startsWith(name,"q4_")&&isfolder(actual)
    rmdir(actual,'s');
end
end
