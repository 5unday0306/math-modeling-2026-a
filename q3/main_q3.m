function result=main_q3(mode)
arguments
    mode (1,1) string = "reproduce"
end
assert(ismember(mode,["reproduce","finish","validation"]), ...
    'CUMCM:q3:InvalidMode','mode must be reproduce, finish or validation.');
root=string(fileparts(mfilename("fullpath")));
addpath(root);
cfg=q3_compute("config",root);
if string(mode)=="reproduce"
    result=run_q3_weighted(cfg);return
end
if string(mode)=="finish"
    result=q3_outputs("finish",cfg);return
end
loaded=load(fullfile(root,"附件","未舍入结果.mat"));
result=q3_validation("solution",cfg,loaded);
end

function result=run_q3_weighted(cfg)
root=cfg.projectRoot;
data=readmatrix(fullfile(cfg.rawDataDir,'附件1.xlsx'),'Sheet','Sheet1','Range','A2:C242');
assert(isequal(size(data),[241,3])&&isequal(data(:,1),(0:60:14400)')&&all(isfinite(data),'all'));
ambient=struct('time_s',data(:,1),'temperature_C',data(:,2),'moisture_kgkg',data(:,3));
selection=struct('temperature',q3_compute("select_window",data(:,1),data(:,2),50,.8,30), ...
    'moisture',q3_compute("select_window",data(:,1),data(:,3),.05,.8,30));
scenario=struct('id','weighted_variance_08','temperature_C',selection.temperature.time_mean, ...
    'moisture_kgkg',selection.moisture.time_mean,'switch_time_s',14400, ...
    'temperature_window_s',[selection.temperature.start_min*60,14400], ...
    'moisture_window_s',[selection.moisture.start_min*60,14400]);
fprintf('温度窗口 %.0f—240 min，%.0f点，时间均值 %.15g ℃。\n', ...
    selection.temperature.start_min,selection.temperature.sample_count,scenario.temperature_C);
fprintf('水分窗口 %.0f—240 min，%.0f点，时间均值 %.15g kg/kg。\n', ...
    selection.moisture.start_min,selection.moisture.sample_count,scenario.moisture_kgkg);
loaded=load(fullfile(root,'data','q3_start_4h.mat'),'upstream','prefix','earlyMaximum');
p=q3_compute("parameters",cfg,loaded.upstream);timer=tic;
stamp=string(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));caseRoot="计算缓存/"+stamp;
candidate=q3_compute("solve_event_case",cfg,loaded.upstream,p,ambient,scenario,caseRoot,timer,loaded.prefix);
history=q3_compute("tail_maximum_history",root,loaded.earlyMaximum,candidate);
computeSeconds=toc(timer);
save(fullfile(root,'附件','未舍入结果.mat'),'candidate','selection','history','ambient','computeSeconds','-v7');
result=q3_outputs("finish",cfg);
fprintf('\n加权选窗方案完成：%.4f小时，%d行；计算耗时%.1f秒。\n', ...
    result.stop_time_h,result.output_rows,result.compute_seconds);
end
