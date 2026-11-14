function varargout=q3_validation(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q3:InvalidValidationAction','Invalid validation action.');
switch action
    case "solution"
        varargout{1}=validate_q3_solution(varargin{:});
    otherwise
        error('CUMCM:q3:InvalidValidationAction','Unknown validation action: %s',action);
end
end

function report=validate_q3_solution(cfg,loaded)
root=cfg.projectRoot;
candidate=loaded.candidate;event=candidate.event;fields=candidate.fields;
strict=event.strict_nominal_stop&&event.maximum_C<0.15;
bracket=event.bracket_left_s<event.stop_time_s&&event.stop_time_s<event.bracket_right_s;
endMs=round(event.stop_time_s*1000);
lattice=endMs<=flintmax&&abs(endMs/1000-event.stop_time_s)<=32*eps(max(1,event.stop_time_s));
fieldsOk=isequal(size(fields.moisture_kgkg),[numel(fields.time_s),21])&& ...
    all(isfinite(fields.moisture_kgkg)&fields.moisture_kgkg>0,'all');
tracesOk=true;traceColumns=zeros(numel(candidate.trace_paths),1);
for k=1:numel(candidate.trace_paths)
    path=fullfile(root,candidate.trace_paths(k));
    tracesOk=tracesOk&&isfile(path)&&h5readatt(path,'/accepted','complete')==1;
    info=h5info(path,'/accepted');traceColumns(k)=info.Dataspace.Size(2);
end
historyMonotone=all(diff(loaded.history(:,1))>0)&&all(isfinite(loaded.history),'all');
selectionOk=abs(loaded.selection.temperature.alpha-0.8)<eps && ...
    abs(loaded.selection.moisture.alpha-0.8)<eps;
status='PASS_SCREEN';
if ~(strict&&bracket&&lattice&&fieldsOk&&tracesOk&&historyMonotone&&selectionOk)
    status='FAIL_SCREEN';
end
report=struct('strict_nominal_stop',strict,'bracket_orientation',bracket, ...
    'stop_on_millisecond_lattice',lattice,'fields_finite_positive',fieldsOk, ...
    'traces_complete',tracesOk,'trace_columns',traceColumns, ...
    'maximum_history_monotone',historyMonotone,'selection_alpha_08',selectionOk, ...
    'stop_time_h',event.stop_time_h,'stop_time_s',event.stop_time_s, ...
    'maximum_moisture_at_stop',event.maximum_C,'event_status',event.status, ...
    'output_rows',numel(fields.time_s),'status',status);
end
