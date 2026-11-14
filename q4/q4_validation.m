function varargout=q4_validation(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q4:ValidationAction','Invalid validation action.');
switch action
    case "environment"
        varargout{1}=validate_q4_environment(varargin{:});
    case "answer"
        varargout{1}=validate_q4_answer(varargin{:});
    otherwise
        error('CUMCM:q4:ValidationAction','Unknown validation action: %s',action);
end
end

function checks=validate_q4_environment(in,p)
a=in.ambient;stats=in.environment_selection;cfg=p.selectionConfig;
names=["temperature","moisture"];expectedFirst=[131,151];expectedCount=[80,60];
expectedPass=[75,39];expectedStart=[136,195];expectedN=[105,46];
expectedMean=[49.994057142857145,0.049989782608695646];
errors=zeros(2,4);
for c=1:2
    q=stats.(names(c));w=q.windows;
    assert(q.selection.first_exceed_min==expectedFirst(c)&&height(w)==expectedCount(c));
    assert(q.selection.screen_pass_count==expectedPass(c)&&q.selection.start_min==expectedStart(c));
    assert(q.selection.sample_count==expectedN(c)&&abs(q.selection.mean_value-expectedMean(c))<1e-12);
    threshold=[cfg.temperatureThreshold,cfg.moistureThreshold];
    assert(a(expectedFirst(c)+1,c+1)>threshold(c)&&all(a(1:expectedFirst(c),c+1)<=threshold(c)));
    for k=1:height(w)
        start=w.start_min(k);x=a(start+1:end,c+1);n=numel(x);
        assert(n==w.sample_count(k)&&n==241-start&&w.end_s(k)==14400);
        mu=sum(x)/n;v=var(x,0);tm=(start:240)';fit=[ones(n,1),tm-mean(tm)]\x;
        edges=[0,cumsum(w.block_sizes(k,:))];
        bm=zeros(1,3);bv=bm;
        for b=1:3
            block=x(edges(b)+1:edges(b+1));bm(b)=mean(block);bv(b)=var(block,0);
        end
        d=abs(fit(2))*(240-start)/sqrt(v);span=(max(bm)-min(bm))/sqrt(v);ratio=max(bv)/min(bv);
        errors(c,:)=max(errors(c,:),[abs(mu-w.mean_value(k)),abs(v-w.sample_variance(k)), ...
            abs(d-w.standardized_drift(k)),max(abs(bv-w.block_variances(k,:)))]);
        assert(abs(span-w.three_block_mean_span(k))<1e-10&&abs(ratio-w.three_block_variance_ratio(k))<1e-10);
        assert(w.screen_pass(k)==(d<=cfg.driftMax&&span<=cfg.blockSpanMax&&ratio<=cfg.varianceRatioMax));
    end
    expectedLoss=cfg.lambda*log(w.sample_variance/min(w.sample_variance))+ ...
        (1-cfg.lambda)*log(max(w.sample_count)./w.sample_count);
    assert(max(abs(expectedLoss-w.total_loss))<1e-13);
    assert(max(abs(w.total_loss-w.variance_penalty-w.sample_penalty))<1e-13);
    assert(any(a(expectedStart(c)+1:end,c+1)<threshold(c)),'Below-threshold samples must remain.');
end
assert(all(errors(:,1)<1e-12)&&all(errors(:,2)<1e-14)&&all(errors(:,3)<1e-9)&&all(errors(:,4)<1e-14));
assert(isequal(stats.temperature.selection.block_sizes,[35,35,35]));
assert(isequal(stats.moisture.selection.block_sizes,[16,15,15]));
assert(isequal(stats.temperature.sensitivity.weight_selected_start_min,[131,131,131,136,140,140]));
assert(isequal(stats.moisture.sensitivity.weight_selected_start_min,[151,151,151,195,195,195]));
converted=a;converted(:,2)=a(:,2)+273.15;converted(:,3)=a(:,3)*1000;cc=cfg;
cc.temperatureThreshold=cfg.temperatureThreshold+273.15;cc.moistureThreshold=cfg.moistureThreshold*1000;
shift=q4_compute("select_suffix",converted,cc);
assert(shift.temperature.selection.start_min==136&&shift.moisture.selection.start_min==195);
assert(max(abs(shift.temperature.windows.total_loss-stats.temperature.windows.total_loss))<1e-10);
assert(max(abs(shift.moisture.windows.total_loss-stats.moisture.windows.total_loss))<1e-10);
bad=a;bad(:,2)=50;mustReject(bad,cfg,'CUMCM:q4:NoExceedance');
bad=a;bad(:,2)=50;bad(212:end,2)=51;mustReject(bad,cfg,'CUMCM:q4:NoCandidates');
bad=a;bad(:,2)=51;mustReject(bad,cfg,'CUMCM:q4:DegenerateVariance');
bad=a;bad(:,2)=51.1;mustReject(bad,cfg,'CUMCM:q4:DegenerateVariance');
bad=a;bad(2,2)=NaN;mustReject(bad,cfg,'CUMCM:q4:SelectionInput');
bad=a;bad(2,1)=0;mustReject(bad,cfg,'CUMCM:q4:SelectionInput');
bad=a;bad(:,2)=51+(0:240)'/240;mustReject(bad,cfg,'CUMCM:q4:NoStableWindow');
bad=a;bad(end-9:end,2)=50;
degenerate=q4_compute("select_suffix",bad,cfg);
assert(~degenerate.temperature.windows.screen_pass(end)&& ...
    contains(degenerate.temperature.windows.rejection_reason(end),'DEGENERATE_BLOCK'));
[chosen,tol]=q4_compute("choose_suffix",[1;1+8*eps(1)],[true;true],[31;90],64);
assert(chosen==2&&tol>=8*eps(1));
assert(q4_compute("choose_suffix",[1;1.1],[true;true],[31;90],64)==1);
assert(q4_compute("choose_suffix",[0;1],[false;true],[90;31],64)==2);
assert(isempty(q4_compute("choose_suffix",[0;1],[false;false],[90;31],64)));
assert(.8*log(2)>.2*log(2)&&.2*log(100/31)>.2*log(100/90));
scenarios=["unified_stable_suffix","temperature_last","temperature_mean30", ...
    "temperature_mean58","temperature_mean59","temperature_mean60"];
for scenario=scenarios
    for t=[0,60,14399,14400]
        bc=q4_compute("interpolate_ambient",t,in,scenario);expected=interp1(a(:,1),a(:,2:3),t,'linear');
        assert(abs(bc.temperature_K-(expected(1)+273.15))<1e-12&&bc.moisture_kgkg==expected(2));
    end
    bc=q4_compute("interpolate_ambient",14400+eps(14400),in,scenario);
    assert(bc.moisture_kgkg==stats.selection.moisture_kgkg);
end
left=q4_compute("interpolate_ambient",14400,in,p.environment);right=q4_compute("interpolate_ambient",14400+eps(14400),in,p.environment);
assert(right.temperature_K==stats.selection.temperature_K&&left.moisture_kgkg==a(end,3));
try
    q4_compute("interpolate_ambient",14401,in,"temperature_tail_variance");
    error('CUMCM:q4:ExpectedFailure','Obsolete scenario accepted.');
catch ex
    assert(strcmp(ex.identifier,'CUMCM:q4:Ambient'));
end
checks=struct('formula_errors_mean_variance_drift_blockvariance',errors, ...
    'temperature',stats.temperature.selection,'moisture',stats.moisture.selection, ...
    'temperature_jump_K',right.temperature_K-left.temperature_K, ...
    'moisture_jump_kgkg',right.moisture_kgkg-left.moisture_kgkg, ...
    'units_and_threshold_transformed_together',true,'negative_test_count',7, ...
    'degenerate_blocks_excluded',true,'roundoff_ties_tested',true,'all_checks_passed',true);
end

function mustReject(a,cfg,id)
try
    q4_compute("select_suffix",a,cfg);
    error('CUMCM:q4:ExpectedFailure','Invalid data were accepted.');
catch ex
    assert(strcmp(ex.identifier,id),'CUMCM:q4:UnexpectedFailure','Expected %s, got %s',id,ex.identifier);
end
end


function answer=validate_q4_answer(candidate,in,p,mode)
event=candidate.event;
assert(event.strict_nominal_stop&&event.value<p.threshold);
assert(event.left_s<event.root_s&&event.root_s<event.right_s);
assert(event.right_s-event.left_s<=p.rootWidth_s);
geo=q4_compute("interpolate_radius",event.stop_time_s,in,p.geometry);
bc=q4_compute("interpolate_ambient",event.stop_time_s,in,p.environment);
maximum=q4_compute("spatial_maximum",event.stop_state,candidate.mesh,geo,bc,p);
assert(maximum.maximum_C<p.threshold&&abs(maximum.maximum_C-event.value)<1e-12);
assert(candidate.accepted_half_steps>0&&candidate.rejected_macro_steps>=0);
answer=struct('stop_time_s',event.stop_time_s,'stop_time_h',event.stop_time_h, ...
    'root_s',event.root_s,'root_bracket_s',[event.left_s,event.right_s], ...
    'maximum_moisture_kgkg',maximum.maximum_C,'surface_moisture_kgkg',event.stop_state(end), ...
    'mean_first_crossing_h',candidate.mean_event.root_s/3600, ...
    'final_radius_cm',100*geo.R_m,'accepted_half_steps',candidate.accepted_half_steps, ...
    'rejected_macro_steps',candidate.rejected_macro_steps,'answer_checks_passed',true);
if string(mode)=="result"
    assert(p.N==160);
    assert(abs(answer.stop_time_h-51.1128)<1e-12);
    assert(abs(answer.root_s-184005.8642578125)<=p.rootWidth_s);
    assert(abs(answer.maximum_moisture_kgkg-0.14999988850250043)<5e-7);
    assert(abs(answer.surface_moisture_kgkg-0.052608260343871151)<5e-7);
    assert(abs(answer.mean_first_crossing_h-34.988049316406247)<5e-6);
    assert(abs(answer.final_radius_cm-1.2)<1e-12);
end
end
