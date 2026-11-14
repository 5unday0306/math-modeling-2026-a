function varargout=q2_validation(action,varargin)
action=string(action);
assert(isscalar(action),'CUMCM:q2:ValidationAction','Invalid validation action.');
switch action
    case "resources"
        varargout{1}=validate_q2_resources(varargin{:});
    case "solution"
        varargout{1}=validate_q2_solution(varargin{:});
    otherwise
        error('CUMCM:q2:ValidationAction','Unknown validation action: %s',action);
end
end

function output=validate_q2_resources(cfg,p,ambient)
started=tic;p.wallSeconds=3480;p.validationTimer=started;
stamp=datetime('now','TimeZone','UTC','Format',"yyyyMMdd'T'HHmmssSSS'Z'");
nonce=char(java.util.UUID.randomUUID());runId="q2_validation_cases_"+string(stamp)+"_"+string(nonce(1:8));
relative="artifacts/q2/validation/"+runId;folder=fullfile(cfg.projectRoot,relative);mkdir(folder);
cache=containers.Map('KeyType','char','ValueType','any');cases=struct([]);traceBytes=0;
priorCharge=0;remainingWall=3600;recovered=[];
report=struct('run_id',runId,'status','RUNNING','precision_complete',false,'stage31_status','NOT_RUN', ...
    'four_component_budget',1.25e-5,'joint_budget',5e-5,'wall_budget_seconds',3600, ...
    'trace_budget_bytes',p.maxTraceBytes,'working_budget_bytes',p.maxWorkingBytes, ...
    'q3_production_restart_authorized',false,'pending_items',strings(0,1));
report.input_bindings=[q2_outputs("validation_io",'record',cfg.projectRoot,"data/raw/附件1.xlsx"), ...
    q2_outputs("validation_io",'record',cfg.projectRoot,"data/raw/附件3/result2.xlsx"), ...
    q2_outputs("validation_io",'record',cfg.projectRoot,"models/q2/model_contract.md"), ...
    q2_outputs("validation_io",'record',cfg.projectRoot,"plan/contracts/q2_implementation_contract.md")];
source=["main_q2.m","q2_compute.m","q2_validation.m","q2_outputs.m","q2_trace_store.m"];
report.code_bindings=struct([]);
for k=1:numel(source)
    item=q2_outputs("validation_io",'record',cfg.projectRoot,source(k));
    if k==1,report.code_bindings=item;else,report.code_bindings(k)=item;end
end
candidate=[];lastError=[];rounds=struct([]);
try
    statePath=fullfile(cfg.projectRoot,'plan','workflow_state.json');
    workflow=jsondecode(fileread(statePath));
    if isfield(workflow.questions.q2,'validation_resume')&&workflow.questions.q2.validation_resume.enabled
        recovered=q2_recover_validation(cfg,p,ambient,workflow.questions.q2.validation_resume);
        cache=recovered.cache;cases=recovered.cases;traceBytes=recovered.trace_bytes;
        priorCharge=recovered.prior_charge_seconds;remainingWall=3600-priorCharge;
        p.wallSeconds=remainingWall-120;
        report.recovery=recovered.audit;
        report.input_bindings=[report.input_bindings(:);recovered.bindings(:)]';
        fprintf('q2 validation: recovered %d completed cases; prior budget charged %.0f s.\n',numel(cases),priorCharge);
    end
    base=runCase(20,1,0,0,ambient,'linear');candidate=base;
    if isempty(recovered)
    report.mvp_identities=validate_q2_solution(cfg,p,ambient,base,'mvp');checkpoint();
    fprintf('q2 validation: frozen-coefficient Bessel reference.\n');
    report.baseline=q2_validation_baselines(p,ambient,base.fields);checkpoint();
    fprintf('q2 validation: pdepe N40/N80 self-refinement.\n');
    cross=validate_q2_pdepe(p,ambient,base.fields);crossFields=cross.fine_fields;
    report.pdepe=rmfield(cross,'fine_fields');clear cross;checkpoint();
    fprintf('q2 validation: finite-cylinder N20x20/N40x40 sensitivity.\n');
    report.finite_cylinder=validate_q2_finite_cylinder(p,ambient,base.fields);checkpoint();
    else
        for name=["mvp_identities","baseline","pdepe","finite_cylinder"]
            assert(isfield(recovered.old_report,name),'CUMCM:q2:RecoveryInput','Missing completed auxiliary evidence.');
            report.(name)=recovered.old_report.(name);
        end
        checkpoint();
    end
    alternative=ambient;ix=1:181;
    alternative.time_s=ambient.time_s(ix);alternative.temperature_C=ambient.temperature_C(ix);
    alternative.moisture_kgkg=ambient.moisture_kgkg(ix);
    alternative.validation_interpolation='PCHIP_SENSITIVITY_ONLY';
    alternative.temperature_pp=pchip(alternative.time_s,alternative.temperature_C);
    alternative.moisture_pp=pchip(alternative.time_s,alternative.moisture_kgkg);
    alt=runCase(20,1,0,0,alternative,'pchip');
    report.interpolation_sensitivity=q2_compute("compare_fields",base.fields,alt.fields);
    report.interpolation_sensitivity.scope='Same N20/dt1; first181 observations only; not a discretization error.';
    report.interpolation_sensitivity.post_3h_extrapolation_used=false;checkpoint();
    N=20;dt=1;roundIndex=0;
    while true
        roundIndex=roundIndex+1;q2_validation_deadline(p);
        fprintf('q2 validation: refinement round %d, N=%d, dt=%.9g s.\n',roundIndex,N,dt);
        anchor=runCase(N,dt,0,0,ambient,'linear');candidate=anchor;
        space=runCase(2*N,dt,0,0,ambient,'linear');
        temporal=runCase(N,dt/2,0,0,ambient,'linear');
        algebra=runCase(N,dt,1,0,ambient,'linear');
        quadrature=runCase(N,dt,0,1,ambient,'linear');
        rr=struct('radial_segments',N,'max_step_s',dt, ...
            'space',q2_compute("compare_fields",anchor.fields,space.fields), ...
            'time',q2_compute("compare_fields",anchor.fields,temporal.fields), ...
            'algebraic',q2_compute("compare_fields",anchor.fields,algebra.fields), ...
            'quadrature',q2_compute("compare_fields",anchor.fields,quadrature.fields));
        errorMatrix=[rr.space.maximum_difference;rr.time.maximum_difference; ...
            rr.algebraic.maximum_difference;rr.quadrature.maximum_difference];
        rr.four_component_errors=errorMatrix;rr.components_pass=all(errorMatrix<1.25e-5,'all');
        rr.joint=struct();rr.joint_completed=false;
        if roundIndex==1,rounds=rr;else,rounds(roundIndex)=rr;end
        checkpoint();
        disp(errorMatrix);
        if rr.components_pass
            joint=runCase(2*N,dt/2,1,1,ambient,'linear');candidate=joint;
            rr.joint=q2_compute("compare_fields",anchor.fields,joint.fields);rr.joint_completed=true;
            rounds(roundIndex)=rr;checkpoint();
            if all(rr.joint.maximum_difference<5e-5)&&rr.joint.rounded_equal
                report.precision_complete=true;break
            end
        end
        refineSpace=any(rr.space.maximum_difference>=1.25e-5);
        refineTime=any(rr.time.maximum_difference>=1.25e-5);
        if rr.components_pass,refineSpace=true;refineTime=true;end
        if ~refineSpace&&~refineTime
            error('CUMCM:q2:PrecisionUnresolved','Algebraic/quadrature error still exceeds the approved budget.');
        end
        N=N*2^refineSpace;dt=dt/2^refineTime;
    end
    if ~exist('crossFields','var')
        cross=validate_q2_pdepe(p,ambient,base.fields);crossFields=cross.fine_fields;
        report.pdepe.selected_cross_fields_recomputed_after_recovery=true;
    end
    report.pdepe.fine_vs_selected=q2_compute("compare_fields",candidate.fields,crossFields);
catch e
    lastError=e;report.status='BLOCKED';report.error_id=e.identifier;
    report.error_message=replace(string(e.message),cfg.projectRoot,'.');
    report.pending_items=["Complete four-component and joint-refinement/rounding confirmation"; ...
        "Independent stage31 review; q3 production handoff remains disabled"];
    if ~contains(e.identifier,'ResourceLimit')&&~contains(e.identifier,'PrecisionUnresolved')
        report.status='FAIL';
    end
end
p.wallSeconds=remainingWall;
report.cases=cases;report.refinement_rounds=rounds;
report.runtime_seconds_before_output=toc(started);report.total_trace_bytes=traceBytes;
if ~isempty(candidate)
    report.selected_case=candidate.case_id;report.selected_settings=candidate.settings;
    report.selected_trace=candidate.trace_relative;
    report.selected_trace_audit=q2_compute("verify_trace",candidate.solution.trace,p,candidate.mesh,ambient);
    report.selected_balance=q2_compute("balance_summary",candidate.solution.trace,p,candidate.mesh,ambient);
    report.output_state='DRAFT_STAGE30_NOT_APPROVED_FOR_Q3';
    fields=candidate.fields;solution=candidate.solution;
    save(fullfile(folder,'selected_candidate.mat'),'fields','solution','-v7.3');
    mechanism=q2_compute("material_properties",fields.moisture_kgkg(:,[1,end]), ...
        fields.temperature_C(:,[1,end])+273.15,p);
    mechanism.temperature_C=fields.temperature_C(:,[1,end]);
    mechanism.moisture_kgkg=fields.moisture_kgkg(:,[1,end]);
    mechanism.logD_C=-p.a./mechanism.moisture_kgkg+p.a/p.C0;
    mechanism.logD_T=-p.B./(mechanism.temperature_C+273.15)+p.B/p.T0_K;
    initial=q2_compute("material_properties",p.C0,p.T0_K,p);
    report.mechanism_additivity_error=max(abs(log(mechanism.D/initial.D)-mechanism.logD_C-mechanism.logD_T),[],'all');
    assert(report.mechanism_additivity_error<1e-12,'CUMCM:q2:Mechanism','logD decomposition differs.');
    paths=q2_outputs("write_csv",cfg,fields,mechanism,report);
    paths(end+1)=q2_outputs("write_result2",cfg,fields);
    paths=[paths,q2_outputs("plot_results",cfg,fields,mechanism,report)];
    if report.precision_complete
        paths=[paths,q2_outputs("write_restart",cfg,candidate.solution, ...
            struct('run_id',runId,'input_bindings',report.input_bindings,'numerical_settings',candidate.settings))];
    end
    report.generated_artifacts=paths;report.full_field_values=453600;report.paper_values=60;
    report.official_template_copy_written=true;report.three_candidate_figures_written=true;
end
report.runtime_seconds=toc(started);report.prior_budget_charge_seconds=priorCharge;
assert(report.runtime_seconds+priorCharge<=3600,'CUMCM:q2:ResourceLimit','Output exceeded the shared validation budget.');
if report.precision_complete,report.status='AUTHOR_CHECKS_COMPLETE_AWAITING_STAGE31';end
q2_outputs("validation_io",'json',fullfile(folder,'validation_diagnostics.json'),report);
q2_outputs("validation_io",'json',fullfile(cfg.resultsDir,'q2','q2_validation_summary.json'),report);
if ~isempty(lastError),rethrow(lastError);end
q2_validation_deadline(p);
output=struct('artifacts',[paths,relative+"/validation_diagnostics.json",relative+"/selected_candidate.mat"], ...
    'metrics',report);
    function checkpoint()
        report.cases=cases;report.refinement_rounds=rounds;report.elapsed_seconds=toc(started);
        q2_outputs("validation_io",'json',fullfile(folder,'validation_diagnostics.json'),report);
    end
    function c=runCase(segments,step,algebraLevel,quadLevel,environment,tag)
        key=sprintf('N%d_dt%.12g_a%d_q%d_%s',segments,step,algebraLevel,quadLevel,tag);
        if isKey(cache,key),c=cache(key);return;end
        mesh=q2_compute("build_mesh",p.radius_m,segments);width=2*(segments+1)+13;
        expectedTrace=8*width*(ceil(p.tFinal/step)+1)+2^20;
        working=2*8*(8*(p.tFinal+1)*42+256*width+100*width)+cache.Count*8*42*(p.tFinal+1);
        disk=java.io.File(char(cfg.projectRoot)).getUsableSpace();predicted=25;
        if ~isempty(cases)
            weights=arrayfun(@(v)v.runtime_seconds/max(1,v.accepted_steps)/max(1,(v.radial_segments/20)^0.7),cases);
            predicted=1.35*max(weights)*ceil(p.tFinal/step)*(segments/20)^0.7;
        end
        reserve=120;
        plan=struct('case_id',key,'predicted_seconds',predicted,'remaining_seconds',remainingWall-toc(started), ...
            'planned_working_bytes',working,'planned_trace_bytes',expectedTrace,'existing_trace_bytes',traceBytes);
        report.next_case_resource_plan=plan;checkpoint();
        assert(toc(started)+predicted+reserve<=remainingWall&&working<=p.maxWorkingBytes&& ...
            traceBytes+expectedTrace<=p.maxTraceBytes&&disk>=expectedTrace+1024^3, ...
            'CUMCM:q2:ResourceLimit','Cannot start %s within total time/memory/trace budget. See resource plan.',key);
        traceRel=relative+"/"+string(key)+"/accepted_states.h5";
        settings=struct('maxStep_s',step,'algebraicLevel',algebraLevel,'quadratureLevel',quadLevel, ...
            'totalTimer',started,'trace_path',fullfile(cfg.projectRoot,traceRel));
        caseTimer=tic;fprintf('q2 case %s (estimated %.1f s).\n',key,predicted);
        solution=q2_compute("solve",p,environment,mesh,settings);
        fields=q2_compute("extract_output_grid",solution,mesh);
        traceAudit=q2_compute("verify_trace",solution.trace,p,mesh,environment);
        secondsUsed=toc(caseTimer);info=dir(settings.trace_path);traceBytes=traceBytes+info.bytes;
        meta=struct('case_id',key,'radial_segments',segments,'max_step_s',step, ...
            'algebraic_level',algebraLevel,'quadrature_level',quadLevel,'boundary',tag, ...
            'runtime_seconds',secondsUsed,'accepted_steps',solution.accepted_steps, ...
            'rejected_steps',solution.rejected_steps,'trace_path',traceRel,'trace_bytes',info.bytes, ...
            'trace_audit',traceAudit);
        if isempty(cases),cases=meta;else,cases(end+1)=meta;end
        c=struct('mesh',mesh,'solution',solution,'fields',fields,'case_id',key, ...
            'trace_relative',traceRel,'settings',rmfield(settings,{'totalTimer','trace_path'}));
        cache(key)=c;checkpoint();
    end
end

function report=validate_q2_solution(cfg,p,ambient,candidate,mode)
assert(string(mode)=="mvp",'CUMCM:q2:NotImplemented','Full validation is stage30.');
mesh=candidate.mesh; sol=candidate.solution; fields=candidate.fields;
assert(isequal(size(fields.temperature_C),[10800,21]) && ...
    isequal(size(fields.moisture_kgkg),[10800,21]) && fields.time_s(end)==10800);
assert(all(isfinite(fields.temperature_C),'all')&&all(fields.moisture_kgkg>0,'all'));
assert(sol.final_time_s==10800&&sol.rejected_steps>=0);
for N=[20,40,80,160,320,640,1280,2560,5120,10240]
    q2_compute("build_mesh",p.radius_m,N);
end
cl=[0.03;0.3;1.8];cr=[0.12;2.1;2.55];tl=[302;305;312];tr=[310;320;322];
f=q2_compute("face_coefficients",cl,cr,tl,tr,p,0);sw=q2_compute("face_coefficients",cr,cl,tr,tl,p,0);
reference=zeros(3,2);
for i=1:3
    C=@(s)(1-s)*cl(i)+s*cr(i); T=@(s)(1-s)*tl(i)+s*tr(i);
    reference(i,1)=integral(@(s)p.k0+p.kC*C(s)./(1+C(s)),0,1,'AbsTol',1e-14,'RelTol',1e-12);
    reference(i,2)=integral(@(s)p.D0*exp(-p.a./C(s)-p.B./T(s)),0,1,'AbsTol',1e-22,'RelTol',1e-12);
end
faceError=max(abs([f.k,f.D]-reference)./reference,[],'all');
swapError=max(abs([f.k-sw.k,f.D-sw.D])./reference,[],'all');
iso=q2_compute("face_coefficients",cl,cr,310*ones(3,1),310*ones(3,1),p,0);
isoReference=zeros(3,1);
for i=1:3
    isoReference(i)=integral(@(c)p.D0*exp(-p.a./c-p.B/310),cl(i),cr(i), ...
        'AbsTol',1e-22,'RelTol',1e-12)/(cr(i)-cl(i));
end
isoError=max(abs(iso.D-isoReference)./isoReference);
assert(max([faceError,swapError,isoError])<1e-9,'CUMCM:q2:FaceCheck','Face identity failed.');
jacError=zeros(2,3);m=mesh.N+1;d=2*m;ti=(1:2:d)';ci=ti+1;
for ccase=1:2
    y=zeros(d,1);y(ti)=303+12*mesh.r/mesh.radius;
    if ccase==1,y(ci)=2.4-1.8*mesh.r/mesh.radius;else,y(ci)=0.16-0.13*mesh.r/mesh.radius;end
    prev=y;prev(ti)=prev(ti)-0.03;prev(ci)=prev(ci)+0.001;
    boundary=struct('temperature_K',318,'moisture_kgkg',0.04);
    [~,J,a]=q2_compute("residual_jacobian",y,prev,mesh,1,boundary,p,0);
    for direction=1:3
        v=sin((1:d)'*direction);v(ci)=0.02*v(ci);h=1e-4;
        rp=q2_compute("residual_jacobian",y+h*v,prev,mesh,1,boundary,p,0);
        rm=q2_compute("residual_jacobian",y-h*v,prev,mesh,1,boundary,p,0);
        numerical=a.row_scale.*((rp-rm)/(2*h));exact=a.row_scale.*(J*v);
        jacError(ccase,direction)=max(abs(numerical-exact))/max(max(abs(exact)),eps);
    end
end
assert(max(jacError,[],'all')<2e-6,'CUMCM:q2:JacobianCheck','Jacobian derivative check failed.');
testC=[0.03;0.15;1;2.55];testT=[301.15;310;315;323.396];h=1e-6;
props=q2_compute("material_properties",testC,testT,p);
plus=q2_compute("material_properties",testC+h,testT,p);minus=q2_compute("material_properties",testC-h,testT,p);
deriv=abs((plus.b-minus.b)/(2*h)-props.b_C)./props.b_C;
assert(max(deriv)<1e-6,'CUMCM:q2:Derivative','Capacity derivative differs.');
short=p;short.tFinal=4;settings=struct('maxStep_s',1,'algebraicLevel',1,'quadratureLevel',0);
zero=short;zero.hc=0;zero.hm=0;
z=q2_compute("solve",zero,ambient,mesh,settings);
equal=ambient;equal.temperature_C(:)=p.T0_K-273.15;equal.moisture_kgkg(:)=p.C0;
steady=q2_compute("solve",short,equal,mesh,settings);
zeroError=[max(abs(z.temperature_K-p.T0_K),[],'all'),max(abs(z.moisture_kgkg-p.C0),[],'all')];
steadyError=[max(abs(steady.temperature_K-p.T0_K),[],'all'),max(abs(steady.moisture_kgkg-p.C0),[],'all')];
assert(max([zeroError,steadyError])<1e-12,'CUMCM:q2:LimitCheck','Uniform limit drifted.');
domainRejected=false;
try
    q2_compute("material_properties",0,301.15,p);
catch e
    domainRejected=strcmp(e.identifier,'CUMCM:q2:MaterialDomain');
end
assert(domainRejected,'CUMCM:q2:NegativeTest','Nonpositive C accepted.');
[endpointT,endpointC]=q2_compute("interpolate_ambient",[14400,14401],ambient);
assert(isequal(endpointT,[50.165,50.165])&&isequal(endpointC,[0.04986,0.04986]));
changed=ambient;changed.temperature_C(182:end)=999;changed.moisture_kgkg(182:end)=999;
times=(0:10800)';[aT,aC]=q2_compute("interpolate_ambient",times,ambient);[bT,bC]=q2_compute("interpolate_ambient",times,changed);
assert(isequal(aT,bT)&&isequal(aC,bC),'CUMCM:q2:BoundaryLeak','Post-3h changed prefix.');
traceAudit=q2_compute("verify_trace",sol.trace,p,mesh,ambient);
report=struct('mode','mvp','radial_segments',mesh.N,'time_horizon_s',10800, ...
    'full_field_values',numel(fields.temperature_C)+numel(fields.moisture_kgkg), ...
    'paper_positions',30,'paper_values',60,'output_interval_s',1, ...
    'temperature_range_C',[min(fields.temperature_C,[],'all'),max(fields.temperature_C,[],'all')], ...
    'moisture_range_kgkg',[min(fields.moisture_kgkg,[],'all'),max(fields.moisture_kgkg,[],'all')], ...
    'temperature_center_10800_C',fields.temperature_C(end,1), ...
    'temperature_surface_10800_C',fields.temperature_C(end,end), ...
    'moisture_center_10800_kgkg',fields.moisture_kgkg(end,1), ...
    'moisture_surface_10800_kgkg',fields.moisture_kgkg(end,end), ...
    'face_scalar_relative_error',faceError,'face_swap_relative_error',swapError, ...
    'isothermal_face_relative_error',isoError,'jacobian_direction_relative_errors',jacError, ...
    'capacity_derivative_relative_error',max(deriv),'zero_exchange_error',zeroError, ...
    'steady_error',steadyError,'nonpositive_state_rejected',domainRejected, ...
    'post_3h_input_exclusion',true,'trace_audit',traceAudit, ...
    'accepted_steps',sol.accepted_steps,'rejected_steps',sol.rejected_steps, ...
    'rejected_intervals',sol.rejected_intervals,'newton_max_iterations',sol.newton_max_iterations, ...
    'backtracks_max',sol.backtracks_max,'dt_min_s',sol.dt_min_s,'dt_max_s',sol.dt_max_s, ...
    'max_correction',sol.max_correction,'max_scaled_residual',sol.max_scaled_residual, ...
    'four_decimal_stability_verified',false,'q3_production_restart_authorized',false, ...
    'figure_generation','NOT_APPLICABLE_IN_MVP','official_xlsx_written',false, ...
    'optimization_or_statistical_test','NOT_APPLICABLE');
assert(isfolder(cfg.projectRoot));
end

function report=q2_validation_baselines(p,~,~)
prop=q2_compute("material_properties",p.C0,p.T0_K,p); outside=[323.15,0.04];
initial=[p.T0_K,p.C0]; diffusivity=[prop.alpha,prop.D];
conductivity=[prop.k,prop.D]; transfer=[p.hc,p.hm];
times=[1,10,60,300,600]';radius=linspace(0,p.radius_m,21);
exact=zeros(numel(times),21,2);truncation=zeros(1,2);errors=zeros(2,2);
for field=1:2
    bi=transfer(field)*p.radius_m/conductivity(field);
    terms=80;
    while true
        q2_validation_deadline(p);terms=2*terms;roots=besselRoots(bi,terms,p);
        coeff=2*besselj(1,roots)./(roots.*(besselj(0,roots).^2+besselj(1,roots).^2));
        modes=besselj(0,roots*(radius/p.radius_m));
        factors=exp(-diffusivity(field)*times*(roots'.^2)/p.radius_m^2);
        exact(:,:,field)=outside(field)+(initial(field)-outside(field))*((factors.*coeff')*modes);
        half=1:terms/2;
        shorter=outside(field)+(initial(field)-outside(field))* ...
            ((factors(:,half).*coeff(half)')*modes(half,:));
        truncation(field)=max(abs(exact(:,:,field)-shorter),[],'all');
        if truncation(field)<1e-9,break;end
        assert(terms*21*8*4<p.maxWorkingBytes,'CUMCM:q2:ResourceLimit','Bessel reference memory budget.');
    end
    seriesTerms(field)=terms;
    for level=1:2
        mesh=q2_compute("build_mesh",p.radius_m,20*2^(level-1));dt=2^(1-level);n=mesh.N+1;
        conduct=diffusivity(field)*mesh.faces/mesh.dr;
        K=spdiags([0;conduct]+[conduct;0],0,n,n)+ ...
            sparse(1:n-1,2:n,-conduct,n,n)+sparse(2:n,1:n-1,-conduct,n,n);
        beta=diffusivity(field)*transfer(field)/conductivity(field);
        K(n,n)=K(n,n)+mesh.radius*beta;
        rhs=zeros(n,1);rhs(n)=mesh.radius*beta*outside(field);
        M=spdiags(mesh.w/dt,0,n,n);factor=decomposition(M+K,'lu');
        u=initial(field)*ones(n,1);result=zeros(numel(times),21);
        for step=1:round(times(end)/dt)
            q2_validation_deadline(p);
            u=factor\(M*u+rhs);ix=find(times==step*dt);
            if ~isempty(ix),result(ix,:)=u(mesh.output_columns)';end
        end
        errors(level,field)=max(abs(result-exact(:,:,field)),[],'all');
    end
end
assert(all(truncation<1e-9),'CUMCM:q2:BaselineTruncation','Bessel series not converged.');
assert(all(errors(2,:)<errors(1,:)),'CUMCM:q2:BaselineConvergence','Frozen reference does not converge.');
report=struct('kind','FROZEN_INITIAL_COEFFICIENTS_ONLY','time_s',times, ...
    'levels',[20,1;40,0.5],'maximum_errors_T_K_C_kgkg',errors, ...
    'bessel_series_doubling_difference',truncation,'series_terms',seriesTerms, ...
    'observed_reduction',errors(1,:)./errors(2,:), ...
    'constant_boundary_K_kgkg',outside,'passed',true, ...
    'scope','Independent constant-coefficient cylindrical reference, not production physics.');
end
function roots=besselRoots(bi,count,p)
fun=@(x)x.*besselj(1,x)-bi*besselj(0,x);
roots=zeros(count,1);a=0;b=0.1;found=0;
while found<count
    q2_validation_deadline(p);
    if fun(a)*fun(b)<0
        found=found+1;roots(found)=fzero(fun,[a,b]);
    end
    a=b;b=b+0.1;
end
end

function report=validate_q2_pdepe(p,ambient,referenceFields)
started=tic;grids=[40,80];data=cell(2,1);target=linspace(0,p.radius_m,21);
for level=1:2
    r=linspace(0,p.radius_m,grids(level)+1);
    options=odeset('RelTol',10^(-5-2*level),'AbsTol',10^(-7-2*level), ...
        'MaxStep',10/2^(level-1));
    sol=pdepe(1,@equations,@initial,@boundary,r,(0:p.tFinal)',options);
    assert(size(sol,1)==p.tFinal+1&&all(isfinite(sol),'all')&&all(sol>0,'all'), ...
        'CUMCM:q2:PdepeIncomplete','Incomplete pdepe solution.');
    cols=1+(0:20)*(grids(level)/20);
    data{level}=struct('temperature_C',sol(2:end,cols,1)-273.15, ...
        'moisture_kgkg',sol(2:end,cols,2));
end
report=struct('radial_segments',grids,'maximum_step_s',[10,5], ...
    'relative_tolerance',[1e-7,1e-9],'absolute_tolerance',[1e-9,1e-11], ...
    'self_refinement_difference',difference(data{1},data{2}), ...
    'coarse_vs_main',difference(referenceFields,data{1}), ...
    'fine_vs_main',difference(referenceFields,data{2}), ...
    'runtime_seconds',toc(started),'scope','Independent pdepe spatial implementation; no experimental truth claim.', ...
    'output_radius_m',target,'completed',true);
report.fine_fields=data{2};
    function [c,f,s]=equations(~,~,u,dudr)
        q2_validation_deadline(p);
        m=q2_compute("material_properties",u(2),u(1),p);
        c=[m.b;1];f=[m.k*dudr(1);m.D*dudr(2)];s=[0;0];
    end
    function u=initial(~)
        u=[p.T0_K;p.C0];
    end
    function [pl,ql,pr,qr]=boundary(~,~,~,ur,t)
        [Ta,Ce]=q2_compute("interpolate_ambient",t,ambient);
        pl=[0;0];ql=[1;1];pr=[p.hc*(ur(1)-Ta-273.15);p.hm*(ur(2)-Ce)];qr=[1;1];
    end
end
function d=difference(a,b)
d=[max(abs(a.temperature_C-b.temperature_C),[],'all'), ...
    max(abs(a.moisture_kgkg-b.moisture_kgkg),[],'all')];
end

function report=validate_q2_finite_cylinder(p,ambient,referenceFields)
started=tic;levels=[20,20;40,40];values=cell(2,1);
for level=1:2
    Nr=levels(level,1);Nz=levels(level,2);mesh=q2_compute("build_mesh",p.radius_m,Nr);
    z=linspace(0,p.length_m/2,Nz+1);dz=z(2)-z(1);wz=dz*ones(1,Nz+1);wz([1,end])=dz/2;
    nr=Nr+1;nz=Nz+1;n=nr*nz;volume=mesh.w*wz;
    ids=reshape(1:n,nr,nz);left=[reshape(ids(1:end-1,:),[],1);reshape(ids(:,1:end-1),[],1)];
    right=[reshape(ids(2:end,:),[],1);reshape(ids(:,2:end),[],1)];
    geometry=[reshape((mesh.faces/mesh.dr)*wz,[],1);reshape(repmat(mesh.w/dz,1,Nz),[],1)];
    incidence=sparse([left;right],[(1:numel(left))';(1:numel(left))'], ...
        [-ones(numel(left),1);ones(numel(left),1)],n,numel(left));
    area=zeros(nr,nz);area(end,:)=mesh.radius*wz;area(:,end)=area(:,end)+mesh.w;
    adjacency=speye(n)+sparse([left;right],[right;left],1,n,n);
    pattern=kron(spones(adjacency),ones(2));
    y0=repmat([p.T0_K;p.C0],n,1);
    options=odeset('RelTol',1e-7,'AbsTol',1e-9,'MaxStep',10,'JPattern',pattern);
    [times,states]=ode15s(@rhs,(0:p.tFinal)',y0,options);
    assert(times(end)==p.tFinal&&size(states,1)==p.tFinal+1&&all(isfinite(states)&states>0,'all'), ...
        'CUMCM:q2:FiniteCylinderIncomplete','Incomplete 2D sensitivity.');
    radial=mesh.output_columns;
    mid=ids(radial,1);ending=ids(radial,end);
    values{level}=struct('mid_temperature_C',states(2:end,2*mid-1)-273.15, ...
        'mid_moisture_kgkg',states(2:end,2*mid), ...
        'end_temperature_C',states(2:end,2*ending-1)-273.15, ...
        'end_moisture_kgkg',states(2:end,2*ending), ...
        'temperature_range_C',[min(states(:,1:2:end),[],'all'),max(states(:,1:2:end),[],'all')]-273.15, ...
        'moisture_range_kgkg',[min(states(:,2:2:end),[],'all'),max(states(:,2:2:end),[],'all')]);
    clear states
end
v=values{2};report=struct('geometry','AXISYMMETRIC_HALF_CYLINDER', ...
    'end_boundary_assumption','End hc and hm equal lateral values; not measured end-face data.', ...
    'radial_axial_segments',levels,'runtime_seconds',toc(started), ...
    'self_refinement_mid',diffFields(values{1},v,'mid'), ...
    'self_refinement_end',diffFields(values{1},v,'end'), ...
    'fine_mid_vs_radial',[max(abs(v.mid_temperature_C-referenceFields.temperature_C),[],'all'), ...
        max(abs(v.mid_moisture_kgkg-referenceFields.moisture_kgkg),[],'all')], ...
    'fine_end_vs_radial',[max(abs(v.end_temperature_C-referenceFields.temperature_C),[],'all'), ...
        max(abs(v.end_moisture_kgkg-referenceFields.moisture_kgkg),[],'all')], ...
    'temperature_range_C',v.temperature_range_C,'moisture_range_kgkg',v.moisture_range_kgkg, ...
    'no_official_geometry_error_threshold',true,'completed',true);
    function dy=rhs(tNow,y)
        q2_validation_deadline(p);
        T=y(1:2:end);C=y(2:2:end);material=q2_compute("material_properties",C,T,p);
        faces=q2_compute("face_coefficients",C(left),C(right),T(left),T(right),p,0);
        [Ta,Ce]=q2_compute("interpolate_ambient",tNow,ambient);
        heat=incidence*(geometry.*faces.k.*(T(left)-T(right)))+area(:)*p.hc.*(Ta+273.15-T);
        wet=incidence*(geometry.*faces.D.*(C(left)-C(right)))+area(:)*p.hm.*(Ce-C);
        dy=zeros(2*n,1);dy(1:2:end)=heat./(volume(:).*material.b);dy(2:2:end)=wet./volume(:);
    end
end
function e=diffFields(a,b,location)
e=[max(abs(a.([location,'_temperature_C'])-b.([location,'_temperature_C'])),[],'all'), ...
    max(abs(a.([location,'_moisture_kgkg'])-b.([location,'_moisture_kgkg'])),[],'all')];
end

function recovered=q2_recover_validation(cfg,p,ambient,request)
path=string(request.checkpoint);assert(startsWith(path,"artifacts/q2/validation/")&& ...
    ~contains(path,"..")&&~contains(path,":"),'CUMCM:q2:RecoveryInput','Unsafe checkpoint.');
old=jsondecode(fileread(fullfile(cfg.projectRoot,path)));
assert(strcmp(q2_outputs("validation_io",'hash',fullfile(cfg.projectRoot,path)),request.checkpoint_sha256), ...
    'CUMCM:q2:RecoveryInput','Checkpoint changed after recovery selection.');
for item=reshape([old.input_bindings(:);old.code_bindings(:)],1,[])
    actual=string(item.path);
    if actual=="src/q2/validate_q2_resources.m"
        actual=string(request.archived_dispatcher);
    elseif actual=="src/q2/q2_recover_validation.m"&&isfield(request,'archived_recovery_helper')
        actual=string(request.archived_recovery_helper);
    end
    assert(strcmp(q2_outputs("validation_io",'hash',fullfile(cfg.projectRoot,actual)),item.sha256), ...
        'CUMCM:q2:RecoveryInput','Bound source/input differs: %s',item.path);
end
recovered=struct('cache',containers.Map('KeyType','char','ValueType','any'), ...
    'cases',old.cases(:)','trace_bytes',0,'old_report',old,'prior_charge_seconds',request.prior_charge_seconds);
freshBudget=isfield(request,'new_budget_authorized')&&isequal(request.new_budget_authorized,true);
if freshBudget
    assert(isfield(request,'current_run_budget_seconds')&&request.current_run_budget_seconds==3600&& ...
        request.prior_charge_seconds==0&&isfield(request,'authorized_at_utc')&& ...
        strlength(string(request.authorized_at_utc))>0&&isfield(request,'budget_authorization_text')&& ...
        string(request.budget_authorization_text)=="允许新的 60 分钟预算", ...
        'CUMCM:q2:RecoveryInput','A new budget requires the explicit recorded user authorization.');
else
    assert(request.prior_charge_seconds>=old.elapsed_seconds&&request.prior_charge_seconds<3480, ...
        'CUMCM:q2:RecoveryInput','Invalid remaining budget.');
end
bindings=q2_outputs("validation_io",'record',cfg.projectRoot,path);
bindings(end+1)=q2_outputs("validation_io",'record',cfg.projectRoot,request.archived_dispatcher);
if isfield(request,'archived_recovery_helper')
    bindings(end+1)=q2_outputs("validation_io",'record',cfg.projectRoot,request.archived_recovery_helper);
end
for k=1:numel(old.cases)
    meta=old.cases(k);assert(meta.rejected_steps==0&&meta.trace_audit.complete_trace, ...
        'CUMCM:q2:RecoveryInput','Only fully recorded zero-rejection cases can be recovered.');
    mesh=q2_compute("build_mesh",p.radius_m,meta.radial_segments);d=2*(mesh.N+1);
    tracePath=fullfile(cfg.projectRoot,meta.trace_path);info=dir(tracePath);
    assert(info.bytes==meta.trace_bytes&&h5readatt(tracePath,'/accepted','complete')==1);
    shape=h5info(tracePath,'/accepted');count=shape.Dataspace.Size(2);
    assert(count==meta.accepted_steps+1&&shape.Dataspace.Size(1)==d+13);
    temperature=zeros(10801,21);moisture=temperature;filled=false(10801,1);
    maxDiagnostics=zeros(1,12);last=[];lastTime=-Inf;
    dtMin=Inf;dtMax=0;
    for first=1:512:count
        n=min(512,count-first+1);block=h5read(tracePath,'/accepted',[1,first],[d+13,n]);
        assert(all(isfinite(block),'all')&&all(block(2:d+1,:)>0,'all'));
        t=block(1,:);assert(all(diff([lastTime,t])>0));lastTime=t(end);
        atInteger=t==floor(t);row=t(atInteger)+1;
        temperature(row,:)=block(2*mesh.output_columns,atInteger)';
        moisture(row,:)=block(2*mesh.output_columns+1,atInteger)';filled(row)=true;
        diagnostics=block(d+2:end,:)';maxDiagnostics=max(maxDiagnostics,max(abs(diagnostics),[],1));
        positive=diagnostics(:,1)>0;
        if any(positive),dtMin=min(dtMin,min(diagnostics(positive,1)));dtMax=max(dtMax,max(diagnostics(:,1)));end
        last=block(:,end);
    end
    assert(lastTime==10800&&all(filled)&&all(temperature(1,:)==p.T0_K)&&all(moisture(1,:)==p.C0));
    solution=struct('time_s',(0:10800)','temperature_K',temperature,'moisture_kgkg',moisture, ...
        'final_state',last(2:d+1),'final_time_s',10800,'radius_m',mesh.r, ...
        'trace',struct('path',tracePath),'streaming',true, ...
        'accepted_steps',meta.accepted_steps,'rejected_steps',0,'rejected_intervals',zeros(0,2), ...
        'newton_max_iterations',maxDiagnostics(2),'backtracks_max',NaN, ...
        'max_correction',maxDiagnostics(3:4),'max_scaled_residual',maxDiagnostics(5:6), ...
        'max_balance',maxDiagnostics(7:8),'max_balance_sum_difference',[NaN,NaN], ...
        'dt_min_s',dtMin,'dt_max_s',dtMax,'runtime_seconds',meta.runtime_seconds);
    fields=q2_compute("extract_output_grid",solution,mesh);
    settings=struct('maxStep_s',meta.max_step_s,'algebraicLevel',meta.algebraic_level, ...
        'quadratureLevel',meta.quadrature_level);
    c=struct('mesh',mesh,'solution',solution,'fields',fields,'case_id',meta.case_id, ...
        'trace_relative',string(meta.trace_path),'settings',settings);
    recovered.cache(char(meta.case_id))=c;
    recovered.trace_bytes=recovered.trace_bytes+info.bytes;
    bindings(end+1)=q2_outputs("validation_io",'record',cfg.projectRoot,string(meta.trace_path));
end
recovered.bindings=bindings;
recovered.audit=struct('source_checkpoint',path,'source_run_id',old.run_id, ...
    'completed_cases_recovered',numel(old.cases),'partial_case_reused',false, ...
    'prior_budget_charge_seconds',request.prior_charge_seconds, ...
    'prior_charge_is_allocation_not_measured_runtime',true, ...
    'original_formal_exit_code',[],'original_completed_author_receipt_exists',false, ...
    'unknown_fields_not_invented',{{'backtracks_max','max_balance_sum_difference'}}, ...
    'scope','Continue validation from complete accepted-state evidence only; final/public reproduction never uses this mode.');
recovered.audit.new_budget_authorized=freshBudget;
recovered.audit.source_elapsed_seconds=old.elapsed_seconds;
recovered.audit.source_prior_budget_allocation_seconds=[];
if isfield(old,'recovery')&&isfield(old.recovery,'prior_budget_charge_seconds')
    recovered.audit.source_prior_budget_allocation_seconds=old.recovery.prior_budget_charge_seconds;
end
if freshBudget
    recovered.audit.authorized_at_utc=request.authorized_at_utc;
    recovered.audit.budget_note='New user-authorized 3600s run budget; historical elapsed/allocation retained separately, not erased or claimed as zero runtime.';
end
assert(ambient.time_s(1)==0);
end

function q2_validation_deadline(p)
if isfield(p,'validationTimer')
    assert(toc(p.validationTimer)<=p.wallSeconds,'CUMCM:q2:ResourceLimit', ...
        'The shared q2 validation wall-clock budget is exhausted.');
end
end
