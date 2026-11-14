function varargout = q1_validation(action,varargin)
action=string(action);
assert(isscalar(action),"CUMCM:q1:InvalidValidationAction","验证动作必须是标量字符串。");
switch action
    case "resources"
        varargout{1}=validate_q1_resources(varargin{:});
    case "solution"
        varargout{1}=validate_q1_solution(varargin{:});
    otherwise
        error("CUMCM:q1:InvalidValidationAction","未知验证动作：%s",action);
end
end

function output = validate_q1_resources(cfg,p,ambient)
folder=fullfile(cfg.projectRoot,"artifacts","q1","validation");
stamp=datetime("now","TimeZone","UTC"); stamp.Format="yyyyMMdd'T'HHmmssSSS'Z'";
uuid=char(java.util.UUID.randomUUID()); id="q1_validation_"+string(stamp)+"_"+string(uuid(1:8));
traceDir=fullfile(folder,"internal_traces",id); mkdir(traceDir);
reportPath=fullfile(folder,id+"_diagnostics.json");
relativeReport="artifacts/q1/validation/"+id+"_diagnostics.json";
p.totalTimer=tic;
report=struct("scope","q1","case_id",id,"status","BLOCKED","phase","BASELINE", ...
    "limits",struct("wall_seconds",p.wallSeconds,"working_bytes",p.maxWorkingBytes, ...
    "trace_bytes",p.maxTraceBytes,"max_N",[],"integer_index_limit",p.maxSegments, ...
    "spatial_limit_policy","RESOURCE_ONLY","min_step",p.minStep), ...
    "cases",[],"space_tests",[],"time_tests",[],"joint_tests",[], ...
    "resource_predictions",[],"completed",strings(0,1), ...
    "required",["limits_and_baselines","space","time","algebraic","quadrature", ...
    "joint_rounding","pdepe","finite_cylinder","interpolation_sensitivity","formal_outputs"], ...
    "error_id","","error_message","");
cache=containers.Map('KeyType','char','ValueType','char');
try
    previous=runCase(20,1,1,0); report.baselines=q1_validation_baselines(p,ambient,previous);
    done("limits_and_baselines"); report.phase="SPACE"; N=20; dt=1;
    while true
        assert(2*N<=p.maxSegments,"CUMCM:q1:ResourceLimit","Spatial integer index cannot be represented.");
        current=runCase(2*N,dt,1,0); item=compare(previous,current);
        item.N_coarse=N; item.N_fine=2*N; item.dt_s=dt;
        report.space_tests=[report.space_tests;item]; checkpoint();
        fprintf('Space %d -> %d: heat %.12g; moisture %.12g; budget %.12g\n', ...
            N,2*N,item.maximum_change(1),item.maximum_change(2),p.displayDelta/4);
        N=2*N; previous=current;
        if all(item.maximum_change<p.displayDelta/4), break; end
    end
    done("space"); report.phase="TIME";
    while true
        assert(dt/2>=2*p.minStep,"CUMCM:q1:ResourceLimit","No temporal level left for joint refinement.");
        current=runCase(N,dt/2,1,0); item=compare(previous,current);
        item.N=N; item.dt_coarse_s=dt; item.dt_fine_s=dt/2;
        report.time_tests=[report.time_tests;item]; checkpoint();
        fprintf('Time N%d %.12g -> %.12g: heat %.12g; moisture %.12g; rounded %d/%d\n', ...
            N,dt,dt/2,item.maximum_change,item.rounding_mismatches);
        dt=dt/2; previous=current;
        if all(item.maximum_change<p.displayDelta/4), break; end
    end
    done("time");
    while true
        report.phase="ALGEBRAIC";
        loose=runCase(N,dt,0,0); tight=runCase(N,dt,1,0);
        report.algebraic=compare(loose,tight); checkpoint();
        assert(all(report.algebraic.maximum_change<p.displayDelta/4), ...
            "CUMCM:q1:ConvergenceBudgetFailure","Algebraic discrepancy exceeds budget.");
        done("algebraic"); report.phase="QUADRATURE";
        fineQuadrature=runCase(N,dt,1,1);
        report.quadrature=compare(tight,fineQuadrature); checkpoint();
        assert(all(report.quadrature.maximum_change<p.displayDelta/4), ...
            "CUMCM:q1:ConvergenceBudgetFailure","Quadrature discrepancy exceeds budget.");
        done("quadrature"); report.phase="JOINT_ROUNDING";
        assert(2*N<=p.maxSegments && dt/2>=p.minStep,"CUMCM:q1:ResourceLimit","Joint grid exceeds limits.");
        current=runCase(2*N,dt/2,1,1); item=compare(fineQuadrature,current);
        item.N_coarse=N; item.N_fine=2*N; item.dt_coarse_s=dt; item.dt_fine_s=dt/2;
        report.joint_tests=[report.joint_tests;item]; checkpoint();
        if all(item.maximum_change<p.displayDelta)&&all(item.rounding_mismatches==0), break; end
        N=2*N; dt=dt/2;
    end
    N=2*N; dt=dt/2; done("joint_rounding"); report.phase="INDEPENDENT_PDE";
    report.pdepe=validate_q1_pdepe(p,ambient,current);
    errors=report.space_tests(end).maximum_change+report.time_tests(end).maximum_change ...
        +report.algebraic.maximum_change+report.quadrature.maximum_change;
    independent=report.pdepe.self_changes;
    report.pdepe.explanation_budget=errors+2*independent(end,:);
    assert(all(independent(end,:)<=independent(1,:)) && ...
        all(report.pdepe.difference_from_main(end,:)<=report.pdepe.explanation_budget), ...
        "CUMCM:q1:ConvergenceBudgetFailure","pdepe differences not explained by observed refinements.");
    report.pdepe=rmfield(report.pdepe,["temperature_C","moisture_kgkg"]); done("pdepe");
    report.phase="FINITE_CYLINDER";
    report.finite_cylinder=validate_q1_finite_cylinder(p,ambient,current); done("finite_cylinder");
    report.phase="INTERPOLATION";
    report.interpolation=interpolationSensitivity(N,dt,current); done("interpolation_sensitivity");
    report.phase="FORMAL_OUTPUTS"; report.radial_segments=N; report.max_step_s=dt;
    report.full_field_values=75600; report.status="PASS";
    paths=q1_outputs("write_csv",cfg,current,report);
    paths=[paths,q1_outputs("write_result1",cfg,current),q1_outputs("plot_results",cfg,current,report)];
    assert(toc(p.totalTimer)<p.wallSeconds,"CUMCM:q1:ResourceLimit","Outputs exceeded total time budget.");
    done("formal_outputs"); report.phase="AUTHOR_COMPLETE_AWAITING_31"; checkpoint();
    output=struct("artifacts",[paths,relativeReport],"metrics",report);
catch ex
    report.status="FAIL";
    if any(contains(string(ex.identifier),["ResourceLimit","RoundingInstability"])), report.status="BLOCKED"; end
    report.error_id=string(ex.identifier);
    report.error_message=replace(string(ex.message),string(cfg.projectRoot),".");
    checkpoint(); fprintf('Validation diagnostics: %s\n',relativeReport); rethrow(ex);
end

    function fields=runCase(n,step,a,q)
        key=sprintf('pair_N%d_dt%d_a%d_q%d',n,round(1/step),a,q);
        if isKey(cache,key)
            saved=load(cache(key),"fields"); fields=saved.fields; return
        end
        elapsed=toc(p.totalTimer); predicted=0;
        if ~isempty(report.cases)
            last=report.cases(end);
            predicted=last.runtime_seconds*(n/last.N)*(last.dt_s/step);
        end
        reserve=7*predicted+120;
        entries=dir(fullfile(traceDir,"*.h5")); occupied=sum([entries.bytes]);
        plannedTrace=2*(8*(n+3)*(round(p.tFinal/step)+1)+1024^2);
        plannedMemory=2*8*(2*(p.tFinal+1)*(n+1)+512*(n+3)+100*n+8*1800*21);
        prediction=struct("case_id",string(key),"elapsed_seconds",elapsed, ...
            "estimated_case_seconds",predicted,"reserved_seconds",reserve, ...
            "planned_trace_bytes",plannedTrace,"occupied_trace_bytes",occupied, ...
            "planned_working_bytes",plannedMemory);
        report.resource_predictions=[report.resource_predictions;prediction]; checkpoint();
        assert(elapsed+predicted+reserve<p.wallSeconds,"CUMCM:q1:ResourceLimit", ...
            "Next %s estimated %.1fs plus %.1fs reserve; %.1fs of %.1fs already used.", ...
            key,predicted,reserve,elapsed,p.wallSeconds);
        assert(occupied+plannedTrace<=p.maxTraceBytes && plannedMemory<=p.maxWorkingBytes, ...
            "CUMCM:q1:ResourceLimit","Next %s exceeds storage/memory budget.",key);
        fprintf('Case %s started; elapsed %.1fs; forecast %.1fs.\n',key,elapsed,predicted);
        started=tic; trial=p;
        trial.trace_path=fullfile(traceDir,string(key)+"_heat.h5");
        h=q1_compute("solve_heat",trial,ambient,n,step,a);
        hb=q1_compute("verify_heat_trace",trial.trace_path,trial,h.mesh);
        trial.trace_path=fullfile(traceDir,string(key)+"_moisture.h5");
        c=q1_compute("solve_moisture",trial,ambient,n,step,a,q);
        cb=q1_compute("verify_moisture_trace",trial.trace_path,trial,c.mesh);
        assert(hb.steps==round(p.tFinal/step)&&cb.steps==c.accepted_steps, ...
            "CUMCM:q1:BalanceFailure","Trace does not contain all accepted steps.");
        assert(hb.maximum_balance_J<=10*2*pi*p.length_m*p.rho*p.cp*sum(h.mesh.w)*p.algebraicTolerance(a+1) ...
            &&cb.maximum_balance_kgkg_m2<=10*sum(c.mesh.w)*p.algebraicTolerance(a+1), ...
            "CUMCM:q1:BalanceFailure","Independent full-trace balance failed.");
        fields=q1_compute("extract_output_grid",h,c,n);
        summary=struct("case_id",string(key),"N",n,"dt_s",step,"algebraic_level",a, ...
            "quadrature_level",q,"runtime_seconds",toc(started),"heat_balance",hb, ...
            "moisture_balance",cb,"newton_max_iterations",max(c.iterations), ...
            "rejected_steps",c.rejected_steps,"rejected_intervals",c.rejected_intervals, ...
            "dt_min_s",c.dt_min_s,"dt_max_s",c.dt_max_s, ...
            "heat_range",[min(h.state,[],"all"),max(h.state,[],"all")], ...
            "moisture_range",[min(c.state,[],"all"),max(c.state,[],"all")]);
        report.cases=[report.cases;summary];
        pathValue=fullfile(traceDir,string(key)+"_outputs.mat"); save(pathValue,"fields","summary");
        cache(key)=char(pathValue); checkpoint();
        fprintf('Case %s finished in %.1fs; balances %.5g J / %.5g kgkg*m2.\n', ...
            key,summary.runtime_seconds,hb.maximum_balance_J,cb.maximum_balance_kgkg_m2);
    end
    function result=interpolationSensitivity(n,step,reference)
        dense=ambient; dense.time_s=(0:step:p.tFinal)';
        dense.temperature_C=pchip(ambient.time_s,ambient.temperature_C,dense.time_s);
        dense.moisture_kgkg=pchip(ambient.time_s,ambient.moisture_kgkg,dense.time_s);
        trial=p; trial.trace_path=fullfile(traceDir,"pchip_heat.h5");
        h=q1_compute("solve_heat",trial,dense,n,step,1);
        trial.trace_path=fullfile(traceDir,"pchip_moisture.h5");
        c=q1_compute("solve_moisture",trial,dense,n,step,1,1);
        assert(c.rejected_steps==0,"CUMCM:q1:ConvergenceBudgetFailure", ...
            "PCHIP scenario rejected steps; exact scenario evaluation needs refinement.");
        result=compare(reference,q1_compute("extract_output_grid",h,c,n));
        result.interpretation="Reconstruction sensitivity only; main boundary remains linear.";
    end
    function done(name)
        report.completed=unique([report.completed;string(name)],"stable"); checkpoint();
    end
    function checkpoint()
        report.elapsed_seconds=toc(p.totalTimer);
        list=dir(fullfile(traceDir,"*.h5")); report.trace_bytes=sum([list.bytes]);
        report.unperformed=setdiff(report.required,report.completed,"stable");
        fid=fopen(reportPath,"w","n","UTF-8"); assert(fid>=0);
        cleanup=onCleanup(@() fclose(fid)); fprintf(fid,"%s\n",jsonencode(report,"PrettyPrint",true));
    end
end

function item=compare(coarse,fine)
names=["temperature_C","moisture_kgkg"]; maximum=zeros(1,2); count=zeros(1,2); locations=cell(1,2);
for k=1:2
    difference=abs(coarse.(names(k))-fine.(names(k))); maximum(k)=max(difference,[],"all");
    mismatch=round(coarse.(names(k)),4)~=round(fine.(names(k)),4); count(k)=nnz(mismatch);
    [t,r]=find(mismatch,20); locations{k}=[t,(r-1)/10];
end
item=struct("maximum_change",maximum,"rounding_mismatches",count, ...
    "first_20_mismatch_time_s_radius_cm",{locations});
end

function result=q1_validation_baselines(p,ambient,fields)
short=p; short.tFinal=2; short.hc=0; short.hm=0;
zh=q1_compute("solve_heat",short,ambient,20,1,1); zc=q1_compute("solve_moisture",short,ambient,20,1,1,1);
result.zero_exchange=[max(abs(zh.state-p.theta0),[],"all"),max(abs(zc.state-p.C0),[],"all")];
short=p; short.tFinal=2; equal=ambient;
equal.temperature_C(:)=p.theta0; equal.moisture_kgkg(:)=p.C0;
zh=q1_compute("solve_heat",short,equal,20,1,1); zc=q1_compute("solve_moisture",short,equal,20,1,1,1);
result.steady=[max(abs(zh.state-p.theta0),[],"all"),max(abs(zc.state-p.C0),[],"all")];
assert(all([result.zero_exchange,result.steady]<=p.algebraicTolerance(2)), ...
    "CUMCM:q1:BalanceFailure","Exact limiting scenario failed.");
alpha=p.k/(p.rho*p.cp); D=p.D0*exp(-p.a/p.C0);
result.scales=struct("alpha_m2s",alpha,"D_initial_m2s",D, ...
    "Bi_heat",p.hc*p.radius_m/p.k,"Bi_moisture_initial",p.hm*p.radius_m/D);
times=(1:120)'; rad=(0:20)/20;
heatExact=cylinder(alpha,p.hc*p.radius_m/p.k,p.theta0,p.theta0+10);
moistureExact=cylinder(D,p.hm*p.radius_m/D,p.C0,2);
result.constant_cylinder_diffusivity_m2s=D;
errors=zeros(2,2); boundary=ambient;
boundary.temperature_C(:)=p.theta0+10; boundary.moisture_kgkg(:)=2;
short=p; short.tFinal=120;
for level=0:1
    n=20*2^level; h=q1_compute("solve_heat",short,boundary,n,2^-level,1);
    constant=short; constant.D0=D; constant.a=0;
    c=q1_compute("solve_moisture",constant,boundary,n,2^-level,1,1);
    [~,rows]=ismember(times,h.time_s); [~,crows]=ismember(times,c.time_s);
    errors(level+1,:)=[max(abs(h.state(rows,1:n/20:end)-heatExact),[],"all"), ...
        max(abs(c.state(crows,1:n/20:end)-moistureExact),[],"all")];
end
result.constant_cylinder_errors=errors;
assert(all(errors(2,:)<errors(1,:)),"CUMCM:q1:ConvergenceBudgetFailure","Analytic cylinder refinement failed.");
Ta=interp1(ambient.time_s,ambient.temperature_C,(0:1800)');
lump=zeros(1801,1); lump(1)=p.theta0; rate=2*p.hc/(p.rho*p.cp*p.radius_m);
for k=2:1801, lump(k)=(lump(k-1)+rate*Ta(k))/(1+rate); end
result.lumped_temperature_difference_C=max(abs(fields.temperature_C-lump(2:end)),[],"all");
result.lumped_interpretation="Scale comparison only; Bi not small, not an accuracy certificate.";
    function exact=cylinder(diffusivity,Bi,initial,external)
        fun=@(x) x.*besselj(1,x)-Bi*besselj(0,x);
        brackets=linspace(1e-8,641*pi,21000); values=fun(brackets); roots=[];
        for j=1:numel(brackets)-1
            if values(j)*values(j+1)<0, roots(end+1)=fzero(fun,brackets(j:j+1)); end
            if numel(roots)==640, break; end
        end
        assert(numel(roots)==640,"CUMCM:q1:SolverNonConvergence","Bessel root count.");
        coefficient=2*besselj(1,roots)./(roots.*(besselj(0,roots).^2+besselj(1,roots).^2));
        terms=(exp(-times*(diffusivity*roots.^2/p.radius_m^2)).*coefficient);
        exact=external+(initial-external)*terms*besselj(0,roots'*rad);
        tail=(initial-external)*terms(:,321:end)*besselj(0,roots(321:end)'*rad);
        assert(max(abs(tail),[],"all")<p.displayDelta/4, ...
            "CUMCM:q1:ConvergenceBudgetFailure","Analytic series truncation not small enough.");
    end
end

function report=validate_q1_finite_cylinder(p,ambient,referenceFields)
report=struct("assumption","End-face hc/hm equal side-face coefficients; symmetry at z=0. Not official end data.", ...
    "cases",[],"interpretation","Geometry scenario, no official physical PASS tolerance.");
previous=[];
for Nz=[10,20]
    started=tic; N=20; mesh=q1_compute("build_mesh",p.radius_m,N); nr=N+1; nz=Nz+1;
    dz=(p.length_m/2)/Nz; wz=dz*ones(nz,1); wz([1,end])=dz/2;
    W=mesh.w*wz'; volume=W(:); count=numel(W);
    grid=reshape(1:count,nr,nz);
    left=[reshape(grid(1:end-1,:),[],1);reshape(grid(:,1:end-1),[],1)];
    right=[reshape(grid(2:end,:),[],1);reshape(grid(:,2:end),[],1)];
    area=[reshape(mesh.faces*wz',[],1);repmat(mesh.w,Nz,1)];
    spacing=[mesh.dr*ones(N*nz,1);dz*ones(nr*Nz,1)];
    exterior=zeros(nr,nz); exterior(end,:)=p.radius_m*wz'; exterior(:,end)=exterior(:,end)+mesh.w;
    E=exterior(:); g=p.k*area./spacing;
    K=sparse([left;left;right;right],[left;right;left;right],[g;-g;-g;g],count,count);
    Kh=K+spdiags(p.hc*E,0,count,count);
    Jheat=-spdiags(1./(p.rho*p.cp*volume),0,count,count)*Kh;
    pattern=spones(K+speye(count));
    options=odeset("RelTol",1e-8,"AbsTol",1e-10,"MaxStep",10,"Jacobian",Jheat);
    [~,H]=ode15s(@heat,(0:p.tFinal)',p.theta0*ones(count,1),options);
    options=odeset("RelTol",1e-8,"AbsTol",1e-10,"MaxStep",10,"JPattern",pattern);
    [~,C]=ode15s(@moisture,(0:p.tFinal)',p.C0*ones(count,1),options);
    centre=cat(3,H(2:end,grid(:,1)),C(2:end,grid(:,1)));
    ends=cat(3,H(2:end,grid(:,end)),C(2:end,grid(:,end)));
    assert(all(isfinite(H),"all")&&all(isfinite(C)&C>0,"all"), ...
        "CUMCM:q1:RangeOrFiniteFailure","Invalid finite-cylinder state.");
    change=[NaN,NaN];
    if ~isempty(previous), change=reshape(max(abs(centre-previous),[],[1,2]),1,2); end
    item=struct("radial_segments",N,"axial_half_segments",Nz, ...
        "midplane_difference",difference(centre,referenceFields), ...
        "endplane_difference",difference(ends,referenceFields), ...
        "axial_refinement_change",change,"runtime_seconds",toc(started));
    report.cases=[report.cases;item]; previous=centre;
end
    function du=heat(t,u)
        guard(); [Ta,~]=q1_compute("interpolate_ambient",t,ambient);
        du=Jheat*(u-p.theta0)+p.hc*E*(Ta-p.theta0)./(p.rho*p.cp*volume);
    end
    function du=moisture(t,u)
        guard(); [~,Ce]=q1_compute("interpolate_ambient",t,ambient);
        flux=zeros(size(left));
        radial=spacing==mesh.dr; axial=~radial;
        flux(radial)=q1_compute("moisture_face_flux",u(left(radial)),u(right(radial)),mesh.dr,p,1);
        flux(axial)=q1_compute("moisture_face_flux",u(left(axial)),u(right(axial)),dz,p,1);
        F=area.*flux;
        du=-(accumarray(left,F,[count,1])-accumarray(right,F,[count,1])+p.hm*E.*(u-Ce))./volume;
    end
    function guard()
        if isfield(p,"totalTimer")
            assert(toc(p.totalTimer)<p.wallSeconds,"CUMCM:q1:ResourceLimit","Finite-cylinder time limit.");
        end
    end
end
function d=difference(u,reference)
d=[max(abs(u(:,:,1)-reference.temperature_C),[],"all"), ...
    max(abs(u(:,:,2)-reference.moisture_kgkg),[],"all")];
end

function report=validate_q1_pdepe(p,ambient,referenceFields)
previous=[]; changes=[]; differences=[]; runtimes=[];
for N=[20,40,80]
    started=tic;
    options=odeset("RelTol",1e-9,"AbsTol",1e-11,"MaxStep",10);
    u=pdepe(1,@equations,@initial,@boundary,linspace(0,p.radius_m,N+1), ...
        (0:p.tFinal)',options);
    fields=cat(3,u(2:end,1:N/20:end,1),u(2:end,1:N/20:end,2));
    assert(all(isfinite(fields),"all")&&all(fields(:,:,2)>0,"all"), ...
        "CUMCM:q1:RangeOrFiniteFailure","Invalid pdepe reference.");
    differences=[differences;comparePdepe(fields,referenceFields)];
    if ~isempty(previous), changes=[changes;max(abs(fields-previous),[],[1,2])]; end
    previous=fields; runtimes(end+1)=toc(started);
end
report=struct("N",[20,40,80],"self_changes",reshape(changes,[],2), ...
    "difference_from_main",differences,"runtime_seconds",runtimes, ...
    "temperature_C",previous(:,:,1),"moisture_kgkg",previous(:,:,2), ...
    "interpretation","Independent numerical comparison, not experimental accuracy.");
    function [c,f,s]=equations(~,~,u,dudr)
        guard();
        assert(u(2)>0,"CUMCM:q1:NonPositiveTrial","Nonpositive pdepe state.");
        c=[p.rho*p.cp;1]; f=[p.k*dudr(1);p.D0*exp(-p.a/u(2))*dudr(2)]; s=[0;0];
    end
    function u=initial(~)
        u=[p.theta0;p.C0];
    end
    function [pl,ql,pr,qr]=boundary(~,~,~,ur,t)
        [Ta,Ce]=q1_compute("interpolate_ambient",t,ambient);
        pl=[0;0]; ql=[1;1]; pr=[p.hc*(ur(1)-Ta);p.hm*(ur(2)-Ce)]; qr=[1;1];
    end
    function guard()
        if isfield(p,"totalTimer")
            assert(toc(p.totalTimer)<p.wallSeconds,"CUMCM:q1:ResourceLimit","pdepe time limit.");
        end
    end
end
function d=comparePdepe(fields,reference)
d=[max(abs(fields(:,:,1)-reference.temperature_C),[],"all"), ...
    max(abs(fields(:,:,2)-reference.moisture_kgkg),[],"all")];
end

function report = validate_q1_solution(cfg,p,ambient,candidate,mode)
assert(string(mode)=="mvp","CUMCM:q1:ModelContractConflict","Full validation not implemented in MVP.");
h=candidate.heat; c=candidate.moisture; f=candidate.fields;
assert(all(h.state(1,:)==p.theta0)&&all(c.state(1,:)==p.C0), ...
    "CUMCM:q1:RangeOrFiniteFailure","Initial states changed.");
assert(all(isfinite(h.state),"all")&&all(isfinite(c.state),"all")&&all(c.state>0,"all"), ...
    "CUMCM:q1:RangeOrFiniteFailure","Invalid full-grid state.");
assert(all(h.state>=28 & h.state<=41.513,"all") && ...
    all(c.state>=0.01963 & c.state<=2.55,"all"), ...
    "CUMCM:q1:RangeOrFiniteFailure","Full-grid range failure.");
[Ta,Ce]=q1_compute("interpolate_ambient",ambient.time_s,ambient);
assert(isequal(Ta,ambient.temperature_C)&&isequal(Ce,ambient.moisture_kgkg), ...
    "CUMCM:q1:InputStructureMismatch","Boundary knot interpolation failed.");
ptest=p; ptest.tFinal=2; ptest.hc=0; ptest.hm=0;
zh=q1_compute("solve_heat",ptest,ambient,20,1,1);
zc=q1_compute("solve_moisture",ptest,ambient,20,1,1,0);
zeroHeat=max(abs(zh.state-p.theta0),[],"all");
zeroC=max(abs(zc.state-p.C0),[],"all");
assert(zeroHeat<=p.algebraicTolerance(2)&&zeroC<=p.algebraicTolerance(2), ...
    "CUMCM:q1:BalanceFailure","Zero-exchange exact limit failed.");
eq=ambient; eq.temperature_C(:)=p.theta0; eq.moisture_kgkg(:)=p.C0;
ptest=p; ptest.tFinal=2;
eh=q1_compute("solve_heat",ptest,eq,20,1,1); ec=q1_compute("solve_moisture",ptest,eq,20,1,1,0);
eqHeat=max(abs(eh.state-p.theta0),[],"all"); eqC=max(abs(ec.state-p.C0),[],"all");
assert(eqHeat<=p.algebraicTolerance(2)&&eqC<=p.algebraicTolerance(2), ...
    "CUMCM:q1:BalanceFailure","Uniform steady limit failed.");
L=[0.01963;1;1;2.55]; R=[0.020;1+1e-10;1;1.2];
[flux,dL,dR]=q1_compute("moisture_face_flux",L,R,0.001,p,0);
reverse=q1_compute("moisture_face_flux",R,L,0.001,p,0);
assert(max(abs(flux+reverse))<1e-18 && flux(3)==0 && all(dL>0&dR<0), ...
    "CUMCM:q1:SolverNonConvergence","Flux symmetry/endpoint derivative failed.");
phi=@(x) p.D0*(x.*exp(-p.a./x)-p.a*expint(p.a./x));
reference=-(phi(R)-phi(L))/0.001;
assert(max(abs(flux-reference))<1e-16, ...
    "CUMCM:q1:SolverNonConvergence","Flux identity mismatch.");
nonpositiveRejected=false; extrapolationRejected=false;
try
    q1_compute("moisture_face_flux",0,1,0.001,p,0);
catch ex
    nonpositiveRejected=strcmp(ex.identifier,"CUMCM:q1:NonPositiveTrial");
end
try
    q1_compute("interpolate_ambient",1801,ambient);
catch ex
    extrapolationRejected=strcmp(ex.identifier,"CUMCM:q1:InputStructureMismatch");
end
assert(nonpositiveRejected&&extrapolationRejected, ...
    "CUMCM:q1:InputStructureMismatch","Invalid-input probes were not rejected.");
assert(isequal(f.temperature_C(f.paper_time_s,f.paper_columns), ...
    h.state(f.paper_time_s+1,f.paper_columns)), ...
    "CUMCM:q1:OutputMappingFailure","Paper mapping mismatch.");
report=struct("stage","MVP_ONLY","parameters",p,"input_rows",31, ...
    "full_field_values",numel(f.temperature_C)+numel(f.moisture_kgkg), ...
    "theta_min_C",min(h.state,[],"all"),"theta_max_C",max(h.state,[],"all"), ...
    "C_min_kgkg",min(c.state,[],"all"),"C_max_kgkg",max(c.state,[],"all"), ...
    "heat_balance_max_J",max(abs(h.balance)), ...
    "moisture_balance_max_kgkg_m2",max(abs(c.balance)), ...
    "heat_max_correction_C",h.max_correction,"heat_max_scaled_residual_C",h.max_scaled_residual, ...
    "moisture_max_correction_kgkg",c.max_correction, ...
    "moisture_max_scaled_residual_kgkg",c.max_scaled_residual, ...
    "newton_max_iterations",max(c.iterations),"rejected_steps",c.rejected_steps, ...
    "dt_min_s",c.dt_min_s,"dt_max_s",c.dt_max_s,"radial_segments",20, ...
    "zero_exchange_heat_error_C",zeroHeat,"zero_exchange_moisture_error_kgkg",zeroC, ...
    "steady_heat_error_C",eqHeat,"steady_moisture_error_kgkg",eqC, ...
    "flux_reversal_max",max(abs(flux+reverse)),"flux_phi_max",max(abs(flux-reference)), ...
    "nonpositive_trial_rejected",nonpositiveRejected,"extrapolation_rejected",extrapolationRejected, ...
    "figure_generation","NOT_APPLICABLE_IN_MVP","four_decimal_stability_verified",false, ...
    "unperformed",["refinement","pdepe","finite_cylinder","formal_figures","result1_export"], ...
    "random_seed",cfg.randomSeed);
end
