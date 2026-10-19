function state = costar_paper_geometry(g,c,G,U,nu,o)
%COSTAR_PAPER_GEOMETRY Coordinated working set, NNLS, and inward direction.
% Implements the threshold loop in COSTAR Algorithm 2. All reported
% multipliers, singular values, and inward inequalities use the original
% (unscaled) constraint gradients.
if nargin<6; o=struct(); end
o=local_defaults(o); g=g(:); c=c(:); d=numel(g); q=numel(c);
assert(isequal(size(G),[d,q]) && U>=0 && U<=q && nu>0);
state=struct('active',false(q,1),'lambda',zeros(q,1), ...
    'zeta',g,'Q',zeros(d,0),'gt',g,'rg',norm(g)/max(1,norm(g)), ...
    'r_lambda',0,'kkt',inf,'tau',o.tau_bar, ...
    'sinr_threshold',NaN,'power_threshold',NaN,'sigma_min',inf, ...
    'full_column_rank',true,'positive_multipliers',true, ...
    'normal_inward',zeros(d,1),'success',false, ...
    'message','threshold coordination budget exhausted', ...
    'nnls_residual',inf,'inward_info',struct(), ...
    'work',struct('nnls_calls',0,'threshold_refinements',0, ...
        'inward_qp_calls',0,'inward_qp_sweeps',0));
tau=o.tau_bar;
for refinement=0:o.threshold_max_iter
    tg=min(o.sinr_active_tol,tau); tp=min(o.power_active_tol,tau);
    active=[c(1:U)>-tg;c(U+1:end)>-tp]; A=G(:,active);
    lambda=zeros(q,1); sigma_min=inf; fullrank=true;
    Q=zeros(d,0); nnls_ok=true; nnls_residual=0;
    if ~isempty(A)
        scales=vecnorm(A).';
        safe=max(scales,realmin); As=A./safe.';
        [mu,nnls_ok,nnls_residual]=local_nnls(As,-g,o.nnls_tol);
        lambda(active)=mu./safe;
        state.work.nnls_calls=state.work.nnls_calls+1;
        [Qa,Sa,~]=svd(A,'econ'); singular=diag(Sa);
        rank_threshold=o.rank_tol*max(singular);
        keep=singular>rank_threshold;
        Q=Qa(:,keep);
        fullrank=size(Q,2)==size(A,2);
        if fullrank; sigma_min=min(singular); else; sigma_min=0; end
    end
    zeta=g+G*lambda; gt=g-Q*(Q'*g);
    rg=norm(zeta)/max(1,norm(g));
    rl=norm(lambda.*c,inf)/max(1,norm(lambda,inf));
    inward=costar_paper_inward(A,o);
    state.work.inward_qp_calls=state.work.inward_qp_calls+~isempty(A);
    state.work.inward_qp_sweeps=state.work.inward_qp_sweeps+inward.sweeps;
    state.active=active; state.lambda=lambda; state.zeta=zeta;
    state.Q=Q; state.gt=gt; state.rg=rg; state.r_lambda=rl;
    state.kkt=max([0;c;rg;rl]); state.tau=tau;
    state.sinr_threshold=tg; state.power_threshold=tp;
    state.sigma_min=sigma_min; state.full_column_rank=fullrank;
    state.positive_multipliers=all(lambda(active)>0);
    state.normal_inward=inward.direction;
    state.nnls_residual=nnls_residual; state.inward_info=inward;
    if nnls_ok && inward.success && ...
            tau<=o.kappa_tau*max(rg,nu) && ...
            tau*norm(inward.direction)^2<=1+10*eps
        state.success=true; state.message='working-set coordination passed';
        return;
    end
    if refinement==o.threshold_max_iter; break; end
    tau=tau/2; state.work.threshold_refinements= ...
        state.work.threshold_refinements+1;
    if tau<=realmin || ~isfinite(tau)
        state.message='threshold precision exhausted; stationarity uncertified';
        return;
    end
end
if ~nnls_ok
    state.message='NNLS numerical KKT test failed; stationarity uncertified';
elseif ~inward.success
    state.message=['inward QP unresolved; stationarity uncertified: ' inward.message];
end
end

function [x,ok,residual]=local_nnls(A,b,tolerance)
% Lawson--Hanson active-set NNLS without Optimization Toolbox dependence.
n=size(A,2); x=zeros(n,1); passive=false(n,1);
tol=tolerance*max(1,norm(b)); maxiter=max(30,10*n*n);
for iteration=1:maxiter
    w=A'*(b-A*x); candidates=find(~passive & w>tol);
    if isempty(candidates); break; end
    [~,jj]=max(w(candidates)); passive(candidates(jj))=true;
    for inner=1:maxiter
        z=zeros(n,1); AP=A(:,passive);
        [V,S,R]=svd(AP,'econ'); sv=diag(S);
        keep=sv>max(size(AP))*eps*max([1;sv]);
        z(passive)=R(:,keep)*((V(:,keep)'*b)./sv(keep));
        if all(z(passive)>0); x=z; break; end
        bad=passive & z<=0;
        denominator=x(bad)-z(bad);
        valid=denominator>0;
        if ~any(valid)
            passive(bad)=false; x(bad)=0; break;
        end
        xb=x(bad); alpha=min(xb(valid)./denominator(valid));
        x=x+alpha*(z-x);
        remove=passive & x<=10*eps*max(norm(x,inf),realmin);
        passive(remove)=false; x(remove)=0;
    end
end
r=A*x-b; dual=A'*r;
dual_bad=max([0;-dual]);
comp=norm(x.*dual,inf)/max(1,norm(x,inf));
residual=max(dual_bad,comp)/max(1,norm(b));
ok=all(isfinite(x)) && all(x>=0) && residual<=max(10*tolerance,100*eps);
end

function o=local_defaults(o)
names={'tau_bar','sinr_active_tol','power_active_tol','kappa_tau', ...
    'rank_tol','nnls_tol','inward_tol','inward_max_sweeps','threshold_max_iter'};
values={1e-4,1e-4,1e-5,0.5,1e-11,1e-13,1e-9,2000,80};
for j=1:numel(names)
    if ~isfield(o,names{j}); o.(names{j})=values{j}; end
end
assert(o.tau_bar>0 && o.sinr_active_tol>0 && o.power_active_tol>0 && ...
    o.kappa_tau>0 && o.kappa_tau<1 && o.rank_tol>0 && o.nnls_tol>0 && ...
    o.inward_tol>0 && o.inward_max_sweeps>=0 && o.threshold_max_iter>=0);
end
