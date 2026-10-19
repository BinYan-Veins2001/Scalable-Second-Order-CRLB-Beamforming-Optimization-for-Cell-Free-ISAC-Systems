function [v,info] = costar_guard_scaled_inward(G,q,margin,o)
%COSTAR_GUARD_SCALED_INWARD Minimum-norm directional curvature correction.
% Solves min .5*||v||^2 subject to G'*v <= -q-margin.  The signed entries
% of q are retained: negative constraint curvature need not be corrected.
% The nonnegative dual uses normalized columns and a small Gram matrix.
if nargin<4; o=struct(); end
if ~isfield(o,'inward_tol'); o.inward_tol=1e-9; end
if ~isfield(o,'inward_max_sweeps'); o.inward_max_sweeps=2000; end
if isfield(o,'guard_scaled_max_sweeps'); budget=o.guard_scaled_max_sweeps;
else; budget=o.inward_max_sweeps; end
[d,a]=size(G); v=zeros(d,1); q=q(:);
info=struct('success',false,'message','directional QP budget exhausted', ...
    'multipliers',zeros(a,1),'sweeps',0,'primal_residual',inf, ...
    'scaled_primal_residual',inf,'stationarity_residual',inf, ...
    'complementarity_residual',inf,'dual_residual',0);
assert(isreal(G) && numel(q)==a && all(isfinite(G(:))) && ...
    all(isfinite(q)) && isreal(q));
assert(isreal(margin) && all(isfinite(margin(:))) && ...
    all(margin(:)>=0) && (isscalar(margin) || numel(margin)==a));
assert(isscalar(budget) && budget>=0 && budget==floor(budget) && ...
    isscalar(o.inward_tol) && o.inward_tol>0);
b=q+margin(:); % MATLAB scalar expansion also covers a scalar margin.
if a==0
    % Adding a scalar to an empty column remains empty in supported MATLAB.
    finish([],0,'empty working set',true); return;
end
scales=vecnorm(G).';
if any(~isfinite(scales))
    info.message='nonfinite working-gradient norm'; return;
end
keep=scales>realmin;
if any(~keep & b>0)
    info.message='positive correction required for a zero gradient'; return;
end
if all(b<=0)
    finish(zeros(nnz(keep),1),0,'zero direction is the QP minimizer',true); return;
end
Ak=G(:,keep)./scales(keep).';
bk=b(keep)./scales(keep);
bscale=max(abs(bk));
if ~isfinite(bscale) || bscale<=realmin
    info.message='directional right-side scaling is numerically unresolved'; return;
end
rhs=bk/bscale;
C=Ak'*Ak; C=(C+C')/2; diagonal=diag(C);
if any(~isfinite(C(:))) || any(diagonal<=0) || ~isfinite(bscale)
    info.message='directional QP scaling is numerically unresolved'; return;
end
y=zeros(nnz(keep),1); [~,order]=sort(rhs,'descend');
for sweep=1:budget
    for jj=1:numel(order)
        j=order(jj);
        y(j)=max(0,y(j)+(rhs(j)-C(j,:)*y)/diagonal(j));
    end
    if sweep<=3 || mod(sweep,5)==0
        support=find(y>10*eps*max(norm(y,inf),realmin));
        if ~isempty(support)
            CP=C(support,support); [V,E]=eig(CP); ev=real(diag(E));
            positive=ev>10*eps*max(size(CP))*max([1;abs(ev)]);
            yp=V(:,positive)*((V(:,positive)'*rhs(support))./ev(positive));
            if norm(CP*yp-rhs(support),inf)<= ...
                    max(1e-13,0.01*o.inward_tol)*max(1,norm(rhs(support),inf))
                candidate=zeros(size(y)); candidate(support)=yp;
                negative=candidate<0; alpha=1;
                if any(negative)
                    alpha=min(y(negative)./(y(negative)-candidate(negative)));
                end
                candidate=max(0,y+alpha*(candidate-y));
                old=.5*y'*(C*y)-rhs'*y;
                fresh=.5*candidate'*(C*candidate)-rhs'*candidate;
                if fresh<=old+100*eps*max(1,abs(old)); y=candidate; end
            end
        end
    end
    % Check the dual KKT conditions in the Gram space; a d-vector is only
    % formed on convergence or exit, so each ordinary sweep costs O(a^2).
    Cy=C*y; residual=b(keep)-bscale*scales(keep).*Cy;
    primal=max([0;residual./max(1,abs(b(keep)))]);
    norm2=bscale^2*max(0,y'*Cy);
    comp=norm((bscale*y./scales(keep)).*residual,inf)/max(1,norm2);
    if primal<=o.inward_tol && comp<=o.inward_tol
        finish(y,sweep,'minimum-norm directional KKT tests passed',true); return;
    end
    if any(~isfinite(y)) || norm(y,inf)>1e100
        finish(y,sweep,'directional dual is unbounded or unresolved',false); return;
    end
end
finish(y,budget,'directional QP budget exhausted',false);

    function finish(dual,sweeps,message,proposed_success)
        mu=zeros(a,1);
        if a>0 && ~isempty(dual)
            if all(b<=0); factor=1; else; factor=bscale; end
            mu(keep)=factor*dual./scales(keep);
        end
        v=-G*mu; r=G'*v+b;
        info.multipliers=mu; info.sweeps=sweeps; info.message=message;
        info.primal_residual=max([0;r]);
        info.scaled_primal_residual=max([0;r./max(1,abs(b))]);
        info.stationarity_residual=norm(v+G*mu)/max(1,norm(v));
        info.complementarity_residual=norm(mu.*r,inf)/max(1,norm(v)^2);
        info.dual_residual=max([0;-mu]);
        info.success=proposed_success && all(isfinite(v)) && ...
            info.scaled_primal_residual<=o.inward_tol && ...
            info.stationarity_residual<=o.inward_tol && ...
            info.complementarity_residual<=o.inward_tol && ...
            info.dual_residual<=o.inward_tol;
        if proposed_success && ~info.success
            info.message='original-coordinate directional KKT check failed';
        end
    end
end
