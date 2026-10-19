function result = costar_paper_inward(G,o)
%COSTAR_PAPER_INWARD Minimum-norm direction satisfying G'*h <= -1.
% The nonnegative dual is solved after column and right-hand-side scaling.
% Coordinate minimization is polished on its current positive support.
if nargin<2; o=struct(); end
o=local_defaults(o);
[d,a]=size(G);
result=struct('direction',zeros(d,1),'multipliers',zeros(a,1), ...
    'success',false,'message','inward QP budget exhausted', ...
    'sweeps',0,'primal_residual',inf,'stationarity_residual',inf, ...
    'complementarity_residual',inf,'dual_residual',0);
if a==0
    result.success=true; result.message='empty working set';
    result.primal_residual=0; result.stationarity_residual=0;
    result.complementarity_residual=0; return;
end
scales=vecnorm(G).';
if any(~isfinite(scales)) || any(scales<=realmin)
    result.message='zero or nonfinite working gradient'; return;
end
A=G./scales.'; b=1./scales;
bscale=max(b); rhs=b/bscale;
C=A'*A; C=(C+C')/2;
y=zeros(a,1); [~,order]=sort(rhs,'descend');
for sweep=1:o.inward_max_sweeps
    for jj=1:a
        j=order(jj);
        y(j)=max(0,y(j)+(rhs(j)-C(j,:)*y)/C(j,j));
    end
    % A Newton solve on the positive dual support accelerates correlated
    % constraint normals. A singular or inconsistent support is skipped.
    if sweep<=3 || mod(sweep,5)==0
        support=find(y>10*eps*max(norm(y,inf),realmin));
        if ~isempty(support)
            CP=C(support,support); [V,E]=eig(CP); ev=real(diag(E));
            keep=ev>10*eps*max(size(CP))*max([1;abs(ev)]);
            yp=V(:,keep)*((V(:,keep)'*rhs(support))./ev(keep));
            if norm(CP*yp-rhs(support),inf)<= ...
                    max(1e-13,0.01*o.inward_tol)*max(rhs(support))
                candidate=zeros(a,1); candidate(support)=yp;
                negative=candidate<0;
                alpha=1;
                if any(negative)
                    alpha=min(y(negative)./(y(negative)-candidate(negative)));
                end
                candidate=max(0,y+alpha*(candidate-y));
                old_obj=0.5*y'*(C*y)-rhs'*y;
                new_obj=0.5*candidate'*(C*candidate)-rhs'*candidate;
                if new_obj<=old_obj+100*eps*max(1,abs(old_obj))
                    y=candidate;
                end
            end
        end
    end
    dual_scaled=bscale*y;
    h=-A*dual_scaled;
    mu=dual_scaled./scales;
    constraint=G'*h+1;
    result.direction=h; result.multipliers=mu; result.sweeps=sweep;
    result.primal_residual=max([0;constraint]);
    result.stationarity_residual=norm(h+G*mu)/max(1,norm(h));
    result.complementarity_residual=norm(mu.*constraint,inf)/max(1,norm(h)^2);
    result.dual_residual=max([0;-mu]);
    if result.primal_residual<=o.inward_tol && ...
            result.stationarity_residual<=o.inward_tol && ...
            result.complementarity_residual<=o.inward_tol && ...
            result.dual_residual<=o.inward_tol
        result.success=true; result.message='minimum-norm inward KKT tests passed';
        return;
    end
    if any(~isfinite(h)) || norm(y,inf)>1e100
        result.message='inward QP dual is unbounded or numerically unresolved';
        return;
    end
end
end

function o=local_defaults(o)
if ~isfield(o,'inward_tol'); o.inward_tol=1e-9; end
if ~isfield(o,'inward_max_sweeps'); o.inward_max_sweeps=2000; end
end
