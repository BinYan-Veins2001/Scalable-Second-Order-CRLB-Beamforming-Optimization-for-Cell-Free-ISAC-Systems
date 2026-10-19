function [f,c,actual,noise,out]=costar_guard_path_evaluate(cache,t)
%COSTAR_GUARD_PATH_EVALUATE O(U+M) exact-path backtracking evaluation.
% actual=f(0)-f(t) uses polynomial increments and log1p, avoiding subtraction
% of nearby objective values. noise is a scale heuristic, not an interval
% bound. Recheck the accepted rounded matrix W+t*P+t^2*V externally.
assert(isscalar(t) && isreal(t) && isfinite(t),'t must be a finite real scalar.');
powers=[t,t*t,t*t*t,(t*t)*(t*t)];
dg=compensated_rows(cache.gcoef(:,2:5).*powers);
g=cache.base_statistics+dg;
c=compensated_rows(cache.ccoef.*[1,powers]);
mag=cache.gmagnitude(:,2:5)*abs(powers(:));
constraint_noise=cache.noise_factor*(cache.cmagnitude*[1,abs(powers)]');
S=cache.base_sum; D=cache.base_det;
dS=dg(1)+dg(2);
dD=compensated_rows([cache.base_statistics(2)*dg(1), ...
    cache.base_statistics(1)*dg(2),dg(1)*dg(2), ...
    -2*cache.base_statistics(3)*dg(3),-dg(3)^2]);
St=S+dS; Dt=D+dD;
identifiable=isfinite(St) && isfinite(Dt) && St>0 && Dt>cache.det_floor;
out=struct('ga',g(1),'gb',g(2),'gab',g(3),'sum_ab',St,'det',Dt, ...
    'statistics_increment',dg,'identifiable',identifiable, ...
    'feasible',all(isfinite(c)) && all(c<=0),'crlb',inf, ...
    'constraint_noise',constraint_noise,'error_bound_is_validated',false);
f=inf; actual=-inf; noise=inf;
if ~identifiable || ~(isfinite(S) && isfinite(D) && S>0 && D>0); return; end
x=dS/S; y=dD/D;
if strcmpi(cache.objective_type,'log_crlb_ratio')
    f=log(St)-log(Dt);
    actual=log1p(y)-log1p(x);
else
    f=St/Dt;
    actual=(S/D)*(y-x)/(1+y);
end
g0=cache.base_statistics;
sensitivity=(mag(1)+mag(2))/S+ ...
    (abs(g0(2))*mag(1)+abs(g0(1))*mag(2)+2*abs(g0(3))*mag(3)+ ...
    abs(dg(2))*mag(1)+abs(dg(1))*mag(2)+2*abs(dg(3))*mag(3))/D;
if ~strcmpi(cache.objective_type,'log_crlb_ratio'); sensitivity=sensitivity*S/D; end
noise=cache.noise_factor*sensitivity;
out.crlb=cache.crlb_scale*St/Dt;
out.relative_sum_change=x; out.relative_det_change=y;
end

function s=compensated_rows(X)
% Neumaier compensation handles cancellation between polynomial terms.
s=zeros(size(X,1),1); correction=s;
for j=1:size(X,2)
    v=X(:,j); next=s+v;
    larger=abs(s)>=abs(v);
    correction(larger)=correction(larger)+(s(larger)-next(larger))+v(larger);
    smaller=~larger;
    correction(smaller)=correction(smaller)+(v(smaller)-next(smaller))+s(smaller);
    s=next;
end
s=s+correction;
end
