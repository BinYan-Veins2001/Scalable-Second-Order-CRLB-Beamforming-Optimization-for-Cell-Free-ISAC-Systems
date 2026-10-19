function [actual,noise_estimate,detail] = costar_sosp_stable_decrease(cfg,smat,W,DW)
%COSTAR_SOSP_STABLE_DECREASE Evaluate f(W)-f(W+DW) without close subtraction.
% noise_estimate is a conservative scale heuristic, NOT an interval bound.
beta=cfg.radio.effective_bandwidth_sq(:).';
A={smat.Aa,smat.Ab,smat.Aab}; g=zeros(3,1); dg=zeros(3,1); mag=zeros(3,1);
for j=1:3
    AW=A{j}*W; AD=A{j}*DW;
    g(j)=sum(beta.*real(sum(conj(W).*AW,1)));
    terms=2*real(conj(DW).*AW)+real(conj(DW).*AD);
    dg(j)=sum(beta.*sum(terms,1));
    mag(j)=sum(beta.*sum(2*abs(conj(DW).*AW)+abs(conj(DW).*AD),1));
end
S=g(1)+g(2); det=g(1)*g(2)-g(3)^2;
dS=dg(1)+dg(2);
dD=g(2)*dg(1)+g(1)*dg(2)+dg(1)*dg(2)-2*g(3)*dg(3)-dg(3)^2;
assert(S>0 && det>0 && S+dS>0 && det+dD>0,'Nonpositive FIM domain.');
x=dS/S; y=dD/det;
if strcmpi(cfg.alm.objective_type,'log_crlb_ratio')
    actual=log1p(y)-log1p(x);
else
    actual=(S/det)*(y-x)/(1+y);
end
sensitivity=(mag(1)+mag(2))/S+ ...
    (abs(g(2))*mag(1)+abs(g(1))*mag(2)+2*abs(g(3))*mag(3)+ ...
    abs(dg(2))*mag(1)+abs(dg(1))*mag(2)+2*abs(dg(3))*mag(3))/det;
if ~strcmpi(cfg.alm.objective_type,'log_crlb_ratio')
    sensitivity=sensitivity*S/det;
end
noise_estimate=64*eps*(size(W,1)+size(W,2))*sensitivity;
detail=struct('statistics',g,'statistics_increment',dg, ...
    'relative_sum_change',x,'relative_det_change',y, ...
    'error_bound_is_validated',false);
end
