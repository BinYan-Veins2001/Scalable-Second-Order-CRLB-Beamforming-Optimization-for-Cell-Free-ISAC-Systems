function info=paperfig_efficiency_dense_evd(p,W,lambda,Q,sigma)
%PAPERFIG_EFFICIENCY_DENSE_EVD Entire explicit tangent-Hessian EVD workflow.
% Q (working normals), W, lambda and the scalar test threshold are shared
% inputs. Dense Hessian assembly, tangent-basis completion, projection and
% the complete eigensystem are recomputed on every call.
H=paperfig_efficiency_dense_hessian(p,W,lambda); d=size(H,1);
m=size(Q,2);
if m==0; Z=eye(d);
elseif m==d; Z=zeros(d,0);
else; [basis,~]=qr(Q); Z=basis(:,m+1:end); end
Ht=Z'*(H*Z); Ht=(Ht+Ht')/2;
if isempty(Ht)
    info=struct('negative_inertia',0,'lambda_min',inf,'certified',true, ...
        'negative_direction',[],'rayleigh',NaN,'eigenpair_residual',0);
    return;
end
[V,e]=eig(Ht,'vector'); e=real(e); [minimum,index]=min(e);
count=nnz(e+sigma<0); direction=[]; rayleigh=NaN; residual=NaN;
if count>0
    direction=real(Z*V(:,index)); direction=direction/norm(direction);
    rayleigh=direction'*H*direction;
    residual=norm(Ht*V(:,index)-minimum*V(:,index))/max(1,norm(Ht,1));
end
info=struct('negative_inertia',count,'lambda_min',minimum, ...
    'certified',count==0,'negative_direction',direction,'rayleigh',rayleigh, ...
    'eigenpair_residual',residual);
end
