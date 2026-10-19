function [c,G] = costar_sosp_constraints(cfg,channels,W)
%COSTAR_SOSP_CONSTRAINTS Algebraically unchanged, stable interference sums.
U=cfg.sys.U; M=cfg.sys.M; Nt=cfg.sys.Nt; d=2*numel(W);
noise=max(cfg.comm.noise_power,realmin); gamma=cfg.comm.Gamma;
Y=channels.h_stack'*W;
c=zeros(U+M,1);
if nargout>1; G=zeros(d,U+M); end
for u=1:U
    others=[1:u-1 u+1:U];
    interference=sum(abs(Y(u,others)).^2);
    c(u)=(gamma*interference-abs(Y(u,u))^2)/noise+gamma;
    if nargout>1
        signs=gamma*ones(1,U); signs(u)=-1;
        gradient=(2/noise)*channels.h_stack(:,u)*(signs.*Y(u,:));
        G(:,u)=complex_matrix_to_real_vector(gradient);
    end
end
for m=1:M
    rows=(m-1)*Nt+(1:Nt);
    c(U+m)=sum(abs(W(rows,:)).^2,'all')/cfg.power.Pmax-1;
    if nargout>1
        gradient=zeros(size(W)); gradient(rows,:)=2*W(rows,:)/cfg.power.Pmax;
        G(:,U+m)=complex_matrix_to_real_vector(gradient);
    end
end
end
