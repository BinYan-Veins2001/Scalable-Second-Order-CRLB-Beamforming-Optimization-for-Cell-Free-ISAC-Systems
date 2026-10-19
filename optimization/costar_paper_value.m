function [f,out]=costar_paper_value(cfg,smat,W)
%COSTAR_PAPER_VALUE Unclipped smooth log-CRLB; no derivative allocation.
beta=cfg.radio.effective_bandwidth_sq(:).';
A={smat.Aa,smat.Ab,smat.Aab}; q=zeros(3,1);
for j=1:3; q(j)=sum(beta.*real(sum(conj(W).*(A{j}*W),1))); end
S=q(1)+q(2); D=q(1)*q(2)-q(3)^2;
out=struct('ga',q(1),'gb',q(2),'gab',q(3),'sum_ab',S,'det',D, ...
    'identifiable',isfinite(S)&&isfinite(D)&&S>0&&D>cfg.alm.det_floor,'crlb',inf);
if ~out.identifiable; f=inf; return; end
f=log(S)-log(D);
out.crlb=cfg.sensing.N0*cfg.radio.c^2/(8*pi^2*cfg.sys.Nr*cfg.radio.frame_duration)*S/D;
end
