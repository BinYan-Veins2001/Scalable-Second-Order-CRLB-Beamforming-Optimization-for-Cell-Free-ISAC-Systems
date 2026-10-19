function cache=costar_guard_path_cache(cfg,channels,smat,W,P,V)
%COSTAR_GUARD_PATH_CACHE Exact quartic data on W(t)=W+t*P+t^2*V.
% One setup uses sensing/channel products. Backtracking evaluation thereafter
% needs only three sensing polynomials and U+M constraint polynomials.
% This cache does not certify rounded endpoints: re-evaluate accepted W(t)
% with the ordinary objective/constraint routines before committing it.
assert(isequal(size(W),size(P),size(V)),'Path matrices must have equal sizes.');
U=cfg.sys.U; M=cfg.sys.M; Nt=cfg.sys.Nt;
assert(size(W,2)==U && size(W,1)==M*Nt,'Inconsistent path dimensions.');
beta=cfg.radio.effective_bandwidth_sq(:).';
assert(numel(beta)==U && all(isfinite(beta)) && all(beta>=0));
A={smat.Aa,smat.Ab,smat.Aab};
gcoef=zeros(3,5); gmag=zeros(3,5);
for j=1:3
    AW=A{j}*W; AP=A{j}*P; AV=A{j}*V;
    % Both cross-products are kept, so the expansion also covers a matrix
    % whose Hermitian symmetry has a small floating-point error.
    terms=cat(3,conj(W).*AW,conj(W).*AP+conj(P).*AW, ...
        conj(P).*AP+conj(W).*AV+conj(V).*AW, ...
        conj(P).*AV+conj(V).*AP,conj(V).*AV);
    terms=terms.*reshape(beta,1,U,1);
    gcoef(j,:)=pairwise_sum(reshape(real(terms),[],5));
    % Bounds use unsummed products, before cancellation of cross-terms.
    mags=cat(3,abs(conj(W).*AW),abs(conj(W).*AP)+abs(conj(P).*AW), ...
        abs(conj(P).*AP)+abs(conj(W).*AV)+abs(conj(V).*AW), ...
        abs(conj(P).*AV)+abs(conj(V).*AP),abs(conj(V).*AV));
    gmag(j,:)=pairwise_sum(reshape(mags.*reshape(beta,1,U,1),[],5));
end
Z0=channels.h_stack'*W; Z1=channels.h_stack'*P; Z2=channels.h_stack'*V;
ec=abs_square_coefficients(Z0,Z1,Z2);
em=cat(3,abs(Z0).^2,2*abs(Z0).*abs(Z1), ...
    abs(Z1).^2+2*abs(Z0).*abs(Z2), ...
    2*abs(Z1).*abs(Z2),abs(Z2).^2);
noise=max(cfg.comm.noise_power,realmin); gamma=cfg.comm.Gamma;
ccoef=zeros(U+M,5); cmag=zeros(U+M,5);
for u=1:U
    others=[1:u-1,u+1:U];
    % Sum interference separately to avoid total-minus-signal cancellation.
    desired=reshape(ec(u,u,:),1,5);
    desiredmag=reshape(em(u,u,:),1,5);
    interference=pairwise_sum(reshape(ec(u,others,:),[],5));
    interferencemag=pairwise_sum(reshape(em(u,others,:),[],5));
    ccoef(u,:)=(gamma*interference-desired)/noise;
    ccoef(u,1)=ccoef(u,1)+gamma;
    cmag(u,:)=(gamma*interferencemag+desiredmag)/noise;
    cmag(u,1)=cmag(u,1)+abs(gamma);
end
for m=1:M
    rows=(m-1)*Nt+(1:Nt);
    terms=abs_square_coefficients(W(rows,:),P(rows,:),V(rows,:));
    ccoef(U+m,:)=pairwise_sum(reshape(terms,[],5))/cfg.power.Pmax;
    ccoef(U+m,1)=ccoef(U+m,1)-1;
    mags=cat(3,abs(W(rows,:)).^2,2*abs(W(rows,:)).*abs(P(rows,:)), ...
        abs(P(rows,:)).^2+2*abs(W(rows,:)).*abs(V(rows,:)), ...
        2*abs(P(rows,:)).*abs(V(rows,:)),abs(V(rows,:)).^2);
    cmag(U+m,:)=pairwise_sum(reshape(mags,[],5))/cfg.power.Pmax;
    cmag(U+m,1)=cmag(U+m,1)+1;
end
g=gcoef(:,1); S=g(1)+g(2); D=g(1)*g(2)-g(3)^2;
cache=struct('gcoef',gcoef,'gmagnitude',gmag,'ccoef',ccoef, ...
    'cmagnitude',cmag,'base_statistics',g,'base_sum',S,'base_det',D, ...
    'det_floor',cfg.alm.det_floor,'noise_factor',64*eps*(size(W,1)+U), ...
    'crlb_scale',cfg.sensing.N0*cfg.radio.c^2/(8*pi^2*cfg.sys.Nr*cfg.radio.frame_duration), ...
    'objective_type',cfg.alm.objective_type,'path_degree',4, ...
    'error_bound_is_validated',false);
end

function C=abs_square_coefficients(A,B,C2)
C=cat(3,abs(A).^2,2*real(conj(A).*B), ...
    abs(B).^2+2*real(conj(A).*C2), ...
    2*real(conj(B).*C2),abs(C2).^2);
end

function s=pairwise_sum(X)
% Pairwise real summation has logarithmic accumulation depth.
if isempty(X); s=zeros(1,size(X,2)); return; end
while size(X,1)>1
    n=size(X,1); m=floor(n/2);
    Y=X(1:2:2*m,:)+X(2:2:2*m,:);
    if mod(n,2); Y(end+1,:)=X(end,:); end %#ok<AGROW>
    X=Y;
end
s=X;
end
