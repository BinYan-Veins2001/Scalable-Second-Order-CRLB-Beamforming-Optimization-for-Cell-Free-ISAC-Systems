function [info,parts]=costar_sosp_structured_psd( ...
    cfg,channels,smat,W,lambda,epsilon_H,Q,options)
%COSTAR_SOSP_STRUCTURED_PSD Deterministic structured curvature calculation.
% With Q supplied, check the ambient Hessian first, then the exact tangent
% inertia on null(Q'). A failed tangent test recovers a unit negative
% direction using the same block factorizations. No ambient Hessian or
% randomized Krylov iteration is constructed.
% Six-argument calls retain the ambient-only interface of earlier releases.
% A numerically unresolved factorization returns inconclusive.
ambient_only=nargin<7;
if nargin<8 || isempty(options); options=struct(); end
N=cfg.sys.M*cfg.sys.Nt; U=cfg.sys.U; d=2*N*U;
if ambient_only; Q=zeros(d,0); end
assert(epsilon_H>0 && isfinite(epsilon_H) && numel(lambda)==U+cfg.sys.M);
assert(size(Q,1)==d && isreal(Q),'Q must be a real d-by-r orthonormal basis.');
if norm(Q'*Q-eye(size(Q,2)),'fro')>1e-8*max(1,size(Q,2))
    error('costar_sosp_structured_psd:NonorthogonalNormals', ...
        'Supply an orthonormal basis for the working-constraint normal space.');
end
lambda=real(lambda(:)); m=size(Q,2); s=3+m;
max_attempts=option(options,'max_shift_attempts',20);
factor_tol=option(options,'factor_tolerance_multiplier',64);
info=struct('status','structured_inconclusive','mode','inconclusive', ...
    'conditional_psd',false,'certified',false, ...
    'floating_point_validated',false,'negative_direction',[], ...
    'rayleigh',NaN,'ritz_residual',NaN,'eigenvalue_lower_bound',-inf, ...
    'hvp_calls',1,'iterations',0,'failure_probability',0, ...
    'effective_dim',d-m,'ambient',true,'structured',true, ...
    'negative_inertia',NaN,'ambient_negative_inertia',NaN, ...
    'tangent_negative_inertia',NaN,'block_negative_inertia',NaN, ...
    'minimum_block_separation',inf,'minimum_block_eigenvalue',inf, ...
    'schur_separation',NaN,'tangent_schur_separation',NaN, ...
    'shift',NaN,'shift_attempts',0,'normal_rank',m, ...
    'norm_bound',NaN,'tangent_residual',NaN, ...
    'direction_shifted_rayleigh',NaN,'outer_rank',3);
parts=struct('D',{cell(U,1)},'Y',{cell(U,1)},'S',eye(3), ...
    'shift',NaN,'beta',NaN,'unshifted_eigenvalues',{cell(U,1)}, ...
    'eigenvectors',{cell(U,1)},'shifted_eigenvalues',{cell(U,1)}, ...
    'Q',Q,'coupling_factor',[]);
[~,out]=exact_lagrangian_hess_vec(cfg,channels,smat,W,zeros(size(W)),lambda);
alpha=out.objective.alpha; B=out.objective.outer_hessian;
[VB,EB]=eig((B+B')/2,'vector'); EB=real(EB);
if any(~isfinite(EB)) || any(EB<=0)
    info.status='structured_outer_positive_definiteness_unresolved'; return;
end
F=VB.*sqrt(EB).'; parts.coupling_factor=F;
matrices={smat.Aa,smat.Ab,smat.Aab};
beta=cfg.radio.effective_bandwidth_sq(:);
A=zeros(N); AW=cell(3,1);
for j=1:3; A=A+alpha(j)*matrices{j}; AW{j}=matrices{j}*W; end
power=2*kron(lambda(U+1:end)/cfg.power.Pmax,ones(cfg.sys.Nt,1));
H=channels.h_stack; all_e=zeros(N*U,1); block_norm=0; coupling_norm_sq=0;
for u=1:U
    weights=cfg.comm.Gamma*lambda(1:U); weights(u)=-lambda(u);
    Du=2*beta(u)*A+(2/cfg.comm.noise_power)*(H.*weights.')*H'+diag(power);
    Du=(Du+Du')/2;
    [Vu,eu]=eig(Du,'vector'); eu=real(eu);
    if any(~isfinite(eu))
        info.status='structured_block_nonfinite'; return;
    end
    J=zeros(N,3);
    for j=1:3; J(:,j)=2*beta(u)*AW{j}(:,u); end
    Yu=J*F;
    parts.D{u}=Du; parts.Y{u}=Yu;
    parts.eigenvectors{u}=Vu; parts.unshifted_eigenvalues{u}=eu;
    all_e((u-1)*N+(1:N))=eu;
    block_norm=max(block_norm,max(abs(eu)));
    coupling_norm_sq=coupling_norm_sq+sum(abs(Yu(:)).^2);
end
info.minimum_block_eigenvalue=min(all_e);
info.norm_bound=block_norm+coupling_norm_sq; parts.beta=info.norm_bound;
if ~ambient_only && m==d
    info.status='structured_empty_tangent'; info.mode='empty_tangent';
    info.conditional_psd=true; info.certified=true; info.ambient=false;
    info.eigenvalue_lower_bound=inf; info.tangent_negative_inertia=0;
    info.negative_inertia=0; return;
end
shifts=admissible_shifts(all_e,epsilon_H,max_attempts);
last_status='structured_block_near_singular';
for attempt=1:numel(shifts)
    sigma=shifts(attempt); info.shift_attempts=attempt;
    info.shift=sigma; parts.shift=sigma;
    gram=zeros(3); DinvC=zeros(d,3); C=zeros(d,3); nD=0;
    min_separation=inf; valid=true;
    for u=1:U
        eu=parts.unshifted_eigenvalues{u}+sigma;
        sep=factor_tol*eps*N*max(1,max(abs(eu)));
        min_separation=min(min_separation,min(abs(eu)));
        if any(abs(eu)<=sep); valid=false; break; end
        Vu=parts.eigenvectors{u};
        ir=(u-1)*N+(1:N); ii=N*U+ir;
        Cu=parts.Y{u};
        DCu=Vu*((Vu'*Cu)./eu);
        gram=gram+real(Cu'*DCu);
        C(ir,:)=real(Cu); C(ii,:)=imag(Cu);
        DinvC(ir,:)=real(DCu); DinvC(ii,:)=imag(DCu);
        parts.shifted_eigenvalues{u}=eu;
        nD=nD+2*nnz(eu<0);
    end
    info.minimum_block_separation=min_separation;
    if ~valid; last_status='structured_block_near_singular'; continue; end
    gram=(gram+gram')/2;
    Ka=-eye(3)-gram(1:3,1:3);
    [~,ea]=eig((Ka+Ka')/2,'vector'); ea=real(ea);
    info.schur_separation=min(abs(ea));
    if unresolved(ea,factor_tol)
        last_status='structured_schur_near_singular'; continue;
    end
    nA=nD+nnz(ea<0)-3;
    info.block_negative_inertia=nD; info.ambient_negative_inertia=nA;
    if nA<0 || nA>d
        last_status='structured_inertia_inconsistent'; continue;
    end
    if nA==0
        info.status='structured_ambient_psd'; info.mode='ambient';
        info.ambient=true;
        info.conditional_psd=true; info.certified=true;
        info.eigenvalue_lower_bound=-sigma; info.negative_inertia=0;
        info.tangent_negative_inertia=0;
        parts=finish_parts(parts,sigma); return;
    end
    if ambient_only
        info.status='structured_ambient_indefinite'; info.mode='ambient';
        info.negative_inertia=nA; parts=finish_parts(parts,sigma); return;
    end
    % Normal-space solves are needed only after the cheaper ambient test.
    C=[C,Q]; DinvC=[DinvC,zeros(d,m)];
    extended_gram=zeros(s); extended_gram(1:3,1:3)=gram;
    for u=1:U
        ir=(u-1)*N+(1:N); ii=N*U+ir;
        Qu=Q(ir,:)+1i*Q(ii,:); Vu=parts.eigenvectors{u};
        eu=parts.shifted_eigenvalues{u};
        DQu=Vu*((Vu'*Qu)./eu);
        DinvC(ir,4:end)=real(DQu); DinvC(ii,4:end)=imag(DQu);
        extended_gram(1:3,4:end)=extended_gram(1:3,4:end)+real(parts.Y{u}'*DQu);
        extended_gram(4:end,4:end)=extended_gram(4:end,4:end)+real(Qu'*DQu);
    end
    extended_gram(4:end,1:3)=extended_gram(1:3,4:end)';
    gram=(extended_gram+extended_gram')/2;
    K=blkdiag(-eye(3),zeros(m))-gram; K=(K+K')/2;
    [VK,eK]=eig(K,'vector'); eK=real(eK);
    info.tangent_schur_separation=min(abs(eK));
    if unresolved(eK,factor_tol)
        last_status='structured_tangent_schur_near_singular'; continue;
    end
    nT=nD+nnz(eK<0)-3-m;
    info.tangent_negative_inertia=nT; info.negative_inertia=nT;
    info.ambient=false;
    if nT<0 || nT>d-m
        last_status='structured_tangent_inertia_inconsistent'; continue;
    end
    if nT==0
        info.status='structured_tangent_psd'; info.mode='tangent';
        info.conditional_psd=true; info.certified=true;
        info.eigenvalue_lower_bound=-sigma;
        parts=finish_parts(parts,sigma); return;
    end
    [direction,tangent_residual]=recover_direction(parts,Q,C,DinvC,VK,eK);
    if isempty(direction)
        last_status='structured_negative_recovery_unresolved'; continue;
    end
    V=real_vector_to_complex_matrix(direction,N,U);
    Hv=complex_matrix_to_real_vector(exact_lagrangian_hess_vec( ...
        cfg,channels,smat,W,V,lambda));
    info.hvp_calls=info.hvp_calls+1;
    rayleigh=real(direction'*Hv);
    if ~isfinite(rayleigh) || rayleigh>=-epsilon_H/2 || tangent_residual>1e-8
        last_status='structured_negative_direction_unresolved'; continue;
    end
    info.status='negative_curvature'; info.mode='negative_curvature';
    info.negative_direction=direction; info.rayleigh=rayleigh;
    info.direction_shifted_rayleigh=rayleigh+sigma;
    info.tangent_residual=tangent_residual;
    info.ritz_residual=norm(Hv-Q*(Q'*Hv)-rayleigh*direction);
    parts=finish_parts(parts,sigma); return;
end
info.status=last_status; info.mode='inconclusive';
end

function shifts=admissible_shifts(e,epsilon_H,budget)
% Prefer 3*epsilon_H/4; subsequent values lie inside spectral gaps.
lo=epsilon_H/2; hi=3*epsilon_H/4;
roots=-e(:); roots=sort(unique(roots(roots>lo & roots<hi)));
edges=[lo;roots;hi]; widths=diff(edges);
[~,order]=sort(widths,'descend');
shifts=hi;
for j=1:numel(order)
    p=order(j); shifts(end+1)=edges(p)+widths(p)/2; %#ok<AGROW>
    if numel(shifts)>=budget; break; end
end
% Alternate values if a Schur complement is singular at a gap midpoint.
for j=1:budget
    shifts(end+1)=lo+(hi-lo)*(j/(budget+1)); %#ok<AGROW>
end
shifts=unique(shifts,'stable'); shifts=shifts(1:min(budget,numel(shifts)));
end

function bad=unresolved(e,multiplier)
tol=multiplier*eps*numel(e)*max(1,max(abs(e)));
bad=any(~isfinite(e)) || any(abs(e)<=tol);
end

function [x,tangent_residual]=recover_direction(parts,Q,C,DinvC,VK,eK)
% Select s+1 negative eigenvectors of diag(Dhat,K), as in Section IV-F.
U=numel(parts.D); N=size(parts.D{1},1); d=2*N*U;
m=size(Q,2); s=m+3; count=s+1;
values=[]; users=[]; columns=[]; imaginary=[];
for u=1:U
    e=parts.shifted_eigenvalues{u}; indices=find(e<0);
    for j=indices(:)'
        values=[values;e(j);e(j)]; %#ok<AGROW>
        users=[users;u;u]; columns=[columns;j;j]; %#ok<AGROW>
        imaginary=[imaginary;false;true]; %#ok<AGROW>
    end
end
indices=find(eK<0);
values=[values;eK(indices)]; users=[users;zeros(numel(indices),1)];
columns=[columns;indices]; imaginary=[imaginary;false(numel(indices),1)];
[~,order]=sort(values,'ascend');
x=[]; tangent_residual=inf;
if numel(order)<count; return; end
Vx=zeros(d,count); Vaux=zeros(s,count);
for k=1:count
    t=order(k); u=users(t); j=columns(t);
    if u==0
        Vaux(:,k)=VK(:,j);
    else
        v=parts.eigenvectors{u}(:,j);
        if imaginary(t); v=1i*v; end
        ir=(u-1)*N+(1:N); ii=N*U+ir;
        Vx(ir,k)=real(v); Vx(ii,k)=imag(v);
    end
end
Wx=Vx-DinvC*Vaux; Wy=Vaux(1:3,:);
L=[Q'*Wx;Wy-C(:,1:3)'*Wx];
row_scale=max(vecnorm(L,2,2),realmin); L=L./row_scale;
[~,S,V]=svd(L); singular=diag(S);
rankL=nnz(singular>max(size(L))*eps*max(1,max(singular)));
eta=V(:,rankL+1:end);
X=Wx*eta; X=X-Q*(Q'*X);
[UX,SX,~]=svd(X,'econ'); ex=diag(SX);
keep=ex>max(size(X))*eps*max(ex);
if ~any(keep); return; end
UX=UX(:,keep); HUX=zeros(size(UX));
for u=1:U
    ir=(u-1)*N+(1:N); ii=N*U+ir;
    Xu=UX(ir,:)+1i*UX(ii,:);
    Vu=parts.eigenvectors{u}; eu=parts.shifted_eigenvalues{u};
    Bu=Vu*(eu.*(Vu'*Xu));
    HUX(ir,:)=real(Bu); HUX(ii,:)=imag(Bu);
end
Y=C(:,1:3); HUX=HUX+Y*(Y'*UX);
R=UX'*HUX; [VR,ER]=eig((R+R')/2,'vector');
[minimum,j]=min(real(ER));
if ~isfinite(minimum) || minimum>=0; return; end
x=UX*VR(:,j); x=real(x/norm(x)); tangent_residual=norm(Q'*x);
end

function parts=finish_parts(parts,sigma)
for u=1:numel(parts.D)
    parts.D{u}=parts.D{u}+sigma*eye(size(parts.D{u}));
end
end

function value=option(options,name,fallback)
if isfield(options,name); value=options.(name); else; value=fallback; end
end
