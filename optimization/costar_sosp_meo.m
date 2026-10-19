function info = costar_sosp_meo(Hv,dim,project,options)
%COSTAR_SOSP_MEO Random-block matrix-free eigenvalue lower enclosure.
% Let Q contain a uniform random b-dimensional initial subspace, T=Q'*B*Q,
% R=B*Q-Q*T and t_delta^2 be a lower Beta(b/2,(n-b)/2) quantile. Then
% lambda_min(B)>=lambda_min(T)-||R||_F/t_delta with probability >=1-delta.
% The SAME overlap event covers every adaptive subspace containing that
% initial block. This is not a Ritz-vector-residual certificate.
% Assumptions: exact symmetric HVP, orthogonal projection, known rank, and
% ideal independent Gaussian samples. Floating-point errors are not bounded.
o=fill(options,dim);
assert(o.epsilon_H>0 && isfinite(o.norm_bound) && o.norm_bound>=0);
assert(o.failure_probability>0 && o.failure_probability<1);
assert(o.max_iter>=1 && o.effective_dim>=0);
before=rng; cleanup=onCleanup(@() rng(before)); %#ok<NASGU>
rng(o.seed,'twister'); n=o.effective_dim;
kmax=min(n,o.max_iter); b=min([o.probes,n,kmax]);
info=struct('status','inconclusive','conditional_psd',false, ...
    'floating_point_validated',false,'negative_direction',[], ...
    'rayleigh',NaN,'ritz_residual',NaN,'eigenvalue_lower_bound',NaN, ...
    'budget',n,'iterations',0,'hvp_calls',0,'budget_capped',kmax<n, ...
    'failure_probability',o.failure_probability,'norm_bound',o.norm_bound, ...
    'effective_dim',n,'restarts',0,'basis_orthogonality_error',NaN, ...
    'probes',b,'overlap_lower_bound',NaN,'invariance_defect',NaN, ...
    'compression_minimum',NaN,'spectral_error_allowance',NaN);
info.trace=struct('basis_dim',[],'compression_minimum',[], ...
    'invariance_defect',[],'eigenvalue_lower_bound',[],'rayleigh',[]);
if n==0
    info.status='empty_tangent'; info.conditional_psd=true;
    info.eigenvalue_lower_bound=inf; return;
end
if b==n
    overlap=1;
else
    overlap=o.overlap_safety_factor*sqrt( ...
        betaincinv(o.failure_probability,b/2,(n-b)/2));
end
info.overlap_lower_bound=overlap;
Q=zeros(dim,kmax); HQ=zeros(dim,kmax);
G=project(randn(dim,b)); [Q0,S,~]=svd(G,'econ');
if min(diag(S))<1e-10
    info.status='random_block_rank_failure'; return;
end
Q(:,1:b)=Q0(:,1:b); basis_count=b;
checkpoint=min(kmax,max(12,b));
for k=1:kmax
    if k>basis_count; info.status='numerical_basis_exhausted'; return; end
    q=Q(:,k); h=project(Hv(q)); HQ(:,k)=h;
    info.hvp_calls=info.hvp_calls+1;
    r=h;
    for pass=1:2
        r=r-Q(:,1:basis_count)*(Q(:,1:basis_count)'*r);
    end
    if basis_count<kmax && norm(r)>o.breakdown_tol*max(1,norm(h))
        basis_count=basis_count+1;
        Q(:,basis_count)=r/norm(r);
    end
    exhausted=k==basis_count;
    if k==checkpoint || k==kmax || exhausted
        Qk=Q(:,1:k); Hk=HQ(:,1:k);
        T=Qk'*Hk; T=(T+T')/2;
        [V,D]=eig(T,'vector'); [theta,j]=min(real(D));
        v=project(Qk*V(:,j)); v=v/norm(v);
        hv=project(Hv(v)); info.hvp_calls=info.hvp_calls+1;
        rayleigh=real(v'*hv);
        residual=norm(Hk-Qk*T,'fro');
        spectral_error=residual/overlap;
        if k==n; spectral_error=residual; end
        lower=theta-spectral_error;
        info.rayleigh=rayleigh; info.ritz_residual=norm(hv-rayleigh*v);
        info.iterations=k; info.invariance_defect=residual;
        info.compression_minimum=theta; info.spectral_error_allowance=spectral_error;
        info.eigenvalue_lower_bound=lower;
        if o.record_history
            info.trace.basis_dim(end+1)=k;
            info.trace.compression_minimum(end+1)=theta;
            info.trace.invariance_defect(end+1)=residual;
            info.trace.eigenvalue_lower_bound(end+1)=lower;
            info.trace.rayleigh(end+1)=rayleigh;
        end
        if rayleigh<=-o.epsilon_H/2
            info.status='negative_curvature'; info.negative_direction=v;
            return;
        end
        if lower>=-o.epsilon_H
            info.conditional_psd=true;
            info.status='random_block_psd_conditional';
            if k==n; info.status='full_span_psd_conditional'; end
            info.basis_orthogonality_error=norm(Qk'*Qk-eye(k),'fro');
            return;
        end
        checkpoint=min(kmax,max(k+1,2*checkpoint));
    end
    if exhausted && k<kmax
        candidate=project(randn(dim,1));
        for pass=1:2
            candidate=candidate-Q(:,1:basis_count)*(Q(:,1:basis_count)'*candidate);
        end
        if norm(candidate)<=1e-10
            info.status='numerical_basis_exhausted'; return;
        end
        basis_count=basis_count+1;
        Q(:,basis_count)=candidate/norm(candidate);
        info.restarts=info.restarts+1;
    end
end
info.iterations=k;
info.status='budget_exhausted';
end

function o=fill(o,d)
names={'epsilon_H','norm_bound','failure_probability','max_iter', ...
    'effective_dim','seed','breakdown_tol','probes','record_history', ...
    'overlap_safety_factor'};
values={1e-5,NaN,1e-3,inf,d,2026,1e-13,4,false,0.5};
for j=1:numel(names)
    if ~isfield(o,names{j}); o.(names{j})=values{j}; end
end
assert(o.overlap_safety_factor>0 && o.overlap_safety_factor<=1);
end
