function h=costar_join_histories(a,b,offset)
%COSTAR_JOIN_HISTORIES Preserve every recorded iterate of both RTR stages.
fields={'crlb','objective','tangent_grad_norm','delta','active_count', ...
    'active_change_count','accepted','ratio','step_norm','cg_iter', ...
    'tcg_final_res_norm','lambda_min','kkt','stationarity','primal_residual', ...
    'model_error','restoration_defect','restoration_iterations'};
na=numel(a.iter); nb=numel(b.iter);
h=struct();
h.iter=0:(na+nb-1);
h.time_s=[a.time_s(:);offset+b.time_s(:)].';
h.phase=[repmat("fast",1,na),repmat("accuracy",1,nb)];
h.terminal_start_iter=na;
for j=1:numel(fields)
    name=fields{j};
    h.(name)=[values(a,name,na),values(b,name,nb)];
end
h.z_snapshots=[a.z_snapshots,b.z_snapshots];
h.snapshot_iter=h.iter;
h.active_membership=[a.active_membership,b.active_membership];
if na>0 && nb>0
    h.active_change_count(na+1)=nnz(xor(a.active_membership(:,end), ...
        b.active_membership(:,1)));
end
h.step_type=[a.step_type(:);b.step_type(:)].';
end

function v=values(s,name,n)
if isfield(s,name) && numel(s.(name))==n
    v=s.(name)(:).';
else
    v=nan(1,n);
end
end
