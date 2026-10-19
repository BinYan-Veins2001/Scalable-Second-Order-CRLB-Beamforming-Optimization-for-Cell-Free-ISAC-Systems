function saved_path = save_mat_checkpoint_atomic(target_file, payload)
%SAVE_MAT_CHECKPOINT_ATOMIC Save a checkpoint robustly in a synced folder.
%
% The payload is first written to a new file in the target directory. The
% completed temporary file then replaces the target. Replacement is retried
% because cloud-sync clients can briefly lock an existing checkpoint on
% Windows. If the lock persists, the complete temporary checkpoint is kept
% and its path is reported instead of terminating a long experiment.

if ~isstruct(payload) || ~isscalar(payload)
    error('save_mat_checkpoint_atomic:InvalidPayload', ...
        'payload must be a scalar struct.');
end

target_file = char(target_file);
target_dir = fileparts(target_file);
if isempty(target_dir)
    target_dir = pwd;
end
if ~exist(target_dir,'dir')
    mkdir(target_dir);
end

max_attempts = 12;
temp_file = [tempname(target_dir),'.mat'];
temp_is_ready = false;
last_message = '';

for attempt = 1:max_attempts
    if ~temp_is_ready
        try
            save(temp_file,'-struct','payload');
            temp_is_ready = true;
        catch ME
            last_message = ME.message;
        end
    end

    if temp_is_ready
        try
            [moved,move_message] = movefile(temp_file,target_file,'f');
        catch ME
            moved = false;
            move_message = ME.message;
        end
        if moved
            saved_path = target_file;
            return;
        end
        last_message = move_message;
    end

    if attempt < max_attempts
        pause(min(2,0.25*attempt));
    end
end

if temp_is_ready && isfile(temp_file)
    saved_path = temp_file;
    warning('save_mat_checkpoint_atomic:ReplacementDeferred', ...
        ['Checkpoint replacement remained blocked after %d attempts. ' ...
         'The complete latest checkpoint was preserved at:\n%s\nReason: %s'], ...
        max_attempts,temp_file,last_message);
else
    saved_path = '';
    warning('save_mat_checkpoint_atomic:SaveFailed', ...
        ['Checkpoint could not be written after %d attempts. The experiment ' ...
         'will continue, but the current in-memory progress is not saved. ' ...
         'Reason: %s'],max_attempts,last_message);
end
end
