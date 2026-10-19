function idx = block_indices(block_id, block_size)
%BLOCK_INDICES Return vector indices for a block.
idx = (block_id-1)*block_size + (1:block_size);
end
