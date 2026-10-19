function Wmat = real_vector_to_complex_matrix(z, N, U)
%REAL_VECTOR_TO_COMPLEX_MATRIX Convert real vector back to complex N-by-U matrix.
%
% Inverse of complex_matrix_to_real_vector:
%   z = [real(Wmat(:)); imag(Wmat(:))]

z = z(:);
L = N * U;
if numel(z) ~= 2*L
    error('real_vector_to_complex_matrix:DimensionMismatch', ...
        'Expected real vector length %d, got %d.', 2*L, numel(z));
end
Wmat = reshape(z(1:L), N, U) + 1j * reshape(z(L+1:end), N, U);
end
