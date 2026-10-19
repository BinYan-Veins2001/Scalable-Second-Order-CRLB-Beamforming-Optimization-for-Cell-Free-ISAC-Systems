function z = complex_matrix_to_real_vector(Wmat)
%COMPLEX_MATRIX_TO_REAL_VECTOR Convert complex N-by-U matrix to real vector.
%
% The ordering is:
%   z = [real(Wmat(:)); imag(Wmat(:))]
%
% This representation is used by real-domain first-order algorithms.

z = [real(Wmat(:)); imag(Wmat(:))];
end
