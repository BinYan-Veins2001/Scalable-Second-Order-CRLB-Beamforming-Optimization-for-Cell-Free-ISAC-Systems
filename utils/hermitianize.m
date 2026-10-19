function A = hermitianize(A)
%HERMITIANIZE Force a matrix to be Hermitian numerically.
A = 0.5*(A + A');
end
