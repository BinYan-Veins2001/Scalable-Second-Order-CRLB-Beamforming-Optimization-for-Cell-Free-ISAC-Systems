function y_db = pow2db(x)
%POW2DB Convert linear power ratio to dB.
y_db = 10*log10(max(real(x), realmin));
end
