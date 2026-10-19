function y = db2pow(x_db)
%DB2POW Convert dB to linear power ratio.
y = 10.^(x_db/10);
end
