function p_dbm = w2dbm(p_w)
%W2DBM Convert watt to dBm.
p_dbm = 10*log10(max(real(p_w), realmin)) + 30;
end
