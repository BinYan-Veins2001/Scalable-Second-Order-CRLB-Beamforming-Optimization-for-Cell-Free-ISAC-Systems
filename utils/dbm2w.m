function p_w = dbm2w(p_dbm)
%DBM2W Convert dBm to watt.
p_w = 10.^((p_dbm - 30)/10);
end
