function gout = compute_g_values_from_W(cfg, smat, Wcell)
%COMPUTE_G_VALUES_FROM_W Compute ga, gb, gab from lifted matrices W_u.

U = cfg.sys.U;
beta2 = cfg.radio.effective_bandwidth_sq(:);

ga = 0;
gb = 0;
gab = 0;
for u = 1:U
    Wu = Wcell{u};
    ga = ga + beta2(u) * real(trace(smat.Aa * Wu));
    gb = gb + beta2(u) * real(trace(smat.Ab * Wu));
    gab = gab + beta2(u) * real(trace(smat.Aab * Wu));
end

gout.ga = max(real(ga), 0);
gout.gb = max(real(gb), 0);
gout.gab = real(gab);
gout.det = gout.ga*gout.gb - gout.gab^2;
gout.ratio = (gout.ga + gout.gb) / max(gout.det, realmin);
gout.crlb = cfg.sensing.N0 * cfg.radio.c^2 / (8*pi^2*cfg.sys.Nr*cfg.radio.frame_duration) * gout.ratio;
end
