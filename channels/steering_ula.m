function a = steering_ula(N, theta, cfg)
%STEERING_ULA ULA steering vector [1, exp(j2pi d/lambda sin theta), ...]^T.
% The vector is intentionally not normalized, matching the convention in the paper.

n = (0:N-1).';
a = exp(1j * 2*pi * cfg.radio.antenna_spacing/cfg.radio.lambda * n * sin(theta));
end
