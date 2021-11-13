function [J, grad, Hess] = oem_obj_central(theta, z, R, mdlsim)

N = size(z, 1);
ny = size(z, 2);
npar = length(theta);
y = mdlsim(theta);
e = z - y;

Ri = inv(R);

%% Calcula o custo
J = 0;
for k = 1:N
    ek = e(k, :).';
    J = J + 0.5 * ek.' * Ri * ek; 
end

if nargout == 1
    return
end
%% Calcula a saida perturbada

dy_dp = zeros(N, ny, npar);
for i = 1:length(theta)
    dp = max(abs(theta(i)) * 1e-8, 1e-8);
    
    thetapf = theta;
    thetapb = theta;
    
    thetapf(i) = theta(i) + dp;
    thetapb(i) = theta(i) - dp;

    yf = mdlsim(thetapf);
    yb = mdlsim(thetapb);
    dy_dp(:, :, i) = (yf - yb) / 2 / dp;
end

%% Calcula o gradiente
grad = zeros(npar, 1);
for k = 1:N
    ek = e(k, :).';
    
    dyk_dp = reshape(dy_dp(k, :, :), ny, npar);
    grad = grad - dyk_dp.' * Ri * ek; 
end

if nargout == 2
    return
end


%% Calcula a Hessiana
Hess = zeros(npar, npar);
for k = 1:N
    dyk_dp = reshape(dy_dp(k, :, :), ny, npar);
    Hess = Hess + dyk_dp.' * Ri * dyk_dp; 
end

