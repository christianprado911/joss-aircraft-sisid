function [J, grad, Hess] = oem_obj(theta, z, mdlsim, options)
if nargin < 4
    options = oemoptions();
end
%% Inicializa dimensoes
N = size(z, 1);
ny = size(z, 2);
npar = length(theta);

%% Calcula a saida e o erro
y = mdlsim(theta);
e = z - y;

%% Calcula o custo
if isempty(options.R)
    % Estimar a matriz R generica
    R = e.' * e / N;
    S = chol(R, 'lower');
    J = N * sum(log(diag(S)));
    Ri = S.' \ inv(S);
elseif strcmp(options.R, 'diagonal')
    % Estima a matriz R diagonal
    Rd = sum(e.^2, 1);
    R = diag(Rd);
    Ri = diag(1./Rd);
    J = N / 2 * sum(log(Rd));
else
    % Matriz R informada pelo usuario
    Ri = inv(options.R);
    
    J = 0;
    for k = 1:N
        ek = e(k, :).';
        J = J + 0.5 * ek.' * Ri * ek;
    end
end

if nargout == 1
    return
end
%% Calcula a sensibilidade da saida
if strcmp(options.deriv, 'fwd-diff')
    dy_dp = fwd_diff(N, ny, npar, theta, mdlsim, y, options);
elseif strcmp(options.deriv, 'central-diff')
    dy_dp = central_diff(N, ny, npar, theta, mdlsim, y, options);
elseif strcmp(options.deriv, 'complex-step')
    dy_dp = complex_step(N, ny, npar, theta, mdlsim, options);
elseif strcmp(options.deriv, 'user-fcn')
    [~, ~, dy_dp] = mdlsim(theta);
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

end

function dy_dp = fwd_diff(N, ny, npar, theta, mdlsim, y, options)
dy_dp = zeros(N, ny, npar);
for i = 1:length(theta)
    dp = options.diff_step * max(abs(theta(i)), 1);
    thetap = theta;
    thetap(i) = theta(i) + dp;
    
    yp = mdlsim(thetap);
    dy_dp(:, :, i) = (yp - y) / dp;
end
end

function dy_dp = central_diff(N, ny, npar, theta, mdlsim, y, options)
dy_dp = zeros(N, ny, npar);
for i = 1:length(theta)
    dp = options.diff_step * max(abs(theta(i)), 1);
    thetapf = theta;
    thetapb = theta;
    thetapf(i) = theta(i) + dp;
    thetapb(i) = theta(i) - dp;
    
    ypf = mdlsim(thetapf);
    ypb = mdlsim(thetapb);
    dy_dp(:, :, i) = 0.5*(ypf - ypb) / dp;
end
end

function dy_dp = complex_step(N, ny, npar, theta, mdlsim, options)
dy_dp = zeros(N, ny, npar);
for i = 1:length(theta)
    delta = options.diff_step * max(abs(theta(i)), 1);
    thetap = theta;
    thetap(i) = theta(i) + 1j*delta;
    
    yp = mdlsim(thetap);
    dy_dp(:, :, i) = imag(yp) / delta;
end

end