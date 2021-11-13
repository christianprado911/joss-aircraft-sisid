function [J, grad, Hess, y, x] = oem_obj_rev(p, data, model, options)
if nargin < 4
    options = oemoptions();
end
%% Inicializa dimensoes
N = length(data.t);
ny = size(data.z, 2);
nx = length(model.x0);
npar = length(p);

%% Simula o sistema
dt = diff(data.t);
x = zeros(N, nx);
x(1,:) = model.x0.';

for k=1:N-1
    xk = x(k, :).';
    uk = data.u(k, :).';
    
    % Proximo estado, pelo metodo de Euler    
    xnext = xk + model.f(xk, uk, p)*dt(k);
    x(k+1, :) = xnext.';
end

%% Calcula a saida e o erro
y = model.g(x.', data.u.', p).';
e = data.z - y;

%% Calcula o custo
if isempty(model.R)
    % Estimar a matriz R generica
    R = e.' * e / N;
    S = chol(R, 'lower');
    J = N * sum(log(diag(S)));
    Ri = S.' \ inv(S);
elseif strcmp(model.R, 'diagonal')
    % Estima a matriz R diagonal
    Rd = sum(e.^2, 1);
    R = diag(Rd);
    Ri = diag(1./Rd);
    J = N / 2 * sum(log(Rd));
else
    % Matriz R informada pelo usuario
    Ri = inv(model.R);
    
    J = 0;
    for k = 1:N
        ek = e(k, :).';
        J = J + 0.5 * ek.' * Ri * ek;
    end
end

%% Termina se ha' so' um argumento de saida
if nargout == 1
    return
end

%% Calcula o gradiente
y_ = -e * Ri; % Adjunto da saida

% Inicializa o adjunto do estado e dos parametros
[~, dgN_dx, dgN_dp] = model.g(x(N, :).', data.u(N, :).', p);
p_ = dgN_dp.' * y_(N, :).'; % Adjunto dos parametros

x_ = zeros(N, nx); % Inicializa o adjunto do estado
x_(N, :) = (dgN_dx.' * y_(N, :)').';
for k = N-1:-1:1
    [~, dfk_dx, dfk_dp] = model.f(x(k, :).', data.u(k, :).', p);
    [~, dgk_dx, dgk_dp] = model.g(x(k, :).', data.u(k, :).', p);
    
    xk_ = (dfk_dx*dt(k) + eye(nx)).' * x_(k+1, :).' + dgk_dx.' * y_(k, :)';
    x_(k, :) = (xk_).';
    p_ = p_ + dgk_dp.' * y_(k, :).' + dfk_dp.' * x_(k+1, :).'*dt(k);
end

grad = p_;

if nargout == 2
    return
end

%%
Hess = [];
