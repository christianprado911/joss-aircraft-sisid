function [y, x, dy_dp, dx_dp] = euler_sim(x0, t, u, f, g, p)

%% Inicializa variáveis e tamanho dos vetores
N = length(t);
nx = length(x0);
np = length(p);

calc_deriv = nargout > 2; % Decide se e' necessario calcular derivadas

%% Simula o sistema

if calc_deriv
    dx_dp = zeros(nx, np, N); % Derivada dos estados em relacao aos param
end

x = zeros(N, nx); % Estados do sistema
x(1, :) = x0.';
for k=1:N-1
    % Estados e entradas no instante k
    xk = x(k, :).';
    uk = u(k, :).';
    dt = t(k+1) - t(k); % Incremento de tempo (passo de simulacao)
    
    % Calcula a funcao f e suas deriv (se necessario derivadas)
    if calc_deriv
        [fk, dfk_dx, dfk_dp] = f(xk, uk, p);
    else
        fk = f(xk, uk, p);
    end
    
    % Proximo estado, pelo metodo de Euler
    x(k+1, :) = (xk + fk*dt).';
    
    if calc_deriv
        % Derivada do metodo de Euler
        dx_dp(:, :, k+1) = dx_dp(:,:, k) + dfk_dp*dt ...
                            + dfk_dx * dx_dp(:,:, k)*dt;
    end
end

%% Calcula a saida e suas derivadas
y = g(x.', u.', p).';

if calc_deriv
    ny = size(y, 2);
    dy_dp = zeros(ny, np, N); % Inicializa a derivada das saidas
    for k=1:N
        xk = x(k, :).';
        uk = u(k, :).';
        
        % Calcula as derivadas da funcao g
        [~, dgk_dx, dgk_dp] = g(xk, uk, p);
        dy_dp(:,:,k) = dgk_dp + dgk_dx * dx_dp(:,:,k);
    end
    dx_dp = shiftdim(dx_dp, 2);
    dy_dp = shiftdim(dy_dp, 2);
end
