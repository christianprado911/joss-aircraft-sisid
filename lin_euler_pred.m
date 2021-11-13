function ypred = lin_euler_pred(x0, t, u, y, f, g, K, p)

%% Inicializa variáveis e tamanho dos vetores
N = length(t);
nx = length(x0);
np = length(p);

%% Simula o sistema
K_ = K(p);
xpred = zeros(N, nx); % Estados do sistema
ypred = zeros(size(y)); % Saida predita
xpred(1, :) = x0.';
for k=1:N-1
    % Saída predita no instante k
    xpredk = xpred(k,:).';
    uk = u(k, :).';
    ypredk = g(xpredk, uk, p);

    % Corrige o estado no instante k
    ek = y(k, :).' - ypredk;
    xcorrk = xpredk + K_ * ek;
    dt = t(k+1) - t(k); % Incremento de tempo (passo de simulacao)
    
    % Calcula a funcao f
    fk = f(xcorrk, uk, p);
    
    % Salva os dados
    xpred(k+1, :) = (xcorrk + fk*dt).';
    ypred(k, :) = ypredk.';
end
