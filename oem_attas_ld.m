%% Constantes
d2r = pi/180;

%% Faz derivada simbolica
% syms th1 th2 th3 th4 th5 th6 th7 th8 th9 a q de;
% xs = [a; q];
% us = [de];
% thetas = [th1; th2; th3; th4; th5; th6; th7; th8; th9];
% fs = f_attas_sp(xs, us, thetas);
% gs = g_attas_sp(xs, us, thetas);
% df_dths = jacobian(fs, thetas);
% dg_dths = jacobian(gs, thetas);
% df_dxs = jacobian(fs, xs);
% dg_dxs = jacobian(gs, xs);
% 
% fgen = matlabFunction(fs, df_dxs, df_dths, 'Vars', {xs, us, thetas}, ...
%                       'File', 'f_attas_sp_gen');
% ggen = matlabFunction(gs, dg_dxs, dg_dths, 'Vars', {xs, us, thetas}, ...
%                       'File', 'g_attas_sp_gen');

%% Carrega os dados
data = load('data/fAttasAilRud1.mat');  % LD

N = size(data.fAttasAilRud1, 1);
dt = diff(data.fAttasAilRud1(1:2, 1));
t = dt* (0:N-1);

ay = data.fAttasAilRud1(:, 3);
beta = data.fAttasAilRud1(:, 14) * d2r;
p = data.fAttasAilRud1(:, 7) * d2r;
r = data.fAttasAilRud1(:, 9) * d2r;
da_l = data.fAttasAilRud1(:, 28) * d2r;
da_r = data.fAttasAilRud1(:, 29) * d2r;
dr = data.fAttasAilRud1(:, 30) * d2r;
da = 0.5 * (da_r - da_l) * d2r;

u = [da, dr, beta];
z = [p, r, ay];

%% Cria as funções do modelo
x0 = [0; 0];

f = @f_attas_ld;
g = @g_attas_ld;

mdlsim = @(theta) euler_sim(x0, t, u, f, g, theta);

% Parametros iniciais Latero-direcional
theta0 = zeros(20, 1);
theta0([1, 7]) = -1;

%% Chama o otimizador
oemopt = oemoptions();
oemopt.R = 'diagonal';
oemopt.deriv = 'complex-step';
fun = @(theta) oem_obj(theta, z, mdlsim, oemopt);

optopt = optimoptions('fminunc');
optopt.Algorithm = 'trust-region';
optopt.Display = 'iter';
optopt.SpecifyObjectiveGradient = true;
optopt.HessianFcn = 'objective';
thetaopt = fminunc(fun, theta0, optopt);

%%
y0 = mdlsim(theta0);
yopt = mdlsim(thetaopt);

for i=1:2
    figure(i);
    clf();
    plot(t, z(:,i), '.', t, yopt(:,i), t, y0(:, i), '--')
    legend({'medicoes', 'estimado', 'inicial'})
end