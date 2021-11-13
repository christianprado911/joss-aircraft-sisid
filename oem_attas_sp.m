%% Constantes
d2r = pi/180;

%% Faz derivada simbolica
syms th1 th2 th3 th4 th5 th6 th7 th8 th9 a q de;
xs = [a; q];
us = [de];
thetas = [th1; th2; th3; th4; th5; th6; th7; th8; th9];
fs = f_attas_sp(xs, us, thetas);
gs = g_attas_sp(xs, us, thetas);
df_dths = jacobian(fs, thetas);
dg_dths = jacobian(gs, thetas);
df_dxs = jacobian(fs, xs);
dg_dxs = jacobian(gs, xs);

fgen = matlabFunction(fs, df_dxs, df_dths, 'Vars', {xs, us, thetas}, ...
                      'File', 'f_attas_sp_gen');
ggen = matlabFunction(gs, dg_dxs, dg_dths, 'Vars', {xs, us, thetas}, ...
                      'File', 'g_attas_sp_gen');

%% Carrega os dados
data = load('data/fAttasElv1.mat');

t = data.fAttasElv1(:, 1);
u = data.fAttasElv1(:, 22) * d2r;
z = data.fAttasElv1(:, [13, 8]) * d2r;

z_pre = mean(z(t<5,:), 1).'; % z pre excitacao
u_pre = mean(u(t<5,:), 1).'; % u pre excitacao
%% Cria as funções do modelo
x0 = [z_pre(1); 0];

f = @f_attas_sp;
g = @g_attas_sp;
mdlsim = @(theta) euler_sim(x0, t, u, fgen, ggen, theta);

R0 = diag(var(z(t<5, :), 0, 1));
theta0 = [-1; 0; 0; 0; -1; 0; z_pre(1); u_pre(1); 0];

%% Chama o otimizador
oemopt = oemoptions();
oemopt.R = R0;
oemopt.deriv = 'user-fcn';
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