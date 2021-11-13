%% Constantes
d2r = pi/180;

%% Faz derivada simbolica
syms th1 th2 th3 th4 th5 th6 th7 th8 th9 a q de;
xs = [a; q];
us = [de];
ps = [th1; th2; th3; th4; th5; th6; th7; th8; th9];
fs = f_attas_sp(xs, us, ps);
gs = g_attas_sp(xs, us, ps);
df_dps = jacobian(fs, ps);
dg_dps = jacobian(gs, ps);
df_dxs = jacobian(fs, xs);
dg_dxs = jacobian(gs, xs);

fgen = matlabFunction(fs, df_dxs, df_dps, 'Vars', {xs, us, ps}, ...
                      'File', 'f_attas_sp_gen');
ggen = matlabFunction(gs, dg_dxs, dg_dps, 'Vars', {xs, us, ps}, ...
                      'File', 'g_attas_sp_gen');

%% Carrega os dados
fltdata = load('data/fAttasElv1.mat');

t = fltdata.fAttasElv1(:, 1);
u = fltdata.fAttasElv1(:, 22) * d2r;
z = fltdata.fAttasElv1(:, [13, 8]) * d2r;

z_pre = mean(z(t<5,:), 1).'; % z pre excitacao
u_pre = mean(u(t<5,:), 1).'; % u pre excitacao

data = struct;
data.t = t;
data.u = u;
data.z = z;

%% Cria as funções do modelo
x0 = [z_pre(1); 0];

R0 = diag(var(z(t<5, :), 0, 1));
p0 = [-1; 0; 0; 0; -1; 0; z_pre(1); u_pre(1); 0];

model = struct;
model.f = fgen;
model.g = ggen;
model.x0 = x0;
model.R = R0;

oemopt = oemRevOptions();
obj = @(p) oem_obj_rev(p, data, model, oemopt);

%% Chama o otimizador
optopt = optimoptions('fminunc');
optopt.Algorithm = 'quasi-newton';
optopt.Display = 'iter';
optopt.SpecifyObjectiveGradient = true;
popt = fminunc(obj, p0, optopt);

%%
[~, ~, ~, y0] = obj(p0);
[~, ~, ~, yopt] = obj(popt);

for i=1:2
    figure(i);
    clf();
    plot(t, z(:,i), '.', t, yopt(:,i), t, y0(:, i), '--')
    legend({'medicoes', 'estimado', 'inicial'})
end