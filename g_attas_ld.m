function y = g_attas_ld(x, u, theta)

p = x(1, :);
r = x(2, :);

da = u(1, :);
dr = u(2, :);
b = u(3, :);

Lp = theta(1);
Lr = theta(2);
Lda = theta(3);
Ldr = theta(4);
Lb = theta(5);
Np = theta(6);
Nr = theta(7);
Nda = theta(8);
Ndr = theta(9);
Nb = theta(10);
da_bias = theta(11);
dr_bias = theta(12);
p_bias = theta(13);
r_bias = theta(14);
ay_bias = theta(15);
Yp = theta(16);
Yr = theta(17);
Yda = theta(18);
Ydr = theta(19);
Yb = theta(20);

ay = Yp*p + Yr*r + Yda*(da - da_bias) + Ydr*(dr - dr_bias) + Yb*b;

p_meas = p + p_bias;
r_meas = r + r_bias;
ay_meas = ay + ay_bias;


y = [p_meas; r_meas; ay_meas];