function xdot = f_attas_ld(x, u, theta)

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

pdot = Lp*p + Lr*r + Lda*(da - da_bias) + Ldr*(dr - dr_bias) + Lb*b;
rdot = Np*p + Nr*r + Nda*(da - da_bias) + Ndr*(dr - dr_bias) + Nb*b;

xdot = [pdot; rdot];