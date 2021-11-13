function xdot = f_attas_sp(x, u, p)

a = x(1, :);
q = x(2, :);

de = u(1, :);

Za = p(1);
Zq = p(2);
Zde = p(3);
Ma = p(4);
Mq = p(5);
Mde = p(6);
a_eq = p(7);
de_eq = p(8);

adot = Za*(a - a_eq) + Zq*q + Zde*(de - de_eq);
qdot = Ma*(a - a_eq) + Mq*q + Mde*(de - de_eq);

xdot = [adot; qdot];