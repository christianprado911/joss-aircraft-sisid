function K = K_attas_sp(theta)

K11 = theta(10);
K12 = theta(11);
K21 = theta(12);
K22 = theta(13);

K = [K11, K12; K21, K22];
