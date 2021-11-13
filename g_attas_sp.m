function y = g_attas_sp(x, u, theta)

a = x(1, :);
q = x(2, :);

q_bias = theta(9);

a_meas = a;
q_meas = q + q_bias;

y = [a_meas; q_meas];