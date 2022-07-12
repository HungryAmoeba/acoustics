function T = multi_layer_transmission(rho, c, h, theta, freq)
% requires that rho and c are of length 3
% h is a single value

assert(length(rho) == 3)
assert(length(c) == 3)

R_12 = rayleigh(rho(2)/rho(1), c(1)/c(2), theta);
T_12 = transmission_coeff(rho(2)/rho(1), c(1)/c(2), theta);

theta_2 = asin(c(2)/c(1) * sin(theta));

R_23 = rayleigh(rho(3)/rho(2), c(2)/c(3), theta_2);
T_23 = rayleigh(rho(3)/rho(2), c(2)/c(3), theta_2);

T = T_12 * T_23 * exp(2 * 1i * freq/c(2) * h * cos(theta_2));

end