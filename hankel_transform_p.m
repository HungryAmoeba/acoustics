path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')
function p = p_hankel(r, z, z_0)
% r = 100;
% z = 20;
% z_0 = 8;

c = [1500, 1700,  1800];
rho = [1000, 1500, 2000];
freq = 300;
beta = [.0001,.1, .2];
lambda = c/(freq);
alpha = beta./lambda * log(10)/20;
h=15;

k_arr = 2*pi*freq./c + alpha;
k = k_arr(1);

fun = @(k_x) 1/2 * 1i./sqrt(k ^2 - k_x .^2) .* ...
    rayleigh_wavenumber_strat(k_x, rho, k_arr, h) .* ...
    exp(1i .* sqrt(k^2 - k_x .^ 2) .* (z + z_0)) .* ...
    hankel_zero_order_first(k_x .* r) .* k_x;

p = integral(fun,-Inf,Inf);

end
