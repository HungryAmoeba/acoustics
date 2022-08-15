
r = 100;
z = 20;
z_0 = 8;

c = [1500, 1700,  1800];
rho = [1000, 1500, 2000];
freq = 300;
beta = .3;

k = 2 * pi * freq/1500 + .02i;

fun = @(k_x) 1/2 * 1i./sqrt(k ^2 - k_x .^2) * ...
    rayleigh_wavenumber_strat(k_x, rho, c, beta, h, freq) * ...
    exp(1i * sqrt(k^2 - k_x ^ 2) * (z + z_0)) * ...
    hankel_zero_order_first(k_x * r) * k_x;


%function R = rayleigh_wavenumber_strat(k_x, rho_arr, c_arr, beta, h, freq)


%Evaluate the integral from x=0 to x=Inf.
kx = -250:.01:250;
sum = 0;
for k_x = kx
    sum = sum + fun(k_x) * .01;
end
sum


%q = integral(fun,-Inf,Inf)