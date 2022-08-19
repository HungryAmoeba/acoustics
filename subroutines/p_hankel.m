function p = p_hankel(r, z, z_0, freq, c, rho, h, beta)
% r = 100;
% z = 20;
% z_0 = 8;

% this section was folded to add default values for 
% freq, c, rho, h, beta
for d=1
    if ~exist('freq', 'var')
        freq = 300;
    end

    if ~exist('c', 'var')
        c = [1500, 1700,  1800];
    end

    if ~exist('rho', 'var')
        rho = [1000, 1500, 2000];
    end

    if ~exist('h', 'var')
        h = 200;
    end

    if ~exist('beta', 'var')
        beta = [.0001,.1, .2];
    end
end

lambda = c/(freq);
alpha = beta./lambda * log(10)/20;


k_arr = 2*pi*freq./c + 1i * alpha;
k = k_arr(1);

fun = @(k_x) 1/2 * 1i./sqrt(k ^2 - k_x .^2) .* ...
    reflection_coeff(k_x, rho, k_arr, h) .* ...
    exp(1i .* sqrt(k^2 - k_x .^ 2) .* (z + z_0)) .* ...
    hankel_zero_order_first(k_x .* r) .* k_x;

p = integral(fun,-Inf,Inf);

end
