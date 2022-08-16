function R = reflection_coeff(k_x, rho, k, h)
% assume that k_x is an array
assert(length(k) == length(rho))
assert(length(h) == length(rho) - 2)

if ~exist('h','var')
    R = rayleigh_wavenumber(k_x, rho_arr(2)/rho_arr(1), k_arr(1), k_arr(2));
    return;
end

num_layers = length(rho);

if isrow(h)
    h = [0 h 0]; %padarray(h, 1);
else
    h = [0; h; 0];
    h = h';
end

m_arr = zeros(num_layers - 1,1); 
R_nat_arr = zeros(num_layers - 1, length(k_x));

for i = 1:num_layers-1
   m_arr(i) = rho(i+1)/rho(i);
   R_nat_arr(i, :) = rayleigh_wavenumber(k_x, m_arr(i), k(i), k(i + 1));
   %R_nat_arr(i) = rayleigh_wavenumber(k_x, m_arr(i), k_i, k_(i+1) )
end

R_prime = zeros(num_layers - 1, length(k_x));
R_prime(num_layers - 1, :) = R_nat_arr(num_layers - 1, :);

for i = num_layers-2:-1:1
    phase = exp(2 * 1i * sqrt(k(i + 1)^2 - k_x.^2) .* h(i + 1) );
    numerator = R_nat_arr(i, :) + R_prime(i + 1, :) .* phase;
    denominator = 1 + R_nat_arr(i, :) .* R_prime(i + 1, :) .* phase;
    R_prime(i, :) = numerator./denominator;
end

R = R_prime(1,:);

end