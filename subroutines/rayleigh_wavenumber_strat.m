function R = rayleigh_wavenumber_strat(k_x_arr, rho_arr, k_arr, h)

% Uses: function R = rayleigh_wavenumber(k_x, m, k, k_1)
% Require: k_x, rho_arr, c_arr, beta, h, freq

% Returns: R, an array where R(1) is the overall reflection coefficient
% Returns: R, R(x) is the reflection coeffiient associated to the x-th
% layer
h_copy = h;
%vectorize this code
assert(length(k_arr) == length(rho_arr))
assert(length(h) == length(rho_arr) - 2)
R = zeros(size(k_x_arr));

for ind = 1:length(k_x_arr)
    h = h_copy;
    k_x = k_x_arr(ind);

    if ~exist('h','var')
        R = rayleigh_wavenumber(k_x, rho_arr(2)/rho_arr(1), k_arr(1), k_arr(2));
        return;
    end

    num_layers = length(rho_arr);

    %make sure h is a row vector
    if isrow(h)
        h = [0 h 0]; %padarray(h, 1);
    else
        h = [0; h; 0];
    end
    m_arr = zeros(num_layers - 1,1); 

    R_nat_arr = zeros(num_layers - 1, 1);

    %ensure that frequency gets transformed into angular frequency
    %k_arr = 2*pi*freq./c_arr;
    %calculate attenuation
    %lambda = c_arr./(freq);
    %alpha = beta./lambda * log(10)/20; % imaginary part of k
    %k_arr = k_arr + 1i .* alpha;

    for i = 1:num_layers-1
       m_arr(i) = rho_arr(i+1)/rho_arr(i);
       R_nat_arr(i) = rayleigh_wavenumber(k_x, m_arr(i), k_arr(i), k_arr(i + 1));
       %R_nat_arr(i) = rayleigh_wavenumber(k_x, m_arr(i), k_i, k_(i+1) )
    end

    %create R primed array. calculate backwards, but indicies are in correctly
    %in place
    R_prime = zeros(num_layers - 1, 1);

    R_prime(num_layers-1) = R_nat_arr(num_layers-1);


    for i = num_layers-2:-1:1
        phase = exp(2 * 1i * sqrt(k_arr(i + 1)^2 - k_x^2) * h(i + 1) );
        numerator = R_nat_arr(i) + R_prime(i + 1) * phase;
        denominator = 1 + R_nat_arr(i) * R_prime(i + 1) * phase;
        R_prime(i) = numerator/denominator;

    end

    R(ind) = R_prime(1);
end
end