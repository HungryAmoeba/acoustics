function R = rayleigh_strat(rho, c, h, theta, freq)

%This function calculates the reflection coefficient 
%for an arbitrarily horizontally stratified substrate
%Requires: dim(c) = dim(rho) > 2
%dim(c) = dim(h) + 2
%assumes that coefficients are natural

assert(length(c) == length(rho));
assert(length(c) > 2);
assert(length(c) == length(h) + 2);

num_layers = length(rho);

h = padarray(h, 1);

m_arr = zeros(num_layers - 1,1); n_arr = zeros(num_layers - 1,1);
R_nat_arr = zeros(num_layers - 1, 1);

%ensure that frequency gets transformed into angular frequency
k_arr = 2*pi*freq./c;

theta_arr = zeros(num_layers,1);
theta_arr(1) = theta;
% by snells law, all k sin(theta) = k_1 sin(theta_1)
snell_const = k_arr(1) * sin(theta);

for i = 2:num_layers
    theta_arr(i) = asin(snell_const/ k_arr(i) );
end

for i = 1:num_layers-1
   m_arr(i) = rho(i+1)/rho(i);
   n_arr(i) = c(i)/c(i+1);
   R_nat_arr(i) = rayleigh(m_arr(i),n_arr(i), theta_arr(i));
end


%create R primed array. calculate backwards, but indicies are in correctly
%in place
R_prime = zeros(num_layers - 1, 1);

R_prime(num_layers-1) = R_nat_arr(num_layers-1);


for i = num_layers-2:-1:1
    phase = exp(2 * 1i * k_arr(i + 1) * h(i + 1) * cos(theta_arr(i + 1)));
    numerator = R_nat_arr(i) + R_prime(i + 1) * phase;
    denominator = 1 + R_nat_arr(i) * R_prime(i + 1) * phase;
    R_prime(i) = numerator/denominator;
    
end

R = R_prime;
end