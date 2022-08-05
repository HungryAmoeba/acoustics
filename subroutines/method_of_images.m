function P = method_of_images(r, h, z, z_0, k)
% Study the behavior of a point source in a homogeneous
% fluid layer with impenetrable boundaries

% Upper interface is SOFT
% Lower interface is HARD

% see p 80 of Frisk

% use coordinate system of r and z

% REQUIRE: h, z, z_0, k(freq, c)
% INPUT: r
% returns: p(r)

%freq = 200; c = 1500;
%k = 2*pi*freq/c;

%random inputs for now

%h = 10; z = h/10; z_0 = h/5;

%compute sum
sum = 0; sum_limit = 100;
for n = 0:sum_limit
    R_n1 = sqrt(r ^ 2 + (z - z_0 - 2*n*h)^2);
    R_n2 = sqrt(r^2 + (z + z_0 - 2*(n+1)*h)^2);
    R_n3 = sqrt(r^2 + (z + z_0 + 2 * n * h)^2);
    R_n4 = sqrt(r^2 + (z - z_0 + 2 * (n + 1)*h)^2);
    
    sum = sum + (-1)^n * (exp(1i * k * R_n1)/R_n1 + exp(1i * k * R_n2)/R_n2 - exp(1i * k * R_n3)/R_n3 - exp(1i * k * R_n4)/R_n4);
end

P = sum;

end
