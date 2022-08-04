% Experiment with some greens functions

freq = 200; c = 1500;
k = 2*pi*freq/c;

%random inputs for now

h = 10; z = h/10; z_0 = h/2;

r_arr = 0:2:1000;

p_arr = zeros(size(r_arr));
for ind = 1:length(r_arr)
    r = r_arr(ind);
    p_arr(ind) = method_of_images(r, h, z, z_0, k);
end

figure(1); clf;
plot(r_arr, real(p_arr))
title("p(r) for fixed z = " + z)
xlabel("r")
% Look at the z dependence for fixed r = 30

r = 5;
z_arr = 0:h/100:h;
p_arr = zeros(size(z_arr));
for ind = 1:length(z_arr)
    z = z_arr(ind);
    p_arr(ind) = method_of_images(r, h, z, z_0, k);
end

figure(2); clf;
plot(z_arr, real(p_arr))
title("p(z) for r = " + r)
xlabel("z")

%%
%plot a couple of r for some fixed z

