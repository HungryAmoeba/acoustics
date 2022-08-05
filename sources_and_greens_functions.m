% Experiment with some greens functions

freq = 200; c = 1500;
k = 2*pi*freq/c;

%random inputs for now

h = 10; z = h/10; z_0 = h/2;

r_arr = 0:1:1000;

p_arr = zeros(size(r_arr));
for ind = 1:length(r_arr)
    r = r_arr(ind);
    p_arr(ind) = method_of_images(r, h, z, z_0, k);
end

figure(3); clf;
grid on;
subplot(211)
plot(r_arr, 20*log10(abs(p_arr)), '.--')
title("p(r) for fixed z = " + z)
xlabel("r")
subplot(212)
plot(r_arr, angle(p_arr), '.--')


% Look at the z dependence for fixed r = 30

r = 5;
z_arr = 0:h/1000:h;
p_arr = zeros(size(z_arr));
for ind = 1:length(z_arr)
    z = z_arr(ind);
    p_arr(ind) = method_of_images(r, h, z, z_0, k);
end

figure(4); clf;
grid on;
subplot(211)
plot(z_arr, 20*log10(abs(p_arr)), '.--')
title("p(z) for fixed r = " + r)
xlabel("z")
subplot(212)
plot(z_arr, angle(p_arr), '.--')

%%
%plot a couple of r for some fixed z

