addpath subroutines

r_arr = 5:2.5:100;
z = 20;
z_0 = 8;
p_arr  =zeros(size(r_arr));

for ind = 1:length(r_arr)
    p_arr(ind) = p_hankel(r_arr(ind), z, z_0);
end

figure(1); clf;
plot(r_arr, p_arr);
grid on;
xlabel("r")
ylabel("p(r)")

figure(2); clf;
plot(r_arr, abs(p_arr));
grid on;
xlabel("r")
ylabel("abs(p(r))")