%% 
% n = k1/k = c/c1, m = rho_1/rho
clear all;
 
close all;
%% 
% Reflection from a Lower Velocity, Less Dense Half Space (Water-Air interface) 

rho = 1;
rho_1 = 1.3 * 10^-3;
c = 1500;
c1 = 332;
n = c/c1;
m = rho_1/rho;

theta = 0:pi/1000:pi/2;
y = (rayleigh(m,n,theta));
figure(1);
plot(theta, abs(y), 'LineWidth', 2);
hold on;
plot(theta, angle(y));
xlabel('0 < \theta < \pi/2')
ylabel('|R|')
title("Water-Air Interface")
%% 
% Reflection from a lower velocity, more dense half space (the water bottom 
% interface: model 1)

rho_1 = 1.5;
c1 = 1478;
n = c/c1;
m = rho_1/rho;

theta = 0:pi/1000:pi/2;
y = (rayleigh(m,n,theta));
figure(2); clf;
plot(theta, abs(y),'LineWidth',2);
xlabel('0 < \theta < \pi/2')
ylabel('|R|')
title("Water-Bottom Interface Model 1")
%% 
% Reflection from a Higher Velocity, More Dense Half-Space (the water bottom 
% interface: model 2)

rho_1 = 1.8;
c1 = 1800;
n = c/c1;
m = rho_1/rho;

theta = 0:pi/1000:pi/2;
y = (rayleigh(m,n,theta));
figure(3); clf;
plot(theta, abs(y),'LineWidth',2)
xlabel('0 < \theta < \pi/2')
ylabel('|R|')
title("Water-Bottom Interface Model 2")
%%
%test with random coefficients
rho = [1;1.5;1.7;1.9;2];
c = [1500; 1500; 1575; 1650; 1800];
h = [20;30;8];
freq = 200;
theta = pi/8;

R_array = rayleigh_strat(rho,c,h,theta,freq);

freq_range = 200:2000;
R_vec = [];
for freq = 200:2000
    R_array = rayleigh_strat(rho,c,h,theta,freq);
    R_vec = [R_vec, R_array];
    %R_vec(freq-199) = R_array(1);
    %R_vec_second_layer(freq-199) = R_array(2);
    %R_vec_third_layer(freq-199) = R_array(3);
    %R_vec_fourth_layer(freq-199) = R_array(4);
end

figure(4); clf;
plot(freq_range, abs(R_vec));
xlabel("\omega (Hz)")
ylabel('|R|')
title("Reflection coefficient vs. Frequency")
legend('First Layer', 'Second Layer', 'Third Layer', 'Fourth Layer')

figure(5); clf;

for i = 1:4
    time_domain = ifft(abs(R_vec(i,:)));
    plot(abs(time_domain));
    hold on;
end

title('fft of data');
xlabel("Time")
legend('First Layer', 'Second Layer', 'Third Layer', 'Fourth Layer')
hold off;

figure(6); clf;

for i = 1:4
    time_domain = ifft(abs(R_vec(i,:)));
    plot(abs(time_domain(1:50)));
    hold on;
end

title('fft of data (first 50 points)');
xlabel("Time")
legend('First Layer', 'Second Layer', 'Third Layer', 'Fourth Layer')
hold off;

% create plot of angle vs. frequency
freqs = 10:10:1000;
angle_arr = 0:.01:pi/2-.1;
abs_R = zeros(length(freqs), length(angle_arr));
phase_R = zeros(length(freqs), length(angle_arr));

for freq_index = 1:length(freqs)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freqs(freq_index));
        abs_R(freq_index, angle_index) = abs(R_array(1));
        phase_R(freq_index, angle_index) = angle(R_array(1));
        
    end
end

[X,Y] = meshgrid(freqs, angle_arr);
figure(7); clf;
surf(X, Y, abs_R');
title('|R|')
xlabel('Frequency')
ylabel('Angle')

figure(8); clf;
imagesc([0 1000], [0 pi/2],abs_R);
ylabel('Frequency')
xlabel('Angle')
title('|R|')
x = [5 8];
y = [3 6];
colorbar;

figure(9); clf;
surf(X, Y, phase_R');
title('Phase of R')
xlabel('Frequency')
ylabel('Angle')

figure(10); clf;
imagesc(phase_R);
title('Phase of R')
ylabel('Frequency')
xlabel('Angle')
colorbar;

%% 
% try and get the same plots as page 49 of CoA

c = [1500, 1550, 1800];
rho = [1000, 1500, 2000];
h = 2;
freq = 2000;
angle_arr = 0:pi/200:pi/2;

abs_R = zeros(length(angle_arr), 1);
phase_R = zeros(length(angle_arr), 1);
R_arr = zeros(length(angle_arr), 1);

for angle_index = 1:length(angle_arr)
    R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq);
    abs_R(angle_index) = abs(R_array(1));
    phase_R(angle_index) = angle(R_array(1));
    R_arr(angle_index) = R_array(1);
end

figure(11); clf; 
plot(angle_arr, abs_R);

figure(12); clf;
plot(angle_arr, bottom_loss(R_arr));

angle_arr = 0:pi/200:pi/2;
freq_arr = 0:20:2000;

R_vals = zeros(length(freq_arr), length(angle_arr));

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq_arr(freq_index));
        R_vals(freq_index, angle_index) = R_array(1);
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(13); clf;
surf(X, Y, abs(R_vals));
title('|R|')
xlabel('Frequency')
ylabel('Angle')

figure(14); clf;
imagesc(abs(R_vals));
xlabel('Frequency')
ylabel('Angle')
title('|R|')
x = [5 8];
y = [3 6];
colorbar;

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(15); clf;
surf(X, Y, bottom_loss(R_vals));
title('|R|')
xlabel('Frequency')
ylabel('Angle')