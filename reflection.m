%% 
% n = k1/k = c/c1, m = rho_1/rho

% all angles are incident angles (follow Frisks' convention, where 0 is
% normal incident and pi/2 is horizontal)
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
plot(theta * 180/pi, abs(y), 'LineWidth', 2);
hold on;

plot(theta * 180/pi, angle(y), '--');
xlabel('\theta')
title("Water-Air Interface")

legend("|R|", "Phase",'Location','northwest', 'FontSize', 18)
hold off;
saveas(gcf,'Water-Air Interface.png')
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
plot(theta * 180/pi, abs(y),'LineWidth',2); hold on;
yyaxis right
ylim([0 pi])
plot(theta * 180/pi, angle(y),'--')
xlabel('\theta')

legend("|R|", "Phase",'Location','northwest', 'FontSize', 18)
title("Water-Bottom Interface Model 1")
hold off
saveas(gcf,'Water-Bottom Interface Model 1.png')
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
plot(theta * 180/pi, abs(y),'LineWidth',2)
xlabel('\theta')
ylabel('|R|')

yyaxis right

plot(theta * 180/pi, angle(y), '--')
xlabel('\theta')

legend("|R|", "Phase",'Location','northwest', 'FontSize', 18)
title("Water-Bottom Interface Model 2")
saveas(gcf,'Water-Bottom Interface Model 2.png')
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
% saveas(gcf,'Figure 4.png')
% figure(5); clf;
% 
% for i = 1:4
%     time_domain = ifft(abs(R_vec(i,:)));
%     plot(abs(time_domain));
%     hold on;
% end
% 
% title('fft of data');
% xlabel("Time")
% legend('First Layer', 'Second Layer', 'Third Layer', 'Fourth Layer')
% hold off;
% 
% figure(6); clf;

% for i = 1:4
%     time_domain = ifft(abs(R_vec(i,:)));
%     plot(abs(time_domain(1:50)));
%     hold on;
% end
% 
% title('fft of data (first 50 points)');
% xlabel("Time")
% legend('First Layer', 'Second Layer', 'Third Layer', 'Fourth Layer')
% hold off;

% create plot of angle vs. frequency
freqs = 10:10:1000;
angle_arr = 0:.01:pi/2;
abs_R = zeros(length(freqs), length(angle_arr));
phase_R = zeros(length(freqs), length(angle_arr));

for freq_index = 1:length(freqs)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freqs(freq_index));
        abs_R(freq_index, angle_index) = abs(R_array(1));
        phase_R(freq_index, angle_index) = angle(R_array(1));
        
    end
end

% are the axis labeled correctly?
[X,Y] = meshgrid(freqs, angle_arr);
figure(7); clf;
surf(X, Y * 180/pi, abs_R');
title('|R|')
xlabel('Frequency')
ylabel('Angle')
% saveas(gcf,'Figure 7.png')

figure(8); clf;
imagesc([0 90],[0 1000], abs_R);
ylabel('Frequency')
xlabel('Angle')
title('|R|')
colorbar;
% saveas(gcf, 'Figure 8.png')

figure(9); clf;
surf(X, Y * 180/pi, phase_R');
title('Phase of R')
xlabel('Frequency')
ylabel('Angle')
% saveas(gcf, 'Figure 9.png')

figure(10); clf;
imagesc([0 90],[0 1000],phase_R);
title('Phase of R')
ylabel('Frequency')
xlabel('Angle')
colorbar;
% saveas(gcf, 'Figure 10.png')

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
plot(angle_arr * 180/pi, abs_R);
title("|R| for 2 kHz")
xlabel("angle")
% saveas(gcf, 'Figure 11.png')

figure(12); clf;
plot(angle_arr * 180/pi, bottom_loss(R_arr));
title("Bottom loss for 2 kHz")
xlabel("angle")
% saveas(gcf, 'Figure 12.png')

angle_arr = 0:pi/200:pi/2;
freq_arr = 0:20:2000;

R_vals = zeros(length(freq_arr), length(angle_arr));

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq_arr(freq_index));
        R_vals(freq_index, angle_index) = R_array(1);
    end
end

% why does it look like x and y axis are flipped?
[X,Y] = meshgrid(freq_arr, angle_arr);
figure(13); clf;
surf(X, Y * 180/pi, abs(R_vals)');
title('|R|')
xlabel('Frequency')
ylabel('Angle')
% saveas(gcf, 'Figure 13.png')

figure(14); clf;
imagesc([0 90], [0 2000], abs(R_vals));
xlabel('Angle')
ylabel('Frequency')
title('|R|')
colorbar;
% saveas(gcf, 'Figure 14.png')

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(15); clf;
surf(X, Y * 180/pi, bottom_loss(R_vals)');
title('Bottom Loss')
xlabel('Frequency')
ylabel('Angle')
% saveas(gcf, 'Figure 15.png')

figure(16); clf;
imagesc([0 90], [0 2000], bottom_loss(R_vals));
xlabel('Angle')
ylabel('Frequency')
title('Bottom loss')
colorbar;
saveas(gcf, 'Figure 16.png')

figure(18);

% fourier transoform for arrival times
R_by_freq = R_vals(:,20);
R_by_freq(1) = 0;

plot(abs(R_by_freq));
inverse_fft = ifft(abs(R_by_freq), length(R_by_freq) * 4);
fft_of_signal = fft((R_by_freq), length(R_by_freq) * 5);
figure(17); clf;
plot(abs(inverse_fft))
xlabel(' "Time" ')
ylabel('|R|')
title('Inverse fourier transform of R(frequency) in 3 layered medium')

figure(19); clf;
plot(abs(fft_of_signal));
title("FFT of R(frequency) in 3 layered medium)");
% saveas(gcf, 'Figure 17.png')

%% 

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

% mud over sand base
c = [water_speed, mud_speed, sand_speed];
rho = [water_density, mud_density, sand_density];
h_arr = 2:6:62; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(h_arr), length(angle_arr), length(freq_arr)); 


for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        for h_index = 1:length(h_arr)
            R_array = rayleigh_strat(rho,c,h_arr(h_index),angle_arr(angle_index),freq_arr(freq_index));
            R_arr(h_index, angle_index, freq_index) = R_array(1);
        end
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
% figure(18); clf;
for ind = 1:length(h_arr)
    figure(ind)
    surf(X, Y * 180/pi, abs(squeeze(R_arr(ind,:,:))));
end
surf(X, Y * 180/pi, abs(squeeze(R_arr(2,:,:))))

title('|R| for mud-sand')
xlabel('Frequency')
ylabel('Angle')

%% 

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

n_1 = water_speed/sand_speed;
m_1 = sand_density/water_density;
brewster_angle = asin(sqrt((m_1^2 - n_1^2)/(m_1^2 - 1)));

% sand-mud-sand base
c = [water_speed, sand_speed mud_speed, sand_speed];
rho = [water_density, sand_density, mud_density, sand_density];
h = [600,3]; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(angle_arr), length(freq_arr)); 

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq_arr(freq_index));
        R_arr(angle_index, freq_index) = R_array(1);      
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(19); clf;

surf(X, Y * 180/pi, abs(R_arr))

title('|R| for sand-mud-sand')
xlabel('Frequency')
ylabel('Angle')

figure(20); clf;

surf(X, Y * 180/pi, angle(R_arr))

title("Phase for sand-mud-sand")
xlabel("frequency")
ylabel("Angle")

%% 

%water-sand-mud

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

n_1 = water_speed/sand_speed;
m_1 = sand_density/water_density;
brewster_angle = asin(sqrt((m_1^2 - n_1^2)/(m_1^2 - 1)));

% mud over sand base
c = [water_speed, sand_speed, mud_speed];
rho = [water_density, sand_density, mud_density];
h = 20; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(angle_arr), length(freq_arr)); 

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq_arr(freq_index));
        R_arr(angle_index, freq_index) = R_array(1);      
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(21); clf;

surf(X, Y * 180/pi, abs(R_arr))

title('|R| for sand-mud')
xlabel('Frequency')
ylabel('Angle')

figure(22); clf;

surf(X, Y * 180/pi, angle(R_arr))

title("Phase for sand-mud")
xlabel("frequency")
ylabel("Angle")

%% 

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

% mud over sand base
c = [water_speed, sand_speed, water_speed];
rho = [water_density, sand_density, water_density];
h_arr = 5:5:20; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(h_arr), length(angle_arr), length(freq_arr)); 


for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        for h_index = 1:length(h_arr)
            R_array = rayleigh_strat(rho,c,h_arr(h_index),angle_arr(angle_index),freq_arr(freq_index));
            R_arr(h_index, angle_index, freq_index) = R_array(1);
        end
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
% figure(18); clf;
for ind = 1:length(h_arr)
    figure(ind)
    surf(X, Y * 180/pi, abs(squeeze(R_arr(ind,:,:))));
    title(h_arr(ind))
    xlabel('Frequency')
    ylabel('Angle')

end
%surf(X, Y * 180/pi, abs(squeeze(R_arr(2,:,:))))

