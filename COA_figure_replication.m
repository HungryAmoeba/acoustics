%% 
% try and get the same plots as page 49 of CoA

path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')

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


% fourier transoform for arrival times
R_by_freq = R_vals(:,20);
R_by_freq(1) = 0;
inverse_fft = ifft(abs(R_by_freq));
figure(17); clf;
plot(abs(inverse_fft))
xlabel('Time')
ylabel('|R|')
title('Inverse fourier transform of R(frequency) in 3 layered medium')
% saveas(gcf, 'Figure 17.png')
