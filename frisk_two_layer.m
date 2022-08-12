%% 
% n = k1/k = c/c1, m = rho_1/rho

% all angles are incident angles (follow Frisks' convention, where 0 is
% normal incident and pi/2 is horizontal)
path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')

clear all;
 
close all;
saving = 0;
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
if saving == 1
    saveas(gcf,'Water-Air Interface.png')
end
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
if saving == 1
    saveas(gcf,'Water-Bottom Interface Model 1.png')
end
%% 
% Reflection from a Higher Velocity, More Dense Half-Space (the water bottom 
% interface: model 2)

rho_1 = 1.8;
c1 = 1800;
n = c/c1;
m = rho_1/rho;

freq = 300;
% if attenuation is desired
lambda = [c, c1]./(freq);
beta = .3; % typical is .1 - .5
alpha = beta./lambda * log(10)/20; % imaginary part of k
k = 2 * pi * freq ./ [c, c1];
c_complex = 2 * pi * freq./(k + 1i * alpha);

theta = 0:pi/1000:pi/2;
y = (rayleigh(m,n,theta));
y_attenuated = rayleigh(m, c_complex(1)/c1, theta);
figure(3); clf;
plot(theta * 180/pi, abs(y),'LineWidth',2)
xlabel('\theta')
ylabel('|R|')
hold on;
plot(theta * 180/pi, abs(y_attenuated), 'LineWidth', 2)
yyaxis right

%plot(theta * 180/pi, angle(y), '--')
xlabel('\theta')

legend("|R|", "|R| attenuated",'Location','northwest', 'FontSize', 18)
title("Water-Bottom Interface Model 2")

if saving == 1
    saveas(gcf,'Water-Bottom Interface Model 2.png')
end