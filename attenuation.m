% look at attenuation effects
path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')

clear all;
 
close all;
%%
%Example from Frisk
c=1500;
rho=1;
rho_1 = 1.8;
c1 = 1800;
n = c/c1;
m = rho_1/rho;

freq = 300;
% if attenuation is desired
lambda = [c, c1]./(freq);
beta = [.3, .1]; % typical is .1 - .5
alpha = beta./lambda * log(10)/20; % imaginary part of k
k = 2 * pi * freq ./ [c, c1] + 1i * alpha;

%c_complex = 2 * pi * freq./(k + 1i * alpha);

theta = 0:pi/1000:pi/2;
y = (rayleigh(m,n,theta));
y_attenuated = rayleigh(m, c/c1, theta, k);
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

%%
%Look at attenuations effect on the reflection coefficient calculations


