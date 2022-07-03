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