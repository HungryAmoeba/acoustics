function R = rayleigh(m,n,theta)

%This function calculates the reflection coefficient 
%given index of refraction n = k1/k = c/c1 and 
%m=rho_1/rho

R = (m * cos(theta) - sqrt(n^2 - sin(theta)^2))/(m * cos(theta) + sqrt(n^2 - sin(theta)^2));
end