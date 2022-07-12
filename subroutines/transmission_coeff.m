function T = transmission_coeff(m,n,theta)
T = 2*m*cos(theta)/(m * cos(theta) + sqrt(n^2 - sin(theta)^2));
end