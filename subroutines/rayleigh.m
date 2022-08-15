function R = rayleigh(m,n,theta, k_arr)

%This function calculates the reflection coefficient 
%given index of refraction n = k1/k = c/c1 and 
%m=rho_1/rho

if ~exist('k_arr','var')
 % k_arr parameter does not exist, so default it to something
  R = (m * cos(theta) - sqrt(n^2 - sin(theta).^2))./(m * cos(theta) + sqrt(n^2 - sin(theta).^2));

end

if exist('k_arr', 'var')
    assert(length(k_arr) == 2)

    % k_arr parameter has been given, so use it in the calculation
    
    %assume that k_arr is in order 
    R = (m .* k_arr(1) .* cos(theta) - sqrt(k_arr(2).^2 - k_arr(1).^2.* sin(theta).^2))./ ...
        (m * k_arr(1) .* cos(theta) + sqrt(k_arr(2)^2 - k_arr(1)^2* sin(theta).^2));
end

end
