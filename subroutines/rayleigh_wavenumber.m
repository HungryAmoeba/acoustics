function R = rayleigh_wavenumber(k_x, m, k, k_1)

numerator = m .* sqrt(k.^2 - k_x.^2) - sqrt(k_1.^2 - k_x.^2);
denominator =  m .* sqrt(k.^2 - k_x.^2) + sqrt(k_1.^2 - k_x.^2);
R = numerator./denominator;

end

