function hank = hankel_zero_order_first(Z) 

hank = besselj(0, Z) + 1i * bessely(0, Z);

end