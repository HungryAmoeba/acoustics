function BL = bottom_loss(R)

%Calculate bottom loss given R

    BL = -20 * log10(abs(R));
end