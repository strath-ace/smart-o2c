function Ec = EccAnomal(ni,ecc)

Ec = 2*atan2(sqrt((1-ecc)/(1+ecc))*sin(ni/2),cos(ni/2));
