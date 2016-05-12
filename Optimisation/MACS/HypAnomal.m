function F = HypAnomal(ni,ecc)

qq = min(sqrt((ecc-1)/(1+ecc))*tan(ni/2),1-eps);
qq = max(qq,-1+eps);

F = 2*atanh(qq);
