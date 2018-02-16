function t12 = ToFhyp(F1,F2,amax,ecc)

amax = min(amax,0);

t12 = sqrt(-amax^3)*(ecc*(sinh(F2)-sinh(F1)) - F2 + F1);

