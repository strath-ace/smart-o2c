function t12 = ToFell(E1,E2,amax,ecc)

t12 = sqrt(amax^3)*(E2-E1 - ecc*(sin(E2)-sin(E1)));

