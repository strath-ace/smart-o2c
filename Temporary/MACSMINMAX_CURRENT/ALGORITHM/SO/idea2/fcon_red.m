function  [Cf,f]=fcon_red(fworst,C)
%
%  [Cf,f]=fcon_red(fworst,C)
%  
% augmented fitness function for constraint problems
%
%  INPUT
%            fin         : input fitness
%            fworst  :  reference fitness
%            C          : constraint vector
%
%  OUTPUT
%            f          : augmented fitness
%            Cf        : scalar function for constraint violation
%
% (c) Massimiliano Vasile 2003
%
Cf=0;
for i=1:length(C)
    Cf=Cf+max([0 C(i)]);
end
f=fworst+Cf;
return