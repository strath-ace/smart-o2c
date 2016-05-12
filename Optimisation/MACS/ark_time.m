close all
clear all
clc

load ark_test;

ark_size = [10;20;40;80;160;320;640;1280;2560];

times = zeros(size(ark_size,1),2);

for i=1:length(times)
    i/length(i)*100
   
    tic
    ark1 = arch_shrk2(memory,1,3,ark_size(i));
    times(i,1) = toc;
    
    tic
    ark2 = arch_shrk3(memory,1,3,ark_size(i));
    times(i,2) = toc;
    
end

plot(ark_size,times(:,1),'bo',ark_size,times(:,2),'ro')