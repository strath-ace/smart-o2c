close all
clear all
clc

nsamp = 100:100:5000;

for i = 1:length(nsamp)
   
    A = rand(nsamp(i),2);
    tic
    dom = dominance(A,0);
    elapsedTime = toc;
    tic
    dom2 = dominance_new(A,0);
    elapsedTime2 = toc;
    all(dom==dom2)
    plot(nsamp(i),elapsedTime/elapsedTime2,'kx')
    hold on
    
end