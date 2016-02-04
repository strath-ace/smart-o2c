close all
clear all
clc

input.population = zeros(3,3,4);
input.pop(:,:,1) = rand(3);
input.pop(:,:,2) = rand(3);
input.pop(:,:,3) = rand(3);
input.vlb = [-3 -3 -3];
input.vub = [3 3 3];

options = [];

%options.nFeValMax = 100;

memories = optimise(@testfun,input,options);