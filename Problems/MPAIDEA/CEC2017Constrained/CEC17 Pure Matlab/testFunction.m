clear all;
clc
format long e;
x =ones(30,2);
global initial_flag;
for i=1:28
   i 
   initial_flag = 0;
   [f,g,h]=CEC2017(x',i)
end
