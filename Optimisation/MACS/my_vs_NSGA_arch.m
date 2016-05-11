close all
clear all
clc

num_arch = 50;
num_fin = 10;

x = rand(num_arch,1);
memory = [[1:length(x)]' x zeros(size(x,1)) zeros(size(x,1))];
plot(memory(:,2),zeros(size(memory(:,2))),'bo')
hold on

mem2 = arch_shrk3(memory,1,1,num_fin);
mem2 = mem2(:,2);
mem2 = sort(mem2);
plot(mem2,ones(size(mem2)),'r*')

mem3 = arch_shrk4(memory,1,1,num_fin);
mem3 = mem3(:,2);
mem3 = sort(mem3);
plot(mem3,2*ones(size(mem3)),'k^')

str1 = ['Mean distance with energy method: ', num2str(mean(diff(mem2))), ', std dev : ', num2str(std(diff(mem2))) ];
str2 = ['Mean distance with energy method: ', num2str(mean(diff(mem3))), ', std dev : ', num2str(std(diff(mem3))) ];
    
sprintf(str1)
sprintf(str2)