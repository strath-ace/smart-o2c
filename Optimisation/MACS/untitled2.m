close all
clear all
clc

x = 1:99;
y = zeros(length(x),1);

for i=1:length(x)
    y(i) = nchoosek(length(x),i);
end

plot(x,y,'bo')
hold on
plot(x(round(2*length(x)/3)),y(round(2*length(y)/3)),'r*');