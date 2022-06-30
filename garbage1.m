%为了验证一个函数的凹凸性设置的
clear all;
close all;
clc
x=-5:0.01:5;
y=-5:0.01:5;
[X,Y] = meshgrid(x,y);
Z=sqrt(X.^2+Y.^2);
final=1-abs(X./Z);
mesh(X,Y,final);