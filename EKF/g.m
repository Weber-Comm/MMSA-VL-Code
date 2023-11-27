function [theta, d, v, beta] = g(theta0, d0, v0, beta0, deltat)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
nthata = 0.001*randn(1);  % theta的噪声  
nd = 0.002*randn(1);  %d的噪声
nv = 0.002*randn(1); %v的噪声
nbeta = 0.1*randn(1); %beta的噪声
theta = theta0 + (1/d0)*v0*deltat*sin(theta0) + nthata;
d = d0 - v0*deltat*cos(theta0)+ nd;
v = v0 + nv;
beta = beta0*(1+(1/d0)*v0*deltat*cos(theta0)) + nbeta;

end

