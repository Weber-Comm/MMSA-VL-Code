function [theta, d, v, beta] = g(theta0, d0, v0, beta0, deltat)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
nthata = 0.001*randn(1);  % theta������  
nd = 0.002*randn(1);  %d������
nv = 0.002*randn(1); %v������
nbeta = 0.1*randn(1); %beta������
theta = theta0 + (1/d0)*v0*deltat*sin(theta0) + nthata;
d = d0 - v0*deltat*cos(theta0)+ nd;
v = v0 + nv;
beta = beta0*(1+(1/d0)*v0*deltat*cos(theta0)) + nbeta;

end

