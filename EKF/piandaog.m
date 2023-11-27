function [G] = piandaog(v,d,theta,beta,deltat)
%UNTITLED5 此处显示有关此函数的摘要
%   此处显示详细说明
G=ones(4,4);
G(1,1)= 1 + v*deltat*cos(theta)/d;
G(1,2) = -v*deltat*sin(theta)/(d^2);
G(1,3) = deltat*sin(theta)/d;
G(1,4) = 0;
G(2,1) = v*deltat*sin(theta);
G(2,2) = 1;
G(2,3) = -deltat*cos(theta);
G(2,4) = 0;
G(3,1) = 0;
G(3,2) = 0;
G(3,3) = 1;
G(3,4) = 0;
G(4,1) = - beta*v*deltat*sin(theta)/d;
G(4,2) = - beta*v*deltat*cos(theta)/(d^2);
G(4,3) = beta*deltat*cos(theta)/d;
G(4,4) = 1 + v*deltat*cos(theta)/d;
end

