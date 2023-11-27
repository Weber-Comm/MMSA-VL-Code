function [h] = h(v, theta, theta0, f, d, beta, ka, Nt, Nr)
%h 用于获得测量值矩阵h
h=zeros(Nr+2,1);
A = ones(Nt,1);
An = ones(Nt,1);
B = ones(Nr,1);
for k = 1:Nt
    A(k,1) = sqrt(1/Nt) * exp(-pi*cos(theta)*(k-1)*1i);
end

for k = 1:Nt
    An(k,1) = sqrt(1/Nt) * exp(-pi*cos(theta0)*(k-1)*1i);
end

for k = 1:Nr
    B(k,1) = sqrt(1/Nr) * exp(-pi*cos(theta)*(k-1)*1i);
end

nz = 1 * ones(Nr,1) .* randn(1) * 10^(-6);%theta的噪声
nt = 6.7 * 10^(-6) * randn(1); %tau的噪声
nf = 0.1 * 10^(-6) * randn(1);
r = ka * beta * B * (A') * An + nz;
h(1:Nr,1) = r;
h(Nr+1,1) = 2*d/(3*10^8) + nt;
h(Nr+2,1) =  2 * v * cos(theta) * f/ (3*10^8) + nf;
end

