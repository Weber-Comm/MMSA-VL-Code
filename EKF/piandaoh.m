function [H] = piandaoh(theta,theta0,v,f,ka,beta,Nt,Nr)
%H 求测量模型h的雅各比矩阵 theta0是预测值
%   此处显示详细说明
H=zeros(Nr+2,4);
A = ones(Nt,1);
An = ones(Nt,1);
B = ones(Nr,1);
ita_theta = ones(Nr,1);
for k = 1:Nt
    A(k,1) = sqrt(1/Nt) * exp(-pi*cos(theta)*(k-1)*1i);
end

for k = 1:Nt
    An(k,1) = sqrt(1/Nt) * exp(-pi*cos(theta0)*(k-1)*1i);
end

for k = 1:Nr
    B(k,1) = sqrt(1/Nr) * exp(-pi*cos(theta)*(k-1)*1i);
end

for k = 1:Nr
    temp=0;
    for j = 1:Nt
        temp = temp - exp(-1i*pi*((j-1)*cos(theta0)-(j-k)*cos(theta)))*1i*pi*(j-k)*sin(theta); 
    end
    ita_theta(k,1) = (beta/sqrt(Nt))*temp;
end
H(1:Nr,1)=ita_theta;
H(1:Nr,4)= ka*B*A'*An;
H(Nr+1,2)=2/(3*10^8);
H(Nr+2,1)= - 2*v*sin(theta)/(3*10^8);
H(Nr+2,3)= 2*f*cos(theta)/(3*10^8);

end

