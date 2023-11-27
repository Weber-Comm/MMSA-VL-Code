% EKF1
clear;
load("tank2to6_center_angle_gt.mat");
load("tank2to6_center_xy_gt.mat")

nResult = [];
rng(3);


tf=1500; %�������
Nt=8; %��������
Nr=8; %��������
fc=28*10^9; %Ƶ��
deltat=0.1; %֡��ʱ����
v=2.1; %�ٶ�
ka = sqrt(Nt*Nr);  %������������

xaxis=tank2to6_center_xy_gt(1,:);  %ground truth x��
yaxis=tank2to6_center_xy_gt(2,:);  %ground truth y��

xaxis(1200:1400) = xaxis(1200:1400) - 4; 

realangle=angleNormalized(atan2(yaxis,xaxis))-pi/2; %ʵ�ʽǶ�

G = 20; %gain
x = [realangle(1),sqrt(xaxis(1)^2 + (yaxis(1))^2),v,0.5+0.5i]';  %�����ĳ�ʼ�˶�״̬������sqrt(xaxis(1)^2 + (yaxis(1)-44)^2)���շ���֮���ʼ���룬��Ҫ��һ��
qs = [0.02^2 0.02^2 0.05^2 0.01^2];
Qs = diag(qs); %״̬����
qm = [ones(1,Nr)*10^(-2) 6.7*10^(-2) 2*10^(-2)];
Qm= diag(qm);  %����ֵ����
P =eye(4);
% x=zeros(1,tf);
Result=zeros(4,tf);   %ResultΪ���
% x(1,1)=0.1; 
Result(1:4,1)=x;  %����ֵ
z=zeros(1,tf);
xpre1=zeros(4,tf);
xpre2=zeros(4,tf);
%��beta��ʵֵ
betareal=zeros(tf,1);
betareal(1,1) = 0.5 + 0.5i;
for m=2:tf
    betareal(m,1) = betareal(m-1,1)*(1 + (1/xaxis(m-1))*v*deltat*cos(realangle(m-1)));
end
for k = 2 : tf 
    %%%%%%%%%%%%%%%%%%EKF��ʼ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [xpre1(1,k),xpre1(2,k),xpre1(3,k),xpre1(4,k)] = g(Result(1,k-1),Result(2,k-1),Result(3,k-1),Result(4,k-1),deltat);        %Ԥ��ֵ1
    [xpre2(1,k),xpre2(2,k),xpre2(3,k),xpre2(4,k)] = g(xpre1(1,k),xpre1(2,k),xpre1(3,k),xpre1(4,k),deltat);        %Ԥ��ֵ2
     %ÿʱÿ�̵Ĳ���ֵ
    y = zeros(Nr+2,1);
    A = ones(Nt,1);
    An = ones(Nt,1);
    B = ones(Nr,1);
    for kk = 1:Nt
        A(kk,1) = sqrt(1/Nt) * exp(-pi*cos(realangle(k))*(kk-1)*1i);
    end
    for kk = 1:Nt
        An(kk,1) = sqrt(1/Nt) * exp(-pi*cos(xpre1(1,k))*(kk-1)*1i);
    end
    for kk = 1:Nr
        B(kk,1) = sqrt(1/Nr) * exp(-pi*cos(realangle(k))*(kk-1)*1i);
    end
    y(1:Nr,1) = ka*betareal(k)*B*(A')*An+1*10^(-4)*randn(1);
    y(Nr+1,1) = 2*sqrt((xaxis(k))^2+(yaxis(k)-44)^2)/(3*10^8)+ 1*10^(-4)*rand(1);    %((xaxis(k))^2+(yaxis(k)-44)^2)�˴�������շ���֮�����Ծ��룬��Ҫ�޸�һ��
    y(Nr+2,1) = 2*v*cos(realangle(k))*fc/(3*10^8)+1*10^(-4)*rand(1);
    %���Ի�
    G = piandaog(Result(1,k-1),Result(2,k-1),Result(3,k-1),Result(4,k-1),deltat);   %״̬�ݻ����Ÿ��Ⱦ���   
    H = piandaoh(realangle(k),xpre1(1,k),xpre1(3,k),fc,ka,xpre1(4,k),Nt,Nr);   %����ֵ���Ÿ��Ⱦ���
    %MSE Matrix Prediction:
    PP=G*P*G'+Qs;
    %Kalman Gain Calculation:
    if sum(sum(1*(isnan(H*PP*H'+Qm))))~=0
        break
        NAN= [NAN, sim];
    end
    Kk=PP * H'* pinv(H*PP*H'+Qm);  
    %State Tracking:
    Result(:,k) = xpre1(:,k) + Kk*(y-h(v, realangle(k),xpre1(1,k),fc,sqrt((xaxis(k))^2+(yaxis(k)-44)^2),betareal(k),ka,Nt,Nr));   
    P=PP-Kk*H*PP;
end



%%



f = figure;
plot(realangle);
hold on
plot(real(Result(1,:)));



%%

writematrix([(1:1500).',real(Result(1,:)).'-pi/2],"EKF_DFRC_angle_est.csv");
