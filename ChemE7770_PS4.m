%% ChemE 7770 PS4 
clear all; close all; clc;
%% Problem 1
%% 1a)
% See the attached paper for calculation steps.

%% 1b). 
% global k_D k
% [result,fval,exit,output]=fsolve(@Signaling,A0,options);
% k=0.1;
k_D=0.1;k1=0.1;
for i=1:1000  
thetaB(i)=1/(1+k_D(i));
xstar(i)=((5.*thetaB(i).*k1+k1+1-5.*thetaB(i))-((5.*thetaB(i)-5.*thetaB(i).*k1-k1-1).^2-4.*(1-5.*thetaB(i)).*(5.*thetaB(i).*k1)).^(0.5))./(2.*(1-5.*thetaB(i)));
ystar(i)=((1+k1+10.*xstar(i).*k1-10.*xstar(i))-((10.*xstar(i)-10.*xstar(i).*k1-k1-1).^2-4.*(1-10.*xstar(i)).*(10.*xstar(i).*k1)).^(0.5))./(2.*(1-10.*xstar(i)));
k_D(i+1)=k_D(i)+0.1;
end
kd=k_D(1:length(k_D)-1);
figure (1)
plot(1./kd,thetaB,1./kd,xstar,1./kd,ystar);
legend('thetaB','xstar','ystar');
xlabel('1/k_D');
title('1b). k=0.1');
% k=10
k2=10;
k_D=0.1;
for i=1:1000 
thetaB2(i)=1/(1+k_D(i));
xstar2(i)=((5.*thetaB2(i).*k2+k2+1-5.*thetaB2(i))-((5.*thetaB2(i)-5.*thetaB2(i).*k2-k2-1).^2-4.*(1-5.*thetaB2(i)).*(5.*thetaB2(i).*k2)).^(0.5))./(2.*(1-5.*thetaB2(i)));
ystar2(i)=((1+k2+10.*xstar2(i).*k2-10.*xstar2(i))-((10.*xstar2(i)-10.*xstar2(i).*k2-k2-1).^2-4.*(1-10.*xstar2(i)).*(10.*xstar2(i).*k2)).^(0.5))./(2.*(1-10.*xstar2(i)));
k_D(i+1)=k_D(i)+0.05;
end
kd=k_D(1:length(k_D)-1);
xstar2(79)=0.5;
ystar2(79)=0.842;
figure (2)
plot(1./kd,thetaB2,1./kd,xstar2,1./kd,ystar2);
legend('thetaB','xstar','ystar');
xlabel('1/k_D');
title('1b). k=10');

%% 1c) Hill Coeff
% See excel for detail. Hill coeff for 1/kd vs thetaB, xstar, ystar
% respectively:for k=0.1, n=1 3.35 7.08; c1/2=2.22 0.483 0.245;
% for k=10, n=1 1.03 1.057; c1/2=1 0.174 0.02

%% 1d) Change as kd changes
% k=0.1;
KD=1/0.1; k1=0.1;
thetaB1d=1/(1+KD);
xstar1d=((5.*thetaB1d.*k1+k1+1-5.*thetaB1d)-((5.*thetaB1d-5.*thetaB1d.*k1-k1-1).^2-4.*(1-5.*thetaB1d).*(5.*thetaB1d.*k1)).^(0.5))./(2.*(1-5.*thetaB1d));
ystar1d=((1+k1+10.*xstar1d.*k1-10.*xstar1d)-((10.*xstar1d-10.*xstar1d.*k1-k1-1).^2-4.*(1-10.*xstar1d).*(10.*xstar1d.*k1)).^(0.5))./(2.*(1-10.*xstar1d));
KD2=1/0.15;
thetaB1d2=1/(1+KD2);
xstar1d2=((5.*thetaB1d2.*k1+k1+1-5.*thetaB1d2)-((5.*thetaB1d2-5.*thetaB1d2.*k1-k1-1).^2-4.*(1-5.*thetaB1d2).*(5.*thetaB1d2.*k1)).^(0.5))./(2.*(1-5.*thetaB1d2));
ystar1d2=((1+k1+10.*xstar1d2.*k1-10.*xstar1d2)-((10.*xstar1d2-10.*xstar1d2.*k1-k1-1).^2-4.*(1-10.*xstar1d2).*(10.*xstar1d2.*k1)).^(0.5))./(2.*(1-10.*xstar1d2));
thetaB_ch=abs((thetaB1d-thetaB1d2)./thetaB1d).*100;
xstar_ch=abs((xstar1d-xstar1d2)./xstar1d).*100;
ystar_ch=abs((ystar1d-ystar1d2)./ystar1d).*100;
disp(['% change of thetaB, xstar, and ystar for k=0.1 respectivly: ',num2str(thetaB_ch),'%, ',num2str(xstar_ch),'%, and ',num2str(ystar_ch),'%']);

% k=10;
KD=1/0.1; k2=10;
thetaB1d=1/(1+KD);
xstar1d=((5.*thetaB1d.*k2+k2+1-5.*thetaB1d)-((5.*thetaB1d-5.*thetaB1d.*k2-k2-1).^2-4.*(1-5.*thetaB1d).*(5.*thetaB1d.*k2)).^(0.5))./(2.*(1-5.*thetaB1d));
ystar1d=((1+k2+10.*xstar1d.*k2-10.*xstar1d)-((10.*xstar1d-10.*xstar1d.*k2-k2-1).^2-4.*(1-10.*xstar1d).*(10.*xstar1d.*k2)).^(0.5))./(2.*(1-10.*xstar1d));
KD2=1/0.15;
thetaB1d2=1/(1+KD2);
xstar1d2=((5.*thetaB1d2.*k2+k2+1-5.*thetaB1d2)-((5.*thetaB1d2-5.*thetaB1d2.*k2-k2-1).^2-4.*(1-5.*thetaB1d2).*(5.*thetaB1d2.*k2)).^(0.5))./(2.*(1-5.*thetaB1d2));
ystar1d2=((1+k2+10.*xstar1d2.*k2-10.*xstar1d2)-((10.*xstar1d2-10.*xstar1d2.*k2-k2-1).^2-4.*(1-10.*xstar1d2).*(10.*xstar1d2.*k2)).^(0.5))./(2.*(1-10.*xstar1d2));
thetaB_ch2=abs((thetaB1d-thetaB1d2)./thetaB1d).*100;
xstar_ch2=abs((xstar1d-xstar1d2)./xstar1d).*100;
ystar_ch2=abs((ystar1d-ystar1d2)./ystar1d).*100;
disp(['% change of thetaB, xstar, and ystar for k=10 respectivly: ',num2str(thetaB_ch2),'%, ',num2str(xstar_ch2),'%, and ',num2str(ystar_ch2),'%']);

%% e). Importance of parameter tuning
% The output of xstar and ystar has big changes as input has a small changes 
% for k=0.1 while for k=10 the change is less. Thus, it's important to define 
% parameter that fits the model. In this case, if we want the output to be 
% very sensitive, we want k=0.1, on the other hand, if we want the output to 
% has less variation as input varies, we would want k=10.

%% Problem #2
%% 2a)
B0=[4 4 6]';
[Soln,fval,exit,output]=fsolve(@NoInhi,B0);
display(Soln);
%% 2b)
I_1=0.01:100:1000;
I_2=0.01:100:1000;
points=length(I_1);
Aconc=zeros(length(I_1),length(I_2));
for j=1:points
  for i=1:points
    syms A I1 I2
    I1=I_1(j);
    I2=I_2(i);
    eqn=A.*(1+25./((1+I1).*(5+A)-5.*A)+25./((1+I2).*(5+A)-5.*A))==100;
    sol=solve(eqn,A);
    sol=double(sol);
    sol=sol(sol>=0);
    [m,n]=(size(sol));
    if m>1
      val=min(sol);
    else
      val=sol;
    end
    Aconc(j,i)=val;
  end
end
figure (4)
contour3(log(I_1),log(I_2),Aconc,200);
title('3D Plot for K=5');
xlabel('[I_1]');
ylabel('[I_2]');
zlabel('[A]');
% 2c) The system is a OR logic gate.
% 2d)
I_1=0.01:100:1000;
I_2=0.01:100:1000;
points=length(I_1);
Aconc=zeros(length(I_1),length(I_2));
for j=1:points
  for i=1:points
    syms A I1 I2
    I1=I_1(j);
    I2=I_2(i);
    eqn=A.*(1+175/((1+I1).*(35+A)-5.*A)+175./((1+I2).*(35+A)-5.*A))==100;
    sol=solve(eqn,A);
    sol=double(sol);
    sol=sol(sol>=0);
    [m,n]=(size(sol));
    if m>1
      val=min(sol);
    else
      val=sol;
    end
    Aconc(j,i)=val;
  end
end
figure (5)
contour3(log(I_1),log(I_2),Aconc,200);
title('3D Plot for K=35');
xlabel('[I_1]');
ylabel('[I_2]');
zlabel('[A]');
% k=35 gives a slower respond compare to the ones from k=5; also from the
% plot, there are some "unstable" values at the bottle, not very uniform. Thus, this one is "fuzzy"
%% 2e)
% The choices of parameter is very important and it affects how sensitive
% the output with regrad to the input.
