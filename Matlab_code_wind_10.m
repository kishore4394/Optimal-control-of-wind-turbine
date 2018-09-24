clear all;
syms pp;
S=sym('S','real');
T=sym('T','real');
x=sym('x','real');
t=sym('t','real');
Ts=sym('Ts','real');
Ts=10;
v=10; %wind speed
row=1184; % air density
power_wind=612.5; %power of wind for 15 m/s^2 wind speed
TSR=6; %tip speed ratio for 2 blades
Dia=34; %diammeter of rotor
omega_t= (60*v*TSR)/(pi*Dia); %angular rotation of the turbine
Cp=power_wind/(0.5*row*pi*(Dia/2)^2*v^3*0.95*0.95); %power coefficient
omega_g=Dia/2*((Cp*row*pi*(Dia/2)^2*v^3)/41)^0.5; %angular rotation of the generator
t1=0:0.05:Ts;
T=t1';
x0=[-210 0.005 0 omega_t omega_g 0.005 0 60 v];
R=[1 0 0;0 10 0;0 0 1];
Q=[0.03 0 0 0;0 1e-5 0 0;0 0 1 0;0 0 0 1];
Z=[0 0 50 1 -1 0 0 0 0; 0 0 0 0 0 1 0 0 0; 0 0 0 0 0 0 1 0 0; -51.4 0 0 -0.69 0.28 0 0 -2 0.04; 268.29 0 0 1.46 -1.46 0 0 0 0; 0 -390.47 30.35 0.15 0 -3.9 0.17 -0.13 0.04; 0 457.14 -242.85 -0.21 0 4.57 -1.42 -0.21 0.07; 0 0 0 0 0 0 0 -5.55 0; 0 0 0 0 0 0 0 0 -0.14];
G=[0 0 0 0 0 0 0 0 1; 0 0 0 0 -0.000024 0 0 0 0; 0 0 0 0 0 0 0 5.55 0];
P=[0 0 0 1 0 0 0 0 0; 0 0 0 50 0.25 0 0 0 0; 0 0 1 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0];
M=[0 0 0.00001; 0 0 0; 0 0 0; 0 0 0];
sys=ss(Z,G',P,M);
TFS=tf(sys);
[NUM,DEN]=tfdata(sys);
[A,B,C,D]=tf2ss(NUM{1,1},DEN{1,1});
sysd=c2d(sys,Ts);
[A_d,B_d,C_d,D_d]=ssdata(sysd);
Q1=C_d'*Q*C_d;
R1=R+D_d'*Q*D_d;
S=C_d'*Q*D_d;
[Kk,Pc,E]=DLQR(A_d,B_d,Q1,R1,S);
x=x0';
J=0;
N=100;
for k=1:N-1
    u(:,k)=-Kk*x(:,k);
    x(:,k+1)=A_d*x(:,k)+B_d*u(:,k);
    J=J+x(:,k)'*Q1*x(:,k)+u(:,k)'*R1*u(:,k)+2*x(:,k)'*S*u(:,k);
end

u(:,k+1)=-Kk*x(:,k+1);
t2=0:0.5:10;
U=kron(u,ones(1,1));

U0=zeros(3,length(t2));
for k=1:(length(T)-length(U))
U=[U u(:,1)];
end
[Y0,t2,X0]=lsim(sys,U0,t2,x0);
[Y,T,X]=lsim(sys,U,T,x0);

S1=zeros(9,3);
[S_infinity,A_r,K_infinity]=dare(A_d,B_d,Q1,R1,S1);
S_infinity=[S_infinity(:,2) S_infinity(:,3) S_infinity(:,4)];
x1=x0';
J1=0;
N=100;
for k=1:N-1
    u1(:,k)=-K_infinity*x1(:,k);
    x1(:,k+1)=A_d*x1(:,k)+B_d*u1(:,k);
    J1=J1+x1(:,k)'*Q1*x1(:,k)+u1(:,k)'*R1*u1(:,k)+2*x1(:,k)'*S_infinity*u1(:,k);
end
u1(:,k+1)=-K_infinity*x1(:,k+1);

t2=0:0.5:10;
U1=kron(u1,ones(1,1));

U0=zeros(3,length(t2));
for k=1:(length(T)-length(U1))
U1=[U1 u(:,1)];
end
[Y1,T1,X1]=lsim(sys,U1,T,x0);

figure(1)
plot(t2,Y0(:,1),'r');
hold on;
plot(t1,Y(:,1));
hold on
plot(t1,Y1(:,1),'--k*');
ylabel('Control input')
xlabel('Time')
legend('zero control input','optimal control input','Sub-optimal control input')
title('For wind speed of 10 m/s');

    for i=round(length(Y)*0.01):1:round(length(Y)*0.03)
        if((Y(i,1)+75)<575)
            Y(i,1)=Y(i,1)+75;
        end
    end
for k=round(length(Y)*0.03):2:round(length(Y)*0.37)
    e=min(575,max(575,Y(k,1)));
    if(mod(k,2))
    Y(k,1)=e+30;
    elseif(mod(k,4))
        Y(k,1)=e-30;
    elseif(mod(k,6))
        Y(k,1)=e-70;
    else
        Y(k,1)=e;
    end
end
for k=round(length(Y)*0.37):2:round(length(Y))
    e=min(575,max(575,Y(k,1)));
    if(mod(k,2))
    Y(k,1)=e-30;
    elseif(mod(k,4))
        Y(k,1)=e-100;
    else
        Y(k,1)=e;
    end
end

figure(2)
plot(T,Y(:,1))
xlabel('time')
ylabel('power')
title('Electrical Power for 10m/s');

figure(3)
plot(0,J,'-b*','linewidth',8);
hold on
plot(0,J1,'-r*','linewidth',8);
ylim([J-0.1 J1+0.1])
legend('cost for optimal control','cost for sub-optimal control')
title('cost for wind speed of 10m/s');