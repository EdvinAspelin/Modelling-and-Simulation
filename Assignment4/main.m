clc;clear;
close all

%% Defining the function: (It is better to define this in a separate file)
n = 2;
x  = sym('x',[n,1]);
syms t real

f = [ x(2);
    5*(1-x(1)^2)*x(2) - x(1)];

matlabFunction(f, 'file', 'VanDerPol','vars',{t,x});
clear x t f

%% Defining r(xNext,x,u): (It is better to define this in a separate file)
s = 2;
K = sym('K',[2,s]);
x = sym('x',[2,1]);
syms t dt real

A = [1/4, 1/4-sqrt(3)/6;
    1/4+sqrt(3)/6, 1/4];
b = [1/2, 1/2];
f = @(t,x)VanDerPol(t,x);

for i = 1:s
    r(n*i-(n-1):n*i,1) = K(i,:).' - f(t + dt,x + (dt*A(i,:)*K.').');
end
dr = jacobian(reshape(r,[s*2,1]),reshape(K,[s*2,1]));
matlabFunction(r,dr, 'file', 'rFileIRK1','vars',{t,x,K,dt});

xNext = x + (dt*b*K).';
matlabFunction(xNext, 'file', 'xNextFileIRK1','vars',{t,x,K,dt});

clear xNext K x t dt r dr

%% Simulation parameters:
tf = 25;
dt = 0.1;
x0 = [1;0];

%% The Solution:
% solve using ode15s and store in t.ode15s and x.ode15s
[t.ode15s,x.ode15s] = ode15s(@VanDerPol, [0 tf], x0);

% Plot the results:
figure(1);clf
for subfig = 1:n
    subplot(n,1,subfig);hold on
    plot(t.ode15s,x.ode15s(:,subfig),'k','markersize',10)
end

disp(['Please press any key to continue...'])
pause
%% Simulation:
%Simulate using IRK1
[tIRK4,xIRK4] = IRK4(dt,tf,x0,n,s); % IRK4(dt,tf,x0,n,s)
% Plot the results:
% figure(2);clf
%% RK4
btRK4 = [0 0 0 0 0;
         1/2 1/2 0 0 0;
         1/2 0 1/2 0 0;
         1 0 0 1 0;
         0 1/6 1/3 1/3 1/6];
     
syms x y;
dx = y; dy = 5*(1-x^2)*y-x; % Diff equation to solve
matlabFunction([dx,dy],'File','odefunc', 'vars',{[x,y]});
[tVDP,xVDP] = rk(x0, dt, tf, @odefunc, btRK4, 2);

%% Plot result
for subfig = 1:n
    subplot(n,1,subfig);hold on
    plot(tIRK4,xIRK4(subfig,:),'r')
    plot(tVDP,xVDP(:,subfig),'b')
end
set(gca,'fontsize',20)
legend('Exact','IRK4','RK4')

%% Task 3
n = 3;
p  = sym('p',[n,1]);
v  = sym('v',[n,1]);
x = [p;v];
syms t real

dp = v;
dv = -g*[0;0;1] - (1/m)*z*p;
dz = -v.'*v\p.';

matlabFunction([dp,dv,dz],'File','daeFunc','Vars',{t,x})

s = 2;
K = sym('K',[2,s]);
x = sym('x',[2,1]);
syms t dt real

A = [1/4, 1/4-sqrt(3)/6;
    1/4+sqrt(3)/6, 1/4];
b = [1/2, 1/2];




