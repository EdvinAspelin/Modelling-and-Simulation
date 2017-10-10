clc;clear;
close all

%% Defining the function:
lambda = -2;
syms x;
f = @(t,x)(lambda*x);
matlabFunction(f,'File','myFunc','Vars',x);
u = 0;

%% Simulation parameters:
tFinal = 2;
dtEuler = 0.1;
dtRK2 = 0.1;
dtRK4 = 0.1;

x0 = 1;
btEuler = [0 0;
         0 1];
btRK2 = [0 0 0;
         1/2 1/2 0;
         0 0 1];
btRK4 = [0 0 0 0 0;
         1/2 1/2 0 0 0;
         1/2 0 1/2 0 0;
         1 0 0 1 0;
         0 1/6 1/3 1/3 1/6];
%% The exact Solution:
tExact = 0:0.01:tFinal;
xExact = exp(lambda*tExact);

%% Explicit Euler Method:

[tEuler,xEuler] = rk(x0, dtEuler, tFinal, f, btEuler);
%% Runge Kutta Method of order 2:

[tRK2,xRK2] = rk(x0, dtRK2, tFinal, f, btRK2);

%% Runge Kutta Method of order 4:

[tRK4,xRK4] = rk(x0, dtRK4, tFinal, f, btRK4);
%% Compare the results:
plot(tExact,xExact  ,'marker','.','markersize',20)
hold on
plot(tEuler,xEuler  ,'marker','.','markersize',10)
plot(tRK2,xRK2      ,'marker','.','markersize',10)
plot(tRK4,xRK4      ,'marker','.','markersize',10)
xlabel('t');
ylabel('x');
legend('Exact','Euler Method','RK2','RK4')

%% Global error:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
