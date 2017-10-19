clc; clear all
%% Create function
n = 6;
p = sym('p',[3,1]);
v = sym('v',[3,1]);
x  = [p;v];
syms t g z m dt real

f = [x(4:6,1);
     -g*[0;0;1] - (1/m)*z*x(1:3)];

matlabFunction(f,'File','daeFunc','Vars',{t,x,z,g,m});
%%
s = 2;
K = sym('K',[n,s]);
x = sym('x',[n,1]);

A = [1/4, 1/4-sqrt(3)/6;
    1/4+sqrt(3)/6, 1/4];
b = [1/2, 1/2];

f =@(t,x,z,g,m) daeFunc(t,x,z,g,m);

for i = 1:s
    r(n*i-(n-1):n*i,1) = K(:,i) - f(t + dt,x + (dt*A(i,:)*K.').',z,g,m);
end
dr = jacobian(reshape(r,[s*6,1]),reshape(K,[s*6,1]));
matlabFunction(r,dr, 'file', 'rFileIRK1','vars',{t,x,K,dt,g,m,z});

xNext = x + (dt*b*K.').';
matlabFunction(xNext, 'file', 'xNextFileIRK1','vars',{t,x,K,dt,g,m,z});

clear xNext K x t dt r dr g m z dp dv dz
%% Simulation parameters
tf = 2;
dt = 0.1;
x0 = [1;0;0;0;1;0];

%% Simulation

[tIRK4,xIRK4] = IRK4(dt,tf,x0,n,s); % IRK4(dt,tf,x0,n,s)
plot3(xIRK4(1,:),xIRK4(2,:),xIRK4(3,:)),grid on, axis equal



