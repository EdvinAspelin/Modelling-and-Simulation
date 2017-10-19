function [t,x] = IRK4(dt,tf,x0,n,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nsteps = tf/dt;
t = [0:dt:tf];
x = [x0,zeros(n,Nsteps-1)];
K = zeros(n*s,1);

g = 9.82;
m = 1;
syms z;

for k = 1:Nsteps
    %Newton iteration
    iter = true;
    alpha = 1;
    niter = 0;
    while iter
        ztrue = solve(x(1:3,k)'*(-g*[0;0;1]-(1/m)*z*x(1:3,k))+x(4:6,k)'*x(4:6,k),z);
        [r,dr] = rFileIRK1(t,x(:,k),K,dt,g,m,ztrue); % t,x,K,dt,g,m,z
        K = K - dr\r;
        norm(r);
        if norm(r) < 1e-5
            iter = false;
        end
        t(k)
    end
    x(:,k+1) = xNextFileIRK1(t,x(:,k),K,dt);
end
end

