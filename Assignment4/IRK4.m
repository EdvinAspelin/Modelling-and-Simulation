function [t,x] = IRK4(dt,tf,x0,n,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Nsteps = tf/dt;
t = [0:dt:tf];
x = [x0,zeros(n,Nsteps-1)];
K = zeros(2*s,1);

for k = 1:Nsteps
    %Newton iteration
    iter = true;
    alpha = 1;
    niter = 0;
    while iter
        [r,dr] = rFileIRK1(t,x(:,k),K,dt); % (t,x(:,k),K,dt)
        K = K - dr\r;
        norm(r);
        if norm(r) < 1e-5
            iter = false;
        end
    end
    x(:,k+1) = xNextFileIRK1(t,x(:,k),K,dt);
end
end

