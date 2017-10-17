function [r,dr] = rFileRK4(t,in2,in3,dt)
%RFILERK4
%    [R,DR] = RFILERK4(T,IN2,IN3,DT)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    17-Oct-2017 11:00:16

K1 = in3(1,:);
K2 = in3(2,:);
x1 = in2(1,:);
x2 = in2(2,:);
t2 = K1.*dt;
t3 = t2+x1;
t4 = K2.*dt;
t5 = t4+x2;
t6 = t3.^2;
t7 = t6.*5.0;
t8 = t7-5.0;
r = [K1-x2-K2.*dt;K2+t2+x1+t5.*t8];
if nargout > 1
    dr = reshape([1.0,dt+dt.*t3.*t5.*1.0e1,-dt,dt.*t8+1.0],[2,2]);
end
