function f = daeFunc(t,in2,z,g,m)
%DAEFUNC
%    F = DAEFUNC(T,IN2,Z,G,M)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    17-Oct-2017 14:46:30

p1 = in2(1,:);
p2 = in2(2,:);
p3 = in2(3,:);
v1 = in2(4,:);
v2 = in2(5,:);
v3 = in2(6,:);
t2 = 1.0./m;
f = [v1;v2;v3;-p1.*t2.*z;-p2.*t2.*z;-g-p3.*t2.*z];