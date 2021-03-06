function xNext = xNextFileIRK1(t,in2,in3,dt,g,m,z)
%XNEXTFILEIRK1
%    XNEXT = XNEXTFILEIRK1(T,IN2,IN3,DT,G,M,Z)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    17-Oct-2017 14:46:32

K1_1 = in3(1);
K1_2 = in3(7);
K2_1 = in3(2);
K2_2 = in3(8);
K3_1 = in3(3);
K3_2 = in3(9);
K4_1 = in3(4);
K4_2 = in3(10);
K5_1 = in3(5);
K5_2 = in3(11);
K6_1 = in3(6);
K6_2 = in3(12);
x1 = in2(1,:);
x2 = in2(2,:);
x3 = in2(3,:);
x4 = in2(4,:);
x5 = in2(5,:);
x6 = in2(6,:);
xNext = [x1+K1_1.*dt.*(1.0./2.0)+K1_2.*dt.*(1.0./2.0);x2+K2_1.*dt.*(1.0./2.0)+K2_2.*dt.*(1.0./2.0);x3+K3_1.*dt.*(1.0./2.0)+K3_2.*dt.*(1.0./2.0);x4+K4_1.*dt.*(1.0./2.0)+K4_2.*dt.*(1.0./2.0);x5+K5_1.*dt.*(1.0./2.0)+K5_2.*dt.*(1.0./2.0);x6+K6_1.*dt.*(1.0./2.0)+K6_2.*dt.*(1.0./2.0)];
