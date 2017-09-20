clc
clear all
%------------% 1a
syms x y z phi theta L m1 m2 dx dy dz dtheta dphi g ddx ddy ddz ddphi ddtheta
timeDep = [x; % Time depending variables
    y;
    z;
    phi;
    theta;
    dx;
    dy;
    dz;
    dtheta;
    dphi];
dTimeDep = [dx;
    dy;
    dz;
    dphi;
    dtheta;
    ddx;
    ddy;
    ddz;
    ddtheta;
    ddphi];
%-----------%
p1 = [x;
    y;
    z];
dp1 = [dx;
    dy;
    dz];
ddp1 = [ddx;
    ddy;
    ddz];
q = [p1;
    theta;
    phi];
dq = [dp1;
    dtheta;
    dphi];
ddq = [ddp1;
    ddtheta;
    ddphi];
M =  [m1;
    m1;
    m1;
    m2;
    m2];
p2 = [L*sin(q(5))*cos(q(4))+q(1);
    L*sin(q(5))*sin(q(4))+q(2);
    L*cos(q(5))+q(3)];
jp2 = jacobian(p2, q);
dp2 = jp2*dq;

T2 = (1/2)*m2*(dp2.'*dp2);
V2 = m2*g*[0 0 1]*p2;
L2 = T2 - V2;
DL2=jacobian(L2,q)*dq;
DDL2=jacobian(DL2,timeDep)*dTimeDep;
Q = DDL2 - DL2;

b = diag(M)*ddq;
%% 1b
%------------%
syms x1 y1 z1 x2 y2 z2 ...
    phi theta Le m1 m2 ... 
    dx1 dy1 dz1 dx2 dy2 ...
    dz2 dtheta dphi g ddx1 ...
    ddy1 ddz1 ddx2 ddy2 ddz2 ...
    ddphi ddtheta z m1 m2 u ...
    ux uy uz ...
    real
%------------%
% Declare matrixes
p1 = [x1;
    y1;
    z1];
p2 = [x2;
    y2;
    z2];
dp1 = [dx1;
    dy1;
    dz1];
dp2 = [dx2;
    dy2;
    dz2];
ddp1 = [ddx1;
    ddy1;
    ddz1];
ddp2 = [ddx2;
    ddy2;
    ddz2];
q = [p1;
    p2];
dq = [dp1;
    dp2];
ddq = [ddp1;
    ddp2];
e = p1-p2;
de = jacobian(e,q)*dq;
C = (1/2)*((e.'*e)-Le^2);
m = [0;
    0;
    -m1*g;
    0;
    0;
    -m2*g];
M = [m1;
    m1;
    m1;
    m2;
    m2;
    m2];
u=[ux;uy;uz];
Q = u.'*jacobian(p1,q);
% Compute Lagrange
T = (1/2)*(m1*(dp1.'*dp1)+m2*(dp2.'*dp2));
V = g*(m1*p1(3)+m2*p2(3));
L = T - V - z*C;
jdC = jacobian(C, q); 
dC = jdC*dq;
jddC = jacobian(dC, [q;dq]);
ddC = jddC*[dq;ddq];

dL = jacobian(L, q);
ddL = jacobian(jacobian(L, dq), [q;dq])*[dq;ddq];
dyn = ddL - dL.' - Q.';

ddCz = M.^-1.*(m-z.*jdC.');
ddC2 = jddC(1:length(q))*dq + jddC(length(q)+1:length([q;dq]))*ddCz;
accv = M.^-1.*(z*[e;e]+m+Q.'); %ddq
ddCsolve = [e;e].'*accv+[de;de].'*dq;
ztrue = solve(ddCsolve,z);
%% 2
clc
%a
Phi = diag([M;0]) + [zeros(size(diag(M))), jdC.';[jdC,0]];
c = Phi*[ddq;z];
%b
ddqz = Phi^-1*c;
sddqz = simplify(ddqz);

%% 3
clc 
clear all
%--------------%
% Declare variables
syms l p0 m1 m2 m3 z C g Q Z1 Z2 Z3 Z4 Z5... % constants
    L p1 x0 y0 z0 x1 y1 z1 p2 x2 y2 z2 p3 x3 y3 z3 q ... % time dependent
    dL dp1 dx0 dy0 dz0 dx1 dy1 dz1 dp2 dx2 dy2 dz2 dp3 dx3 dy3 dz3 dq ...
    ddL ddx0 ddy0 ddz0 ddp1 ddx1 ddy1 ddz1 ddp2 ddx2 ddy2 ddz2 ddp3 ddx3 ddy3 ddz3 ddq ...
    real

Q = 0;
Z = [Z1;
    Z2;
    Z3;
    Z4;
    Z5];

p0 = [x0;
    y0;
    z0];
p1 = [x1;
    y1;
    z1];
p2 = [x2;
    y2;
    z2];
p3 = [x3;
    y3;
    z3];
dp0 = [dx0;
    dy0;
    dz0];
dp1 = [dx1;
    dy1;
    dz1];
dp2 = [dx2;
    dy2;
    dz2];
dp3 = [dx3;
    dy3;
    dz3];
ddp0 = [ddx0;
    ddy0;
    ddz0];
ddp1 = [ddx1;
    ddy1;
    ddz1];
ddp2 = [ddx2;
    ddy2;
    ddz2];
ddp3 = [ddx3;
    ddy3;
    ddz3];
q = [p1;
    p2;
    p3];
dq = [dp1;
    dp2;
    dp3];
ddq = [ddp1;
    ddp2;
    ddp3];
e1 = [(p0-p1).';
    (p1-p2).';
    (p1-p3).';
    (p2-p3).';
    [0 0 0]];
e2 = [[0 0 0];
    (p1-p2).';
    [0 0 0];
    (p1-p0).';
    (p3-p2).'];
e3 = [[0 0 0];
    [0 0 0];
    (p1-p3).';
    (p0-p1).';
    (p2-p3).'];
de1 = [(dp0-dp1).';
    (dp1-dp2).';
    (dp1-dp3).';
    (dp2-dp3).';
    [0 0 0]];
de2 = [[0 0 0];
    (dp1-dp2).';
    [0 0 0];
    (dp1-dp0).';
    (dp3-dp2).'];
de3 = [[0 0 0];
    [0 0 0];
    (dp1-dp3).';
    (dp0-dp1).';
    (dp2-dp3).'];
e = [e1 e2 e3];
de =[de1 de2 de3];
%--------------%
% Compute Lagrange
G = [0 0 g];
M = [m1;
    m1;
    m1;
    m2;
    m2;
    m2;
    m3;
    m3;
    m3];
m = M.*[G G G].';
V = G*p1 + G*p2 + G*p3;
T = (m1/2)*(dp1.'*dp1) + (m2/2)*(dp2.'*dp2) + (m3/2)*(dp3.'*dp3);
La = T-V-z*C;
C =[(1/2)*((p1-p0).'*(p1-p0)-L^2);
    (1/2)*((p2-p1).'*(p2-p1)-l^2);
    (1/2)*((p3-p1).'*(p3-p1)-l^2);
    (p1-p0).'*(p3-p2);
    (1/2)*((p3-p2).'*(p3-p2)-(2*l)^2)];
%--------------%
% Solve C for z
jdC = jacobian(C, q); 
dC = jdC*dq;
jddC = jacobian(dC, [q;dq]);
ddC = jddC*[dq;ddq];

dLa = jacobian(La, q);
ddLa = jacobian(jacobian(La, dq), [q;dq])*[dq;ddq];
dyn = ddLa - dLa.' - Q.';

ddCz = M.^-1.*(m-(Z.'*jdC).');
ddC2 = jddC(1:length(q))*dq + jddC(length(q)+1:length([q;dq]))*ddCz;
accv = M.^-1.*((Z.'*e).'-m); %ddq

Mq = diag([M;zeros(min(size(jdC)),1)]) + [zeros(size(diag(M))), jdC.';[jdC,zeros(1,min(size(jdC.')))]];