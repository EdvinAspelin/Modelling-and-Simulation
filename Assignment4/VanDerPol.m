function f = VanDerPol(t,in2)
%VANDERPOL
%    F = VANDERPOL(T,IN2)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    17-Oct-2017 11:04:11

x1 = in2(1,:);
x2 = in2(2,:);
f = [x2;-x1-x2.*(x1.^2.*5.0-5.0)];
