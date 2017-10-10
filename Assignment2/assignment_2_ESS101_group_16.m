%% Task a
clc, clear all
%------------------------%
% Declaring input and output data
load('input.mat');
load('output.mat');
N = length(u);
uest = u(1:N/2);        % Estimation values
yest = y(1:N/2);       %
uval = u(N/2+1:end);    % Validation values
yval = y(N/2+1:end);   %
t = linspace(1,2000,2000);
%------------------------%
% Start with defining predictors
% Let's build the matrix H:
H = zeros(N/2,3);           % The matrix has N/2 rows (time instants), 2 columns (2 parameters)
H(1,:) = [0 0 uest(1)];     % The rows then matches Theta-values
H(1,:) = [0 0 uest(2)];
for i=3:N/2
    H(i,1) = -yest(i-1);
    H(i,2) = -yest(i-2);
    H(i,3) = uest(i);
end

th = H'*H\H'*yest;  % Estimation of Theta
a1hat = th(1)       % Values of parameters
a2hat = th(2)       % corresponding to columns in H
b0hat = th(3)

%% 20a task b
yn = yval;  % Defining new copy
un  = uval; % --> easier to switch between validation and estimation
% We'll start with calculating the predication model
ypred = zeros(N/2,1);
for i = 3:N/2
    ypred(i) = -a1hat*yn(i-1) - a2hat*yn(i-2) + b0hat*un(i);
end
predERROR20a = yn-ypred;            % Calculating root mean square
predRMSE20a = rms(predERROR20a)     % of predicted model error
% Now lets do the simulation
ysim = zeros(N/2,1); % Create vector to store simulated output
for i=3:N/2
    ysim(i) = -a1hat*ysim(i-1)  -a2hat*ysim(i-2) +b0hat*un(i);
end
simERROR20a = yn-ysim;                 % Calculating root mean square
simRMSE20a  = rms(simERROR20a)         % of simulated model error
simCOV20a = cov(simERROR20a)

%% 20b task a

% Define predictors

% Let's build the matrix H:
H = zeros(N/2,4); % The matrix has N/2 rows (time instants), 2 columns (2 parameters)
H(1,:) = [0 0 0 uest(1)];
H(2,:) = [0 0 0 uest(2)];
for i=3:N/2
    H(i,1) = -yest(i-1);
    H(i,2) = -yest(i-2);
    H(i,3) = uest(i);
    H(i,4) = uest(i-1);
end

th = H'*H\H'*yest;
a1hat = th(1)
a2hat = th(2)
b0hat = th(3)
b1hat = th(4)

%% 20b task b

yn = yval;  % Defining new copy
un  = uval; % --> easier to switch between validation and estimation
% We'll start with calculating the predication model
ypred = zeros(N/2,1);
for i = 3:N/2
    ypred(i) = b0hat*uest(i)+b1hat*uest(i-1)-a1hat*yest(i-1)-a2hat*yest(i-2);
end
predERROR20b = yn-ypred;            % Calculating root mean square
predRMSE20b = rms(predERROR20a)     % of predicted model error
% Now lets do the simulation
ysim = zeros(N/2,1); % Create vector to store simulated output
for i=3:N/2
    ysim(i) = b0hat*u(i)+b1hat*u(i-1)-a1hat*yest(i-1)-a2hat*ysim(i-2);
end
simERROR20b = yn-ysim;                 % Calculating root mean square
simRMSE20b  = rms(simERROR20b)         % of simulated model error
simCOV20b = cov(simERROR20b)

%% 20c task a

% Define predictors

% Let's build the matrix H:
H = zeros(N/2,4); % The matrix has N/2 rows (time instants), 2 columns (2 parameters)
H(1,:) = [0 0 0 uest(1)];
H(2,:) = [0 0 0 uest(2)];
H(3,:) = [0 0 0 uest(3)];
for i=4:N/2
    H(i,1) = -yest(i-1);
    H(i,2) = -yest(i-2);
    H(i,3) = -yest(i-3);
    H(i,4) = uest(i-1);
end

th = H'*H\H'*yest;
a1hat = th(1)
a2hat = th(2)
a3hat = th(3)
b1hat = th(4)

%% 20c task b

yn = yval;  % Defining new copy
un  = uval; % --> easier to switch between validation and estimation
% We'll start with calculating the predication model
ypred = zeros(N/2,1);
for i = 4:N/2
    ypred(i) = b1hat*uest(i-1)-a1hat*yest(i-1)-a2hat*yest(i-2)-a3hat*yest(i-3);
end
predERROR20b = yn-ypred;            % Calculating root mean square
predRMSE20b = rms(predERROR20a)     % of predicted model error
% Now lets do the simulation
ysim = zeros(N/2,1); % Create vector to store simulated output
for i=4:N/2
    ysim(i) = b1hat*u(i-1)-a1hat*ysim(i-1)-a2hat*ysim(i-2)-a3hat*ysim(i-3);
end
simERROR20c = yn-ysim;                 % Calculating root mean square
simRMSE20c = rms(simERROR20c)         % of simulated model error
simCOV20c = cov(simERROR20c)

%ya(i) = b0*u(t)+e(t) %20a
%yb =@(t) b0*u(t)+b1*u(t-1)+e(t)-a1*y(t-1)-a2*y(t-2) %20b
%yc =@(t) b1*u(t)+e(t)-a1*y(t-1)-a2*y(t-2)-a3*y(t-3) %20c




