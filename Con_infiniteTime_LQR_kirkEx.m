%Author: Geronimo Macias
%last updated on: 10/20/2020
%%% solves an example in a cont. Infinite Time Horizon %%%


clc;clear all;close all;

% note: in the cont. infinite-time method we assume that we are in
% steady-state so we only need to calculate one gain K and apply it for all
% time. In contrast with the cont. finite-time method we have to calculate
% a gain K each time step. 
% ==> much simpler to solve!

%--------------------------------------------------------------------------
% example 3.10-1 from Optimal Control Theory an Intro by Donald E. Kirk
Q = .25*eye(2);
R = 0.05;
A = [-0.0026 0.0539;-0.1078 0.1591];
B = [0.0013;0.0539];
% s.t. inital conditions
X0=[2;1]; 
%--------------------------------------------------------------------------
% time interval
dt = 0.001;
t = 0:dt:100;
%--------------------------------------------------------------------------
% can also use [Kopt,P,eigenValues]=lqr(A,B,Q,R)
P = care(A,B,Q,R);
K_opt = inv(R)*B'*P;
%--------------------------------------------------------------------------
% calling the governing equation to be solved via ode45 with the given
% parameters
[t, X_opt] = ode45(@(t,x) sysFun(t,x,K_opt,A,B),t,X0);
%--------------------------------------------------------------------------
% discarded
% creating anonymous function sys_fun
% sys_fun = @(t,x)[x(2);-K_opt*x]; %(?) system of decoupled first-order ODEs 
% [t, X_opt] = ode45(sys_fun,t,X0);
%--------------------------------------------------------------------------
figure(1);
hold on
plot(t,-K_opt*X_opt','-k','LineWidth',1.5)
hold off
ylabel('control')
xlabel('time')
legend('opt. control')
figure(2);
hold on
plot(t,X_opt(:,1),'-r','LineWidth',1.5)
plot(t,X_opt(:,2),'-b','LineWidth',1.5)
ylabel('states')
xlabel('time')
legend('x(t)','(v(t)')

% governing equation 
function xdot = sysFun(t,x,K,A,B)
    u=-K*x;
    xdot = A*x+B*u;
end
