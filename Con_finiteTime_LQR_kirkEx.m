%Author: Geronimo Macias
%last updated on: 10/20/2020
%%% solves an example in a cont. Finite Time Horizon %%%


clc;clear all;close all;
% continous time finite horizon LQR problem 
t_f = 100; % final time 
dt = 0.001; % resolution
P_N= 100*eye(2); % (or P_f) boundary condition (???)
%--------------------------------------------------------------------------
% example 3.10-1 from Optimal Control Theory an Intro by Donald E. Kirk
Q = .25*eye(2);
R = 0.05;
A = [-0.0026 0.0539;-0.1078 0.1591];
B = [0.0013;0.0539];
%--------------------------------------------------------------------------
N = t_f/dt+1; % number of iterations
%--------------------------------------------------------------------------
t_dis = 0:dt:t_f; % forwards in time
t_res = t_f:-dt:0; % backwards in time
%--------------------------------------------------------------------------
% calculating P_all_res (backwards in time)
P_all_res(:,length(t_res))= P_N(:); 
% we start at the final time and move to the intial time 
% in this case we assign P_N to be the first entry in the series of P
% matrices to follow 
for i = N:-1:2
    % reshape turns the size of the P_all_res matrix from 4x1 to 2x2 (size
    % of A)
    P =reshape(P_all_res(:,i),size(A));
    % calculating value of dPdt at a specific time instant
    dPdt = -(A.'*P + P*A - P*B*inv(R)*B.'*P + Q); 
    % employing euler method 
    P = P - dt*(dPdt);
    % assign P matrix size 4x1 into P_all_res()
    P_all_res(:,i-1)= P(:);
   
end
P_all_res = P_all_res'; % change dim. from 4x1 to 1x4
%--------------------------------------------------------------------------
% applying initial condition and calculating optimal control law (forwards
% in time)
X0=[2;1]; 
X_eul(:,1) = X0;
for i = 2:N
    P_eul = reshape(P_all_res(i-1,:),size(A));
    U_eul(:,i-1) = -inv(R)*B'*P_eul*X_eul(:,i-1);
    X_eul(:,i) = X_eul(:,i-1) + dt* (A*X_eul(:,i-1) + B*U_eul(:,i-1) );  
end
%--------------------------------------------------------------------------
figure(1);
hold on
plot(t_dis(1:end-1),U_eul,'--k','LineWidth',1.5)
ylabel('control')
xlabel('time')
legend('opt. control')

figure(2);
hold on
plot(t_dis,X_eul(1,:),'--r','LineWidth',1.5)
plot(t_dis,X_eul(2,:),'--b','LineWidth',1.5)
ylabel('states')
xlabel('time')
legend('x(t)','v(t)')