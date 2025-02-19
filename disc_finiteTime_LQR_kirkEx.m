%Author: Geronimo Macias
%last updated on: 10/24/2020
%%% solves an example in a discrete Infinite Time Horizon %%%

% (works!!!) on 10/24/2020

clc;clear all;close all; 
% discrete infinite time horizon LQR problem 
t_f = 200; % final time 
dt = 1; % resolution
N = t_f/dt+1; % number of iterations
S{N+1} = 0*eye(2); % S_{N+1} = Q_N boundary condition 
% K{N+1} = [0 0];
%--------------------------------------------------------------------------
% example 3.10-1 from Optimal Control Theory an Intro by Donald E. Kirk
% (discretized version)
Q = dt*[0.25 0;
        0 0.05];
R = dt*0.05;
A = dt*[-0.0026 0.0539;-0.1078 0.1591]+eye(2);
B = dt*[0.0013;0.0539];
%--------------------------------------------------------------------------
t_dis = 0:dt:t_f; % forwards in time
%--------------------------------------------------------------------------
for i = N:-1:1 % (backwards in time)
    % https://en.wikipedia.org/wiki/Linear%E2%80%93quadratic_regulator
    K(i,:) = -inv(R+B'*S{i+1}*B)*(B'*S{i+1}*A);
    S{i} = Q + A'*S{i+1}*A - A'*S{i+1}*B*inv(R+B'*S{i+1}*B)*B'*S{i+1}*A; 
end
X(:,1) = [2;1];
%X_dlqr(:,1) = [2;1];
% P_dlqr = dare(A,B,Q,R);
% K_dlqr = inv(R)*B'*P_dlqr;

for i = 2:N % (forwards in time)
    %u(:,i-1) =-inv(R)*B'*S{i-1}*X(:,i-1);
    u(:,i-1) = K(i,:)*X(:,i-1);
    X(:,i) = A * X(:,i-1) + B*u(:,i-1);
%     U_dlqr(:,i-1) = -K_dlqr*X_dlqr(:,i-1);
%     X_dlqr(:,i) =A * X_dlqr(:,i-1) + B*U_dlqr(:,i-1) ; 
end
%--------------------------------------------------------------------------
figure;
hold on
plot(t_dis(1:length(K(:,1))),K(:,1),'b','LineWidth',1.5);
plot(t_dis(1:length(K(:,2))),K(:,2),'r','LineWidth',1.5);
hold off
title('Feedback gain coefficients')
legend('K_{1}','K_{2}')
figure;
hold on
plot(t_dis,X(1,:),'b','LineWidth',1.5)
plot(t_dis,X(2,:),'r','LineWidth',1.5)
plot(t_dis(1:(length(u))),u,'k','LineWidth',1.5)
hold off
legend('x_1*','x_2*','u*')
% figure;
% hold on
% plot(t_dis(1:end-1),U_dlqr,'b','LineWidth',1.5)
% %plot(t_dis(1:end-1),U,'r','LineWidth',1.5)
% hold off
% ylabel('control')
% xlabel('time')
% legend('dlqr u*')
% figure;
% hold on
% plot(t_dis,X_dlqr,'LineWidth',1.5)
% %plot(t_dis,X,'LineWidth',1.5)
% hold off
% ylabel('states')
% xlabel('time')
% legend('Position','Velocity')



