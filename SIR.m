%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final master's degree dissertation
%   MADOBIS 2022/2023
% Study of differential models applied to biological processes
% Laura Alvarez Valle
% laura.av.1999@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SIR(beta, gamma, I0, R0, N, T)

% We consider the ODES:
%   |I' = beta*I*S - gamma*I
%   |R' = gamma*I
%   |S' = -beta*I*S
%   |I(0) = I0; S(0) = S0; R(0) = R0   con N = S0+I0+R0

% Derivates respect parameters:
%
%   BETA:
%
%   |I_beta' = I*S + beta*I_beta*S + beta*I*S_beta - gamma*I_beta
%   |R' = gamma*I_beta
%   |S' = -I*S - beta*I_beta*S - beta*I*S_beta
%   |I(0) = I0; S(0) = S0; R(0) = R0   con N = S0+I0+R0
%
%   GAMMA:
%
%   |I' = beta*I_gamma*S + beta*I*S_gamma - I - gamma*I_gamma
%   |R' = I + gamma*I_gamma
%   |S' = -beta*I_gamma*S - beta*I*S_gamma
%   |I(0) = I0; S(0) = S0; R(0) = R0   con N = S0+I0+R0

% Evaluation example
%SIR(0.2, 1, 0.4, 0.1, 1, 5)

%% Data
S0 = N - I0 -R0 ;
y0 = [I0; R0; S0];

%% Original ODE

% Equation
f = @(t,y) [(beta*y(3)-gamma)*y(1); ...
    gamma*y(1);...
    -beta*y(1)*y(3)];
    
% Resolution
[times,y] = ode45(f,[0,T],y0);

% Graphics
figure(1)
subplot(1,3,1)
plot(times, y(:,1), 'r-', 'LineWidth', 1.5) % Infected
hold on
plot(times, y(:,2), 'g-', 'LineWidth', 1.5) % Recovered
plot(times, y(:,3), 'b-', 'LineWidth', 1.5) % Susceptible
plot(0,I0, 'c.', 'MarkerSize',18)
plot(0,R0, 'k.', 'MarkerSize',18)
plot(0,S0, 'y.', 'MarkerSize',18)
legend("I","R","S","I0","R0", "S0",  "Location","best")
title('SIR Model'); xlabel('Time'); ylabel('Poblation');
hold off
% % Infected - Recovered
% figure(2)
% plot(y(:,1), y(:,2))
% hold on
% plot(I0,R0, "r.", "MarkerSize",18)
% legend("Orbit", "Initial value", "Location","best")
% title("Orbit"); 
% xlabel("Infected"); ylabel("Recovered")
% hold off
% % Infected - Susceptibles
% figure(3)
% plot(y(:,1), y(:,3))
% hold on
% plot(I0,S0, "r.", "MarkerSize",18)
% legend("Orbit", "Initial value", "Location","best")
% title("Orbit"); 
% xlabel("Infected"); ylabel("Susceptibles")
% hold off
% % Recovered - Susceptibles
% figure(4)
% plot(y(:,2), y(:,3))
% hold on
% plot(R0,S0, "r.", "MarkerSize",18)
% legend("Orbit", "Initial value", "Location","best")
% title("Orbit"); 
% xlabel("Recovered"); ylabel("Susceptibles")
% hold off
shg

%% SDO BETA
%
y0b = [0;0;0]; %Initial derivate point

%Function
fb = @(t,z) [ y(1)*y(3) + beta*z(1)*y(3) + beta*y(1)*z(3) - gamma*z(1);...
              gamma*z(1); ...
             -y(1)*y(3) - beta*z(1)*y(3) - beta*y(1)*z(3)];

[timesb,yb] = ode45(fb,[0,T],y0b);

%Graphics

%figure(5)
subplot(1,3,2)
plot(timesb, yb(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesb, yb(:,2),'g-' ,'LineWidth', 1.5)
plot(timesb, yb(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,T], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate I_\beta","Derivate R_\beta","Derivate S_\beta","I_{\beta 0}, R_{\beta0},S_{\beta0}",  "Location","best")
title('Derivates_\beta SIR'); xlabel('Time');
hold off
shg

%% SDO GAMMA
%
y0g = [0;0;0]; %Initial derivate point

%Function
fg = @(t,z) [ beta*z(1)*y(3) + beta*y(1)*z(3) - y(1)- gamma*z(1);...
              y(1) + gamma*z(1); ...
             - beta*z(1)*y(3) - beta*y(1)*z(3)];

[timesg,yg] = ode45(fg,[0,T],y0g);

%Graphics

%figure(6)
subplot(1,3,3)
plot(timesg, yg(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesg, yg(:,2),'g-' ,'LineWidth', 1.5)
plot(timesg, yg(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,T], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate I_\gamma","Derivate R_\gamma","Derivate S_\gamma","I_{\gamma 0}, R_{\gamma0},S_{\gamma0}",  "Location","best")
title('Derivates_\gamma SIR'); xlabel('Time');
hold off
shg

%% Interpolation

% Original
I = @(t) interp1(times, y(:,1), t) ;
R = @(t) interp1(times, y(:,2), t) ;
S = @(t) interp1(times, y(:,3), t) ;

% Beta

Ib = @(t) interp1(timesb, yb(:,1), t) ;
Rb = @(t) interp1(timesb, yb(:,2), t) ;
Sb = @(t) interp1(timesb, yb(:,3), t) ;

% Gamma

Ig = @(t) interp1(timesg, yg(:,1), t) ;
Rg = @(t) interp1(timesg, yg(:,2), t) ;
Sg = @(t) interp1(timesg, yg(:,3), t) ;


%% SENSITIVITY MATRIX
t = linspace(0,T,5);
t = t(2:end);
lent  = length(t);
np = 2; %number of parameters

% Infected
SI = zeros(lent, 2);

for i = 1:lent
    for j = 1:2
        if j == 1
            SI(i,j) = Ib(t(i));
        else
            SI(i,j) = Ig(t(i));
        end

    end

end

disp("Infected sensivity matrix: ")
disp(SI)

% Tuning method
for i = 1:np
    a = norm(SI(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalue method
A1 = SI'*SI;
[V1,E1] = eig(A1);
disp("SI Eigenvectors Matrix:")
disp(V1) %eigen vectors matrix
disp("SI Eigenvalues Matrix:")
disp(E1) %eigen values matrix

%Correlation method
corSI = corr(A1);
disp("SI Correlation Matrix:")
disp(corSI)

% Recovered
SR = zeros(lent, 2);

for i = 1:lent
    for j = 1:np
        if j == 1
            SR(i,j) = Rb(t(i));
        else
            SR(i,j) = Rg(t(i));
        end

    end

end

disp("Recovered sensivity matrix: ")
disp(SR)

% Tuning method
for i = 1:np
    a = norm(SR(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalue method
A2 = SR'*SR;
[V2,E2] = eig(A2);
disp("SR Eigenvectors Matrix:")
disp(V2) %eigen vectors matrix
disp("SR Eigenvalues Matrix:")
disp(E2) %eigen values matrix

%Correlation method
corSR = corr(A2);
disp("SR Correlation Matrix:")
disp(corSR)

% Susceptibles
SS = zeros(lent, 2);

for i = 1:lent
    for j = 1:np
        if j == 1
            SS(i,j) = Sb(t(i));
        else
            SS(i,j) = Sg(t(i));
        end

    end

end

%disp("Susceptible sensivity matrix: ")
%disp(SS)

% % Tuning method
% for i = 1:np
%     norm(SS(:,i))
% end
%
% % Eigenvalue method
% A3 = SS'*SS;
% [V3,E3] = eig(A3);
%disp(V3) %eigen vectors matrix
%disp(E3) %eigen values matrix

% Infected + Recovered

SIR = zeros(2*lent, np);

%Matrix:
for i = 1:lent
    SIR(2*i-1, :) = SI(i, :);
    SIR(2*i, :) = SR(i, :);
end

disp("Infected and recovered sensitivity matrix:")
disp(SIR)

% Tuning method
for i = 1:np
    a = norm(SIR(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalue method
A4 = SIR'*SIR;
[V4, E4] = eig(A4);
disp("SIR Eigenvectors Matrix:")
disp(V4) %eigen vector matrix
disp("SIR Eigenvalues Matrix:")
disp(E4) %eigen values matrix

%Correlation method
corSIR = corr(A4);
disp("SIR Correlation Matrix:")
disp(corSIR)


% Infected + Recovered + Susceptibles

SIRS = zeros(3*lent, np);

%Matrix:
for i = 1:lent
    SIRS(3*i-2, :) = SI(i, :);
    SIRS(3*i-1, :) = SR(i, :);
    SIRS(3*i, :) = SS(i,:);
end

disp("Complete sensitivity matrix:")
disp(SIRS)

% Tuning method
for i = 1:np
    a = norm(SIRS(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalue method
A5 = SIRS'*SIRS;
[V5, E5] = eig(A5);
disp("SIRS Eigenvectors Matrix:")
disp(V5) %eigenvector matrix
disp("SIRS Eigenvalues Matrix:")
disp(E5) %eigenvalues matrix

%Correlation method
corSIRS = corr(A5);
disp("SIRS Correlation Matrix:")
disp(corSIRS)



end