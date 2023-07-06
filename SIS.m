%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final master's degree dissertation
%   MADOBIS 2022/2023
% Study of differential models applied to biological processes
% Laura Alvarez Valle
% laura.av.1999@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SIS (beta,gamma, N, I0, T)

% We consider the ODES:
%   |I' = beta*I*S - gamma*I
%   |S' = -beta*I*S + gamma*I
%   |I(0) = I0; S(0) = S0   con N = S0+I0

% Then, using that S = N-I, we obtain:
%
%   |I' = I.*(beta*N-beta*I-gamma)
%
% DERIVATES with respect parameters:
%
%   BETA:
%
%   |I_beta' = I*S + beta*Sbeta*I + beta*S*Ibeta- gamma*Ibeta
%   |S_beta' = -I*S - beta*Sbeta*I - beta*S*Ibeta+ gamma*Ibeta
%   |I(0) = I0; S(0) = S0   con N = S0+I0
%
%   GAMMA:
%
%   |I_gamma' = beta*Sgamma*I + beta*S*Igamma - I - gamma*Igamma
%   |S_gamma' = -beta*Sgamma*I - beta*S*Igamma + I + gamma*Igamma
%   |Igamma(0) = 0; Sgamma(0) = 0   
%

% Evaluation example:
% SIS(2.18,0.5, 100, 1, 4)

%% Data
S0 = N - I0;
y0 = [I0; S0]; %Initial values vector

%% Original ODE

% Equation
f = @(t,y) [beta*y(1)*y(2)/N - gamma*y(1); -beta*y(1)*y(2)/N + gamma*y(1) ];
%f = @(t,I) I.*(beta*N-beta*I-gamma);

% Resolution
[times, y] = ode45(f,[0,T],y0);
%[times,I] = ode45(f,[0,T],I0);
%S = N-I;

% Graphics
%figure(1)
subplot(1,3,1)
plot(times, y(:,1), 'b-', 'LineWidth', 1.8)
hold on
%plot(times, S, 'r-','LineWidth', 1.8)
plot(times, y(:,2), 'r--', 'LineWidth', 1.8)
plot(0,I0, 'k.', 'MarkerSize',18)
plot(0,S0, 'c.', 'MarkerSize',18)
legend("I","S","I_0", "S_0",  "Location","best")
title('SIS model', FontSize=14); xlabel('Time', FontSize=13); ylabel('Population', FontSize=13);
hold off
% 
% figure(2)
% plot(y(:,1), y(:,2))
% hold on
% plot(I0,S0, "r.", "MarkerSize",15)
% legend("Orbit", "Initial point", "Location","best")
% title("Orbit"); 
% xlabel("Infected"); ylabel("Susceptibles")
% hold off
% shg

%% SDO BETA
%
y0b = [0;0]; %Initial derivate point

%Function
fb = @(t,z) [ y(1)*y(2)/N + beta*z(2)*y(1)/N + beta*y(2)*z(1)/N - gamma*z(1);...
             -y(1)*y(2)/N - beta*z(2)*y(1)/N - beta*y(2)*z(1)/N+ gamma*z(1)];

[timesb,yb] = ode45(fb,[0,T],y0b);

%Graphics

%figure(3)
subplot(1,3,2)
plot(timesb, yb(:,1), 'b-' ,'LineWidth', 1.8)
hold on
plot(timesb, yb(:,2), 'r--' ,'LineWidth', 1.8)
plot(0,0, 'k.', 'MarkerSize',20)
plot(0,0, 'c.', 'MarkerSize',15)
legend("Derivate I_\beta","Derivate S_\beta","I_{\beta 0}", "S_{\beta0}",  "Location","best")
title('Derivates_\beta SIS', FontSize=14); xlabel('Time', FontSize=13); 
hold off
shg

%% SDO GAMMA

y0g = [0;0]; %Initial derivate point

%Funcion
fg = @(t,z) [ beta*z(2)*y(1)/N + beta*y(2)*z(1)/N - y(1) - gamma*z(1);...
             -beta*z(2)*y(1)/N - beta*y(2)*z(1)/N + y(1) + gamma*z(1)];

[timesg,yg] = ode45(fg,[0,T],y0g);

%Graphics

%figure(4)
subplot(1,3,3)
plot(timesg, yg(:,1), 'b-', 'LineWidth', 1.8)
hold on
plot(timesg, yg(:,2), 'r--', 'LineWidth', 1.8)
plot(0,0, 'k.', 'MarkerSize',20)
plot(0,0, 'c.', 'MarkerSize',15)
legend("Derivate I_\gamma","Derivate S_\gamma","I_{\gamma 0}", "S_{\gamma0}",  "Location","best")
title('Derivates_\gamma SIS', FontSize=14); xlabel('Time', FontSize=13);
hold off
shg

%% Interpolation

%Original
I = @(t) interp1(times, y(:,1), t) ;
S = @(t) interp1(times, y(:,2), t) ;

%Beta
Ib = @(t) interp1(timesb, yb(:,1), t) ;
Sb = @(t) interp1(timesb, yb(:,2), t) ;

%Gamma
Ig = @(t) interp1(timesg, yg(:,1), t) ;
Sg = @(t) interp1(timesg, yg(:,2), t) ;

%% SENSITIVITY MATRIX
t = linspace(0,T,5);
t = t(2:end);
lent  = length(t);

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

disp("Infected sensitivity matrix: ")
disp(SI)

%Tuning method
disp("Column 1 norm = ")
disp(norm(SI(:,1)))
disp("Column 2 norm = ")
disp(norm(SI(:,2)))

%Eigenvalue method
A1 = SI'*SI;
[V1,E1] = eig(A1);
disp("SI Eigenvectors Matrix:")
disp(V1)
disp("SI Eigenvalues Matrix:")
disp(E1)

%Correlation method
corSI = corr(A1);
disp("SI Correlation Matrix:")
disp(corSI)


% Susceptibles
SS = zeros(lent, 2);

for i = 1:lent
    for j = 1:2
        if j == 1
            SS(i,j) = Sb(t(i));
        else
            SS(i,j) = Sg(t(i));
        end

    end

end

% %Eigen value method
% A2 = SS'*SS;
% [V2,E2] = eig(A2);
% disp("SS Eigenvectors Matrix:")
% disp(V2)
% disp("SS Eigenvalues Matrix:")
% disp(E2)
% 
% %"Tuning" method
% disp("Column 1 norm = ")
% disp(norm(SS(:,1)))
% disp("Column 2 norm = ")
% disp(norm(SS(:,2)))

% Infected + Susceptibles

SB = zeros(2*lent, 2);

% Mixed sensitivity matrix
for i = 1:lent
    SB(2*i-1, :) = SI(i, :);  % Odd rows from SI
    SB(2*i, :) = SS(i, :);    % Even rows from SS
end

disp("Mixed sensitivity matrix: ")
disp(SB)

%"Tuning"method
disp("Column 1 norm = ")
disp(norm(SB(:,1)))
disp("Column 2 norm = ")
disp(norm(SB(:,2)))

%Eigenvalues method
A3 = SB'*SB;
[V3,E3] = eig(A3);
disp("SIS Eigenvectors Matrix:")
disp(V3)
disp("SIS Eigenvalues Matrix:")
disp(E3)

%Correlation method
corSB = corr(A3);
disp("SIS Correlation Matrix:")
disp(corSB)

end