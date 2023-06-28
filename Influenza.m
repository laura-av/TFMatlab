%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final master's degree dissertation
%   MADOBIS 2022/2023
% Study of differential models applied to biological processes
% Laura Alvarez Valle
% laura.av.1999@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Influenza

% We consider the ODES:
%   |T' = -beta*T*V
%   |I' = beta*T*V - delta*I
%   |V' = p*I - c*V

% DERIVATES:

%%% 
% beta:
%
%   |Tb' = -T*V -beta*T_b*V -beta*T*V_b
%   |Ib' = T*V + beta*T_b*V + beta*T*V_b - delta*I_b
%   |Vb' = p*I_b - c*V_b

%%% 
% delta:
%
%   |Td' = -beta*T_d*V -beta*T*V_d
%   |Id' = beta*T_d*V + beta*T*V_d -I - delta*I_d
%   |Vd' = p*I_d - c*V_d

%%%
% p:
%
%   |Tp' = -beta*T_p*V -beta*T*V_p
%   |Ip' = beta*T_p*V + beta*T*V_p - delta*I_p
%   |Vp' = I + p*I_p - c*V_p

%%%
% c:
%
%   |Tc' = -beta*T_c*V -beta*T*V_c
%   |Ic' = beta*T_c*V + beta*T*V_c - delta*I_c
%   |Vc' = p*I_c - V - c*V_c

%% Data
T0 = 0.8;
I0 = 2;
V0 = 1.50;
y0 = [T0; I0; V0];
TT = 5;
beta = 0.01;
delta = 0.5;
p = 0.2;
c = 2;

%% Original ODES

% Equation
f = @(t,y) [-beta*y(1)*y(3);  beta*y(1)*y(3) - delta*y(2);  p*y(2)-c*y(3)];
   
% Resolution
[times,y] = ode23(f,[0,TT],y0);

% Graphics
figure(1)
subplot(2,3,1)
plot(times, y(:,1), 'r-', 'LineWidth', 1.5) % Uninfected target cells
hold on
plot(times, y(:,2), 'g-', 'LineWidth', 1.5) % Productively infected cells
plot(times, y(:,3), 'b-', 'LineWidth', 1.5) % Infectious viral titer
plot(0,T0, 'c.', 'MarkerSize',18)
plot(0,I0, 'k.', 'MarkerSize',18)
plot(0,V0, 'y.', 'MarkerSize',18)
legend("T","I","V","T0","I0", "V0",  "Location","best")
title('Influenza Model'); xlabel('Time'); ylabel('Poblation');
hold off

%% SDO BETA
%
y0b = [0;0;0]; %Initial derivate point

% Equation
fb = @(t,z) [-y(1)*y(3)-beta*z(1)*y(3)-beta*y(1)*z(3); ...
    y(1)*y(3)+beta*z(1)*y(3)-beta*y(1)*z(3)-delta*z(2);...
    p*z(2)-c*z(3)];
    
% Resolution
[timesb,yb] = ode45(fb,[0,TT],y0b);

% %Graphics
% 
% figure(2)
subplot(2,3,2)
plot(timesb, yb(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesb, yb(:,2),'g-' ,'LineWidth', 1.5)
plot(timesb, yb(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,TT], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate T_\beta","Derivate I_\beta","Derivate V_\beta","T_{\beta 0}, I_{\beta0},V_{\beta0}",  "Location","best")
title('Derivates_\beta Influenza Model'); xlabel('Time');
hold off
shg

%% SDO DELTA
%
y0d = [0;0;0]; %Initial derivate point

%Function
fd = @(t,z) [ -beta*z(1)*y(3)-beta*y(1)*z(3);...
             beta*z(1)*y(3)+beta*y(1)*z(3)-y(2)-delta*z(2); ...
             p*z(2) - c*z(3)];

[timesd,yd] = ode45(fd,[0,TT],y0d);

% %Graphics
% 
% figure(3)
subplot(2,3,3)
plot(timesd, yd(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesd, yd(:,2),'g-' ,'LineWidth', 1.5)
plot(timesd, yd(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,TT], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate T_\delta","Derivate I_\delta","Derivate V_\delta","T_{\delta 0}, I_{\delta0},V_{\delta0}",  "Location","best")
title('Derivates_\delta Influenza Model'); xlabel('Time');
hold off
shg

%% SDO p
%

%Function
fp = @(t,z) [ -beta*z(1)*y(3)-beta*y(1)*z(3);...
             beta*z(1)*y(3)+beta*y(1)*z(3)-y(2)-delta*z(2); ...
             y(2) + p*z(2) - c*z(3)];

[timesp,yp] = ode45(fp,[0,TT],y0d);
% 
% %Graphics
% 
% figure(4)
subplot(2,3,4)
plot(timesp, yp(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesp, yp(:,2),'g-' ,'LineWidth', 1.5)
plot(timesp, yp(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,TT], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate T_p","Derivate I_p","Derivate V_p","T_{p 0}, I_{p0},V_{p0}",  "Location","best")
title('Derivates_p Influenza Model'); xlabel('Time');
hold off
shg

%% SDO c
%

%Function
fc = @(t,z) [ -beta*z(1)*y(3)-beta*y(1)*z(3);...
             beta*z(1)*y(3)+beta*y(1)*z(3)-delta*z(2); ...
             p*z(2)-y(3)-c*z(3)];

[timesc,yc] = ode45(fc,[0,TT],y0d);

% %Graphics
% 
% figure(5)

subplot(2,3,5)
plot(timesc, yc(:,1),'r-' ,'LineWidth', 1.5)
hold on
plot(timesc, yc(:,2),'g-' ,'LineWidth', 1.5)
plot(timesc, yc(:,3),'b-' ,'LineWidth', 1.5)
plot(0,0, 'y.', 'MarkerSize',18)
plot([0,TT], [0,0], '--k', 'LineWidth', 0.5)
legend("Derivate T_c","Derivate I_c","Derivate V_c","T_{c 0}, I_{c0},V_{c0}",  "Location","best")
title('Derivates_c Influenza Model'); xlabel('Time');
hold off
shg

%% Interpolation

% Original
T = @(t) interp1(times, y(:,1), t) ;
I = @(t) interp1(times, y(:,2), t) ;
V = @(t) interp1(times, y(:,3), t) ;

% beta
Tb = @(t) interp1(timesb, yb(:,1), t) ;
Ib = @(t) interp1(timesb, yb(:,2), t) ;
Vb = @(t) interp1(timesb, yb(:,3), t) ;

% delta
Td = @(t) interp1(timesd, yd(:,1), t) ;
Id = @(t) interp1(timesd, yd(:,2), t) ;
Vd = @(t) interp1(timesd, yd(:,3), t) ;

% p
Tp = @(t) interp1(timesp, yp(:,1), t) ;
Ip = @(t) interp1(timesp, yp(:,2), t) ;
Vp = @(t) interp1(timesp, yp(:,3), t) ;

% c
Tc = @(t) interp1(timesc, yc(:,1), t) ;
Ic = @(t) interp1(timesc, yc(:,2), t) ;
Vc = @(t) interp1(timesc, yc(:,3), t) ;


%% Sensitivity matrix

t = linspace(0,TT,5);
t = t(2:end); %time-points vector
lent  = length(t);
np = 4; % number of parameters


% T: Uninfected target cells

ST = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            ST(i,j) = Tb(t(i));
        elseif j == 2
            ST(i,j) = Td(t(i));
        elseif j == 3
            ST(i,j) = Tp(t(i));
        else
            ST(i,j) = Tc(t(i));
        end

    end

end

disp("Uninfected cells sensivity matrix: ")
disp(ST)

% Eigenvalues
A1 = ST'*ST;
[V1,E1] = eig(A1);
disp("ST Eigenvectors Matrix:")
disp(V1) %eigen vector matrix
disp("ST Eigenvalues Matrix:")
disp(E1) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(ST(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corST = corr(A1);
disp("ST Correlation Matrix:")
disp(corST)

% I: Infected target cells

SI = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            SI(i,j) = Ib(t(i));
        elseif j == 2
            SI(i,j) = Id(t(i));
        elseif j == 3
            SI(i,j) = Ip(t(i));
        else
            SI(i,j) = Ic(t(i));
        end

    end

end

disp("Infected cells sensivity matrix: ")
disp(SI)

% Eigenvalues
A2 = SI'*SI;
[V2,E2] = eig(A2);
disp("SI Eigenvectors Matrix:")
disp(V2) %eigen vector matrix
disp("SI Eigenvalues Matrix:")
disp(E2) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(SI(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSI = corr(A2);
disp("SI Correlation Matrix:")
disp(corSI)

% V: Virus titer

SV = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            SV(i,j) = Vb(t(i));
        elseif j == 2
            SV(i,j) = Vd(t(i));
        elseif j == 3
            SV(i,j) = Vp(t(i));
        else
            SV(i,j) = Vc(t(i));
        end

    end

end

disp("Virus sensivity matrix: ")
disp(SV)

% Eigenvalues
A3 = SV'*SV;
[V3,E3] = eig(A3);
disp("SV Eigenvectors Matrix:")
disp(V3) %eigen vector matrix
disp("SV Eigenvalues Matrix:")
disp(E3) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(SV(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSV = corr(A3);
disp("SV Correlation Matrix:")
disp(corSV)

% Uninfected + Infected
STI = zeros(2*lent, np);

%Matrix:
for i = 1:lent
    STI(2*i-1, :) = ST(i, :);
    STI(2*i, :) = SI(i, :);
end

disp("Uninfected + Infected sensivity matrix: ")
disp(STI)

% Eigenvalues
A4 = STI'*STI;
[V4,E4] = eig(A4);
disp("STI Eigenvectors Matrix:")
disp(V4) %eigen vector matrix
disp("STI Eigenvalues Matrix:")
disp(E4) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(STI(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSTI = corr(A4);
disp("STI Correlation Matrix:")
disp(corSTI)

% Uninfected + Virus
STV = zeros(2*lent, np);

%Matrix:
for i = 1:lent
    STV(2*i-1, :) = ST(i, :);
    STV(2*i, :) = SV(i, :);
end
disp("Uninfected + Virus sensivity matrix: ")
disp(STV)

% Eigenvalues
A5 = STV'*STV;
[V5,E5] = eig(A5);
disp("STV Eigenvectors Matrix:")
disp(V5) %eigen vector matrix
disp("STV Eigenvalues Matrix:")
disp(E5) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(STV(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSTV = corr(A5);
disp("STI Correlation Matrix:")
disp(corSTV)

% Infected + Virus
SIV = zeros(2*lent, np);

%Matrix:
for i = 1:lent
    SIV(2*i-1, :) = SI(i, :);
    SIV(2*i, :) = SV(i, :);
end

disp("Infected + Virus sensivity matrix: ")
disp(SIV)

% Eigenvalues
A6 = SIV'*SIV;
[V6,E6] = eig(A6);
disp("SIV Eigenvectors Matrix:")
disp(V6) %eigen vector matrix
disp("SIV Eigenvalues Matrix:")
disp(E6) %eigen values matrix

% Tuning method
for i = 1:np
    a = norm(SIV(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSIV = corr(A6);
disp("SIV Correlation Matrix:")
disp(corSIV)

% Mixed Matrix: Uninfected + Infected + Virus
STIV = zeros(3*lent, np);

%Matrix:
for i = 1:lent
    STIV(3*i-2, :) = ST(i, :);
    STIV(3*i-1, :) = SI(i, :);
    STIV(3*i, :) = SV(i, :);
end

disp("Mixed sensivity matrix: ")
disp(STIV)

% Eigenvalues
A7 = STIV'*STIV;
[V7,E7] = eig(A7);
disp("STIV Eigenvectors Matrix:")
disp(V7) %eigenvector matrix
disp("STIV Eigenvalues Matrix:")
disp(E7) %eigenvalues matrix

% Tuning method
for i = 1:np
    a = norm(STIV(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

%Correlation method
corSTIV = corr(A7);
disp("STIV Correlation Matrix:")
disp(corSTIV)

end