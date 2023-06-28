%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final master's degree dissertation
%   MADOBIS 2022/2023
% Study of differential models applied to biological processes
% Laura Alvarez Valle
% laura.av.1999@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HIV

% We consider the ODES:
%   |T'   = (p - km*Tm - kw*Tw -kR*Tmw)*T
%   |Tm'  = (pm + km*T -qm*Tw)*Tm + 0.25*kR*Tmw*T
%   |Tw'  = (pw + kw*T - qw*Tm)*Tw + 0.25*kR*Tmw*T
%   |Tmw' = (pmw + 0.5*kR*T)*Tmw + (qm + qw)*Tm*Tw
%
% where:
%
%  (lambda, lambda_m, lambda_w, lambda_mw) = proliferation rates of T, Tm,
%  Tw and Tmw
%
%  (delta, delta_m, delta_w, delta_mw) = death rates
%
%   (p, pm, pw, pmw) = growth rates of T, Tm, Tw and Tmw
%
%   (km, kw, kR) =infection rates mutant, wild-type and recombinant
%
%   (qm, qw) = infection rate recombinant Tw -> Tmw and Tm -> Tmw
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DERIVATES:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RHO:
%
% p:
%
%   |T_p'   = (1 - km*Tm_p - kw*Tw_p -kR*Tmw_p)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_p
%
%   |Tm_p'  = (km*T_p -qm*Tw_p)*Tm  
%               + (pm + km*T -qm*Tw)*Tm_p
%               + 0.25*kR*Tmw_p*T + 0.25*kR*Tmw*T_p
%
%   |Tw_p'  = (kw*T_p - qw*Tm_p)*Tw 
%               + (pw + kw*T - qw*Tm)*Tw_p
%               + 0.25*kR*Tmw_p*T + 0.25*kR*Tmw*T_p
%
%   |Tmw_p' = (0.5*kR*T_p)*Tmw + (pmw + 0.5*kR*T)*Tmw_p
%               + (qm + qw)*Tm_p*Tw + (qm + qw)*Tm*Tw_p
%
%%%
%
% pm:
%
%   |T_pm'   = -(km*Tm_pm + kw*Tw_pm +kR*Tmw_pm)*T
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_pm
%
%   |Tm_pm'  = (1 + km*T_pm -qm*Tw_pm)*Tm 
%               +(pm + km*T -qm*Tw)*Tm_pm  
%               + 0.25*kR*Tmw_pm*T + 0.25*kR*Tmw*T_pm
%
%   |Tw_pm'  = (kw*T_pm - qw*Tm_pm)*Tw 
%               +(pw + kw*T - qw*Tm)*Tw_pm 
%               + 0.25*kR*Tmw_pm*T + 0.25*kR*Tmw*T_pm
%
%   |Tmw_pm' = (0.5*kR*T_pm)*Tmw 
%               + (pmw + 0.5*kR*T)*Tmw_pm
%               + (qm + qw)*Tm_pm*Tw + (qm + qw)*Tm*Tw_pm
%
%%%
%
% pw:
%
%   |T_pw'   = -(km*Tm_pw + kw*Tw_pw + kR*Tmw_pw)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_pw
%
%   |Tm_pw'  = (km*T_pw -qm*Tw_pw)*Tm 
%               + (pm + km*T -qm*Tw)*Tm_pw 
%               + 0.25*kR*Tmw_pw*T + 0.25*kR*Tmw*T_pw
%
%   |Tw_pw'  = (1 + kw*T_pw - qw*Tm_pw)*Tw 
%               + (pw + kw*T - qw*Tm)*Tw_pw
%               + 0.25*kR*Tmw_pw*T + 0.25*kR*Tmw*T_pw
%
%   |Tmw' = (0.5*kR*T_pw)*Tmw + (pmw + 0.5*kR*T)*Tmw_pw
%               + (qm + qw)*Tm_pw*Tw + (qm + qw)*Tm*Tw_pw
%
%%%
%
% pmw:
%
%   |T_pmw'   = -(km*Tm_pmw + kw*Tw_pmw +kR*Tmw_pmw)*T 
%                + (p - km*Tm - kw*Tw -kR*Tmw)*T_pmw
%
%   |Tm_pmw'  = (km*T_pmw -qm*Tw_pmw)*Tm 
%               + (pm + km*T -qm*Tw)*Tm_pmw 
%               + 0.25*kR*Tmw_pmw*T + 0.25*kR*Tmw*T_pmw
%
%   |Tw_pmw'  = (kw*T_pmw - qw*Tm_pmw)*Tw 
%               + (pw + kw*T - qw*Tm)*Tw_pmw 
%               + 0.25*kR*Tmw_pmw*T + 0.25*kR*Tmw*T_pmw
%
%   |Tmw' = (1 + 0.5*kR*T_pmw)*Tmw + (pmw + 0.5*kR*T)*Tmw_pmw 
%               + (qm + qw)*Tm_pmw*Tw + (qm + qw)*Tm*Tw_pmw
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%  k:
%
%%%%
%
% km:
%
%   |T_km'   = -(Tm + km*Tm_km + kw*Tw_km +kR*Tmw_km)*T
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_km
%
%   |Tm_km'  = (T + km*T_km -qm*Tw_km)*Tm 
%               + (pm + km*T -qm*Tw)*Tm_km 
%               + 0.25*kR*Tmw_km*T + 0.25*kR*Tmw*T_km
%
%   |Tw_km'  = (kw*T_km - qw*Tm_km)*Tw 
%               + (pw + kw*T - qw*Tm)*Tw_km 
%               + 0.25*kR*Tmw_km*T + 0.25*kR*Tmw*T_km
%
%   |Tmw_km' = (0.5*kR*T_km)*Tmw + (pmw + 0.5*kR*T)*Tmw_km
%               + (qm + qw)*Tm_km*Tw + (qm + qw)*Tm*Tw_km
%
%%%%
%
%  kw:
%
%   |T_kw'   = -(km*Tm_kw + Tw + kw*Tw_kw +kR*Tmw_kw)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_kw
%
%   |Tm_kw'  = (km*T_kw -qm*Tw_kw)*Tm 
%               + (pm + km*T -qm*Tw)*Tm_kw 
%               +  0.25*kR*Tmw_kw*T + 0.25*kR*Tmw*T_kw
%
%   |Tw_km'  = (T + kw*T_kw - qw*Tm_kw)*Tw 
%               + (pw + kw*T - qw*Tm)*Tw_kw 
%               + 0.25*kR*Tmw*T + 0.25*kR*Tmw*T
%
%   |Tmw_km' = (0.5*kR*T_km)*Tmw + (pmw + 0.5*kR*T)*Tmw_km 
%               + (qm + qw)*Tm_km*Tw + (qm + qw)*Tm*Tw_km
%
%%%%
%
%  kR:
%
%   |T_kR'   = -(km*Tm_kR + kw*Tw_kR + Tmw + kR*Tmw_kR)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_kR
%
%   |Tm_kR'  = (km*T_kR -qm*Tw_kR)*Tm + (pm + km*T -qm*Tw)*Tm_kR 
%               + 0.25*Tmw*T + 0.25*kR*Tmw_kR*T + 0.25*kR*Tmw*T_kR
%
%   |Tw_kR'  = (kw*T_kR - qw*Tm_kR)*Tw + (pw + kw*T - qw*Tm)*Tw_kR 
%               + 0.25*Tmw*T + 0.25*kR*Tmw_kR*T + 0.25*kR*Tmw*T_kR
%
%   |Tmw_kR' = (0.5*T + 0.5*kR*T_kR)*Tmw + (pmw + 0.5*kR*T)*Tmw_kR 
%               + (qm + qw)*Tm_kR*Tw + (qm + qw)*Tm*Tw_kR
%
%%%%%%%%%%%%%%%%%%%%%%%%%
%  q:
%
%%%%
%
%  qm:
%
%   |T_qm'   = -(km*Tm_qm + kw*Tw_qm -kR*Tmw_qm)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_qm
%
%   |Tm_qm'  = (km*T_qm - Tw -qm*Tw_qm)*Tm + (pm + km*T -qm*Tw)*Tm_qm  
%               + 0.25*kR*Tmw_qm*T + 0.25*kR*Tmw*T_qm
%
%   |Tw_qm'  = (kw*T - qw*Tm_qm)*Tw + (pw + kw*T - qw*Tm)*Tw_qm 
%               + 0.25*kR*Tmw_qm*T + 0.25*kR*Tmw*T_qm
%
%   |Tmw_qm' = (0.5*kR*T_qm)*Tmw + (pmw + 0.5*kR*T)*Tmw_qm
%               + Tm*Tw + (qm + qw)*Tm_qm*Tw + (qm + qw)*Tm*Tw_qm
%
%%%
%
%  qw:
%
%   |T_qw'   = -(km*Tm_qw + kw*Tw_qw +kR*Tmw_qw)*T 
%               + (p - km*Tm - kw*Tw -kR*Tmw)*T_qw
%
%   |Tm_qw'  = (km*T_qw -qm*Tw_qw)*Tm + (pm + km*T -qm*Tw)*Tm_qw 
%               + 0.25*kR*Tmw_qw*T + 0.25*kR*Tmw*T_qw
%
%   |Tw_qw'  = (kw*T_qw - Tm - qw*Tm_qw)*Tw + (pw + kw*T - qw*Tm)*Tw_qw 
%               + 0.25*kR*Tmw_qw*T + 0.25*kR*Tmw*T_qw
%
%   |Tmw_qw' = (0.5*kR*T_qw)*Tmw + (pmw + 0.5*kR*T)*Tmw_qw 
%               + Tm*Tw + (qm + qw)*Tm_qw*Tw + (qm + qw)*Tm*Tw_qw
%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data
p = .002;    %Growth rate uninfected
pm = .0017;  %Growth rate mutant infected
pw = 0.034;  %Growth rate wild-type infected
pmw = 0.4;   %Growth rate recombinant infected
%Infection rates
km = 0.003; %Infection rate mutant virus
kw = 0.004; %Infection rate wild-type virus
kR = 0.026; %Infection rate recombinant virus

kr4 = 0.25*kR;
kr2 = 0.5*kR;
%Dual infection rates:
qm = 0.003; 
qw = 0.009;
%

T0 = .50;  %Initial uninfected
Tm0 = .10; %Initial mutant infected
Tw0 = .10; %Initial wild-type infected
Tmw0 = .30; %Initial recombinant infected
TT = 6; %Time

y0 = [T0; Tm0; Tw0; Tmw0]; %Initial values vector

%% Original ODE

% Equation
f = @(t,y) [(p - km*y(2) - kw*y(3) -kR*y(4))*y(1); ...
            (pm + km*y(1) -qm*y(3))*y(2) + kr4*y(4)*y(1); ...
            (pw + kw*y(1) - qw*y(2))*y(3) + kr4*y(4)*y(1); ...
            (pmw + kr2*y(1))*y(4) + (qm + qw)*y(2)*y(3)];


% Resolution
[times,y] = ode45(f,[0,TT],y0);

% Graphics
figure(1)
subplot(2,5,1)
plot(times, y,'-', 'LineWidth', 1.5) % Infected
hold on
plot(0,T0, 'c.', 'MarkerSize',18)
plot(0,Tm0, 'k.', 'MarkerSize',20)
plot(0,Tw0, 'y.', 'MarkerSize',18)
plot(0,Tmw0, 'm.', 'MarkerSize',18)
legend("T","Tm","Tw","Tmw", "T0","Tm0", "Tw0", "Tmw0",  "Location","best")
title('HIV Model '); xlabel('Time'); ylabel('Poblation');
hold off

%% SDO p
y00 = [0;0;0;0]; %Initial derivate point

% Equation
fp = @(t,z) [(1 - km*z(2) - kw*z(3) -kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timesp,yp] = ode45(fp,[0,TT],y00);

% % Graphics
% figure(2)
subplot(2,5,2)
plot(timesp, yp,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_p","Tm'_p","Tw'_p","Tmw'_p", "Location","best")
title('Derivate p HIV Model '); xlabel('Time');
hold off

%% SDO pm

% Equation
fpm = @(t,z) [-(km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (1 + km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timespm,ypm] = ode45(fpm,[0,TT],y00);

% % Graphics
% figure(3)
subplot(2,5,3)
plot(timespm, ypm,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{pm}","Tm'_{pm}","Tw'_{pm}","Tmw'_{pm}", "Location","best")
title('Derivate pm HIV Model '); xlabel('Time');
hold off

%% SDO pw

% Equation
fpw = @(t,z) [-(km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (1 + kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timespw,ypw] = ode45(fpw,[0,TT],y00);

% % Graphics
% figure(4)
subplot(2,5,4)
plot(timespw, ypw,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{pw}","Tm'_{pw}","Tw'_{pw}","Tmw'_{pw}",  "Location","best")
title('Derivate pw HIV Model '); xlabel('Time');
hold off

%% SDO pmw

% Equation
fpmw = @(t,z) [-(km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
                (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
                (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
                (1 + kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timespmw,ypmw] = ode45(fpmw,[0,TT],y00);

% Graphics
% figure(5)
subplot(2,5,5)
plot(timespmw, ypmw,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{pmw}","Tm'_{pmw}","Tw'_{pmw}","Tmw'_{pmw}", "Location","best")
title('Derivate pmw HIV Model '); xlabel('Time');
hold off

%% SDO km

% Equation
fkm = @(t,z) [-(y(2) + km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (y(1) + km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + + (qm + qw)*y(2)*z(3)];


% Resolution
[timeskm,ykm] = ode45(fkm,[0,TT],y00);

% % Graphics
% figure(6)
subplot(2,5,6)
plot(timeskm, ykm,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{km}","Tm'_{km}","Tw'_{km}","Tmw'_{km}", "Location","best")
title('Derivate km HIV Model '); xlabel('Time');
hold off

%% SDO kw

% Equation
fkw = @(t,z) [-(km*z(2) + y(3) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (y(1) + kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timeskw,ykw] = ode45(fkw,[0,TT],y00);

% % Graphics
% figure(7)
subplot(2,5,7)
plot(timeskw, ykw,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{kw}","Tm'_{kw}","Tw'_{kw}","Tmw'_{kw}",  "Location","best")
%"T0_{kw}","Tm0_{kw}", "Tw0_{kw}", "Tmw0_{kw}"
title('Derivate kw HIV Model '); xlabel('Time');
hold off

%% SDO kR

% Equation
fkR = @(t,z) [-(km*z(2) + kw*z(3) + y(4) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + 0.25*y(4)*y(1) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + 0.25*y(4)*y(1) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (0.5*y(1) + kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timeskR,ykR] = ode45(fkR,[0,TT],y00);

% % Graphics
% figure(8)
subplot(2,5,8)
plot(timeskR, ykR,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{kR}","Tm'_{kR}","Tw'_{kR}","Tmw'_{kR}", "Location","best")
title('Derivate kR HIV Model '); xlabel('Time');
hold off

%% SDO qm

% Equation
fqm = @(t,z) [-(km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) - y(3) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + y(2)*y(3) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timesqm,yqm] = ode45(fqm,[0,TT],y00);

% % Graphics
% figure(9)
subplot(2,5,9)
plot(timesqm, yqm,'-', 'LineWidth', 1.8)
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{qm}","Tm'_{qm}","Tw'_{qm}","Tmw'_{qm}", "Location","best")
title('Derivate qm HIV Model '); xlabel('Time');
hold off

%% SDO qw

% Equation
fqw = @(t,z) [-(km*z(2) + kw*z(3) +kR*z(4))*y(1) + (p - km*y(2) - kw*y(3) -kR*y(4))*z(1); ...
            (km*z(1) -qm*z(3))*y(2) + (pm + km*y(1) -qm*y(3))*z(2) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kw*z(1) - y(2) - qw*z(2))*y(3)+ (pw + kw*y(1) - qw*y(2))*z(3) + kr4*z(4)*y(1) + kr4*y(4)*z(1); ...
            (kr2*z(1))*y(4) + (pmw + kr2*y(1))*z(4) + y(2)*y(3) + (qm + qw)*z(2)*y(3) + (qm + qw)*y(2)*z(3)];


% Resolution
[timesqw,yqw] = ode45(fqw,[0,TT],y00);

% Graphics
%figure(10)
subplot(2,5,10)
plot(timesqw, yqw,'-', 'LineWidth', 1.5) 
hold on
plot(0,y00(1), 'k.', 'MarkerSize',18)
plot([0,TT], [0,0], 'k--', 'LineWidth', 0.5)
legend("T'_{qw}","Tm'_{qw}","Tw'_{qw}","Tmw'_{qw}",  "Location","best")
title('Derivate qw HIV Model '); xlabel('Time');
hold off

%% Interpolation

% Original
T = @(t) interp1(times, y(:,1), t) ;
Tm = @(t) interp1(times, y(:,2), t) ;
Tw = @(t) interp1(times, y(:,3), t) ;
Tmw = @(t) interp1(times, y(:,4), t) ;

% p
T_p = @(t) interp1(timesp, yp(:,1), t) ;
Tm_p = @(t) interp1(timesp, yp(:,2), t) ;
Tw_p = @(t) interp1(timesp, yp(:,3), t) ;
Tmw_p = @(t) interp1(timesp, yp(:,4), t) ;

% pm
T_pm = @(t) interp1(timespm, ypm(:,1), t) ;
Tm_pm = @(t) interp1(timespm, ypm(:,2), t) ;
Tw_pm = @(t) interp1(timespm, ypm(:,3), t) ;
Tmw_pm = @(t) interp1(timespm, ypm(:,4), t) ;

%pw
T_pw = @(t) interp1(timespw, ypw(:,1), t) ;
Tm_pw = @(t) interp1(timespw, ypw(:,2), t) ;
Tw_pw = @(t) interp1(timespw, ypw(:,3), t) ;
Tmw_pw = @(t) interp1(timespw, ypw(:,4), t) ;

%pmw
T_pmw = @(t) interp1(timespmw, ypmw(:,1), t) ;
Tm_pmw = @(t) interp1(timespmw, ypmw(:,2), t) ;
Tw_pmw = @(t) interp1(timespmw, ypmw(:,3), t) ;
Tmw_pmw = @(t) interp1(timespmw, ypmw(:,4), t) ;

%km
T_km = @(t) interp1(timeskm, ykm(:,1), t) ;
Tm_km = @(t) interp1(timeskm, ykm(:,2), t) ;
Tw_km = @(t) interp1(timeskm, ykm(:,3), t) ;
Tmw_km = @(t) interp1(timeskm, ykm(:,4), t) ;

%kw
T_kw = @(t) interp1(timeskw, ykw(:,1), t) ;
Tm_kw = @(t) interp1(timeskw, ykw(:,2), t) ;
Tw_kw = @(t) interp1(timeskw, ykw(:,3), t) ;
Tmw_kw = @(t) interp1(timeskw, ykw(:,4), t) ;

%kR
T_kR = @(t) interp1(timeskR, ykR(:,1), t) ;
Tm_kR = @(t) interp1(timeskR, ykR(:,2), t) ;
Tw_kR = @(t) interp1(timeskR, ykR(:,3), t) ;
Tmw_kR = @(t) interp1(timeskR, ykR(:,4), t) ;

%qm
T_qm = @(t) interp1(timesqm, yqm(:,1), t) ;
Tm_qm = @(t) interp1(timesqm, yqm(:,2), t) ;
Tw_qm = @(t) interp1(timesqm, yqm(:,3), t) ;
Tmw_qm = @(t) interp1(timesqm, yqm(:,4), t) ;

%qw
T_qw = @(t) interp1(timesqw, yqw(:,1), t) ;
Tm_qw = @(t) interp1(timesqw, yqw(:,2), t) ;
Tw_qw = @(t) interp1(timesqw, yqw(:,3), t) ;
Tmw_qw = @(t) interp1(timesqw, yqw(:,4), t) ;

%% SENSITIVITY MATRIX

t = linspace(0,TT,5); 
t = t(2:end); %time-points vector
lent  = length(t);
np = 9; % number of parameters

% T: Uninfected Cells

ST = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            ST(i,j) = T_p(t(i));
        elseif j == 2
            ST(i,j) = T_pm(t(i));
        elseif j == 3
            ST(i,j) = T_pw(t(i));
        elseif j == 4
            ST(i,j) = T_pmw(t(i));
        elseif j == 5
            ST(i,j) = T_km(t(i));
        elseif j == 6
            ST(i,j) = T_kw(t(i));
        elseif j == 7
            ST(i,j) = T_kR(t(i));
        elseif j == 8
            ST(i,j) = T_qm(t(i));
        else
            ST(i,j) = T_qw(t(i));
        end

    end

end

disp("Uninfected cells sensivity matrix: ")
disp(ST)

% Tuning method
for i = 1:np
    a = norm(ST(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalues
A1 = ST'*ST;
[V1,E1] = eig(A1);
disp("ST Eigenvectors Matrix:")
disp(V1) %eigen vector matrix
disp("ST Eigenvalues Matrix:")
disp(E1) %eigen values matrix

%Correlation method
corST = corr(A1);
disp("ST Correlation Matrix:")
disp(corST)


% Tm: Mutant-infected Cells

STm = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            STm(i,j) = Tm_p(t(i));
        elseif j == 2
            STm(i,j) = Tm_pm(t(i));
        elseif j == 3
            STm(i,j) = Tm_pw(t(i));
        elseif j == 4
            STm(i,j) = Tm_pmw(t(i));
        elseif j == 5
            STm(i,j) = Tm_km(t(i));
        elseif j == 6
            STm(i,j) = Tm_kw(t(i));
        elseif j == 7
            STm(i,j) = Tm_kR(t(i));
        elseif j == 8
            STm(i,j) = Tm_qm(t(i));
        else
            STm(i,j) = Tm_qw(t(i));
        end

    end

end

disp("Mutant-infected cells sensivity matrix: ")
disp(STm)

% Tuning method
for i = 1:np
   a = norm(STm(:,i));
   fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalues
A2 = STm'*STm;
[V2,E2] = eig(A2);
disp("STm Eigenvectors Matrix:")
disp(V2) %eigen vector matrix
disp("STm Eigenvalues Matrix:")
disp(E2) %eigen values matrix

%Correlation method
corSTm = corr(A2);
disp("STm Correlation Matrix:")
disp(corSTm)


% Tw: Wild type Infected Cells

STw = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            STw(i,j) = Tw_p(t(i));
        elseif j == 2
            STw(i,j) = Tw_pm(t(i));
        elseif j == 3
            STw(i,j) = Tw_pw(t(i));
        elseif j == 4
            STw(i,j) = Tw_pmw(t(i));
        elseif j == 5
            STw(i,j) = Tw_km(t(i));
        elseif j == 6
            STw(i,j) = Tw_kw(t(i));
        elseif j == 7
            STw(i,j) = Tw_kR(t(i));
        elseif j == 8
            STw(i,j) = Tw_qm(t(i));
        else
            STw(i,j) = Tw_qw(t(i));
        end

    end

end

disp("Wild type infected cells sensivity matrix: ")
disp(STw)

% Tuning method
for i = 1:np
    a = norm(STw(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalues
A3 = STw'*STw;
[V3,E3] = eig(A3);
disp("STw Eigenvectors Matrix:")
disp(V3) %eigen vector matrix
disp("STw Eigenvalues Matrix:")
disp(E3) %eigen values matrix

%Correlation method
corSTw = corr(A3);
disp("STw Correlation Matrix:")
disp(corSTw)

% Tmw: Cells infected by mutant viruses and wild-type viruses

STmw = zeros(lent, np);

for i = 1:lent
    for j = 1:np
        if j == 1
            STmw(i,j) = Tmw_p(t(i));
        elseif j == 2
            STmw(i,j) = Tmw_pm(t(i));
        elseif j == 3
            STmw(i,j) = Tmw_pw(t(i));
        elseif j == 4
            STmw(i,j) = Tmw_pmw(t(i));
        elseif j == 5
            STmw(i,j) = Tmw_km(t(i));
        elseif j == 6
            STmw(i,j) = Tmw_kw(t(i));
        elseif j == 7
            STmw(i,j) = Tmw_kR(t(i));
        elseif j == 8
            STmw(i,j) = Tmw_qm(t(i));
        else
            STmw(i,j) = Tmw_qw(t(i));
        end

    end

end

disp("Mutant and wt infected cells sensivity matrix: ")
disp(STmw)

% Tuning method
for i = 1:np
    a = norm(STmw(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalues
A4 = STmw'*STmw;
[V4,E4] = eig(A4);
disp("STmw Eigenvectors Matrix:")
disp(V4) %eigen vector matrix
disp("STmw Eigenvalues Matrix:")
disp(E4) %eigen values matrix

%Correlation method
corSTmw = corr(A4);
disp("STmw Correlation Matrix:")
disp(corSTmw)

% Mixed matrix
SS = zeros(4*lent, np);

%Matrix:
for i = 1:lent
    SS(4*i-3, :) = ST(i, :);
    SS(4*i-2, :) = STm(i, :);
    SS(4*i-1, :) = STw(i, :);
    SS(4*i, :) = STmw(i, :);
end


disp("Mixed sensitivity matrix: ")
disp(SS)

% Tuning method
for i = 1:np
    a = norm(SS(:,i));
    fprintf("Column %1.0f norm = %3.7f \n", i, a)
end

% Eigenvalue method
A5 = SS'*SS;
[V5,E5] = eig(A5);
disp("SS Eigenvectors Matrix:")
disp(V5) %eigen vector matrix
disp("SS Eigenvalues Matrix:")
disp(E5) %eigen values matrix

%Correlation method
corSS = corr(A5);
disp("SS Correlation Matrix:")
disp(corSS)



end