% Unknown variables - Symbolic Functions in t
syms theta(t) theta_dot(t) theta_Ddot(t) x(t) x_dot(t) x_Ddot(t)

% Known Variables
%syms M_t R_t R_tg I_b M_b R_bg R_ti r_b g

M_t = 1.25;    %Mass of the track n kg
R_t = 0.6;    %Outside radius of the track in m
R_ti = 0.5;   %Inside radius of the track in m
R_tg = 0.4;   %Radius of the track centre of mass in m    = (4*R_t^2+4*R_t*R_ti+4*R_ti^2)/(3*pi*(R_t+R_ti))
M_b = 0.25;    %Mass of the ball in kg
r_b = 0.05;    %Radius of the ball in m
I_b = 2/5 * M_b * r_b^2; %Moment of inertia in kg m^2
R_bg = R_ti - r_b;   %Radius from centre of track to centre of ball in m
g = 9.81;  %accel due to gravity in m/s^2


%%%% Equations of Motion
%%%% x is left-right, y is up-down
%%%% Second commented out equations are simplified using small angle approx
%%% Track
x_t = (R_t - R_tg) * sin(x(t)/R_t);
x_t_dot = (R_t - R_tg) * x_dot(t)/R_t * cos(x(t)/R_t);
%x_t = (R_t - R_tg) * x/R_t;
%x_t_dot = (R_t - R_tg) * x_dot /R_t;

y_t = R_tg * (1- cos(x(t)/R_t));
y_t_dot = R_tg / R_t * x_dot(t) *sin(x(t)/R_t);
%y_t = R_tg * (1- cos(x/R_t));
%y_t_dot = R_tg / R_t * x_dot * x/R_t;

%%% Ball
x_b = x + R_bg*sin(theta(t));
x_b_dot = x_dot + R_bg*theta_dot(t)* cos(theta(t));
%x_b = x + R_bg*theta;
%x_b_dot = x_dot + R_bg*theta_dot;

y_b = R_bg * (1- cos(theta(t)));
y_b_dot = R_bg * theta_dot(t) * sin(theta(t));
%y_b = R_bg * (1- cos(theta));
%y_b_dot = R_bg * theta_dot * theta;

%%% Kinetic and potential energies
T = 1/2 * M_t * (x_t_dot^2 + y_t_dot^2) + 1/2 * M_b * (x_b_dot^2 + y_b_dot^2) + 1/2 * I_b *(R_ti/r_b * (theta_dot(t) + (x_dot(t)/R_t)))^2;
V = M_t*g*y_t + M_b*g*y_b;

% ===== OLD EQNS =====
%eq1 = subs(diff(T-V, x_dot), [x_dot, theta_dot], [x_Ddot, theta_Ddot]) - diff(T-V, x);          %LAGRANGE for x
%eq2 = subs(diff(T-V, theta_dot), [x_dot, theta_dot], [x_Ddot, theta_Ddot]) - diff(T-V, theta);  %LAGRANGE for theta


% Updated to find Lagrange Equations exactly, no approximations made.
eq1a = diff(functionalDerivative(T-V, x_dot),t) - functionalDerivative(T-V, x);              %LAGRANGE for x
eq2a = diff(functionalDerivative(T-V, theta_dot), t) - functionalDerivative(T-V, theta);     %LAGRANGE for theta

% Swap notation for derviatives from diff(x(t), t) to x_dot etc.
eq1 = subs(eq1a, [diff(x(t), t), diff(x_dot(t), t), diff(theta(t), t), diff(theta_dot(t), t)], [x_dot(t), x_Ddot(t), theta_dot(t), theta_Ddot(t)]);
eq2 = subs(eq2a, [diff(x(t), t), diff(x_dot(t), t), diff(theta(t), t), diff(theta_dot(t), t)], [x_dot(t), x_Ddot(t), theta_dot(t), theta_Ddot(t)]);

%%%%Manually solved linear matrices using small angle approx/etc.
M = [M_t + 2*R_tg*M_t/R_t + M_t*R_tg^2/R_t^2 + M_b + I_b*R_ti^2/(R_t^2*r_b^2),    I_b*R_ti^2/(R_t*r_b^2)+M_b*R_bg;
    M_b * R_bg + I_b*R_ti^2/(R_t*r_b^2),                                         M_b*R_bg^2+I_b*R_ti^2/r_b^2     ];

K = [g*M_t*R_tg/R_t^2, 0;
    0, g * M_b * R_bg];
