% Unknown variables
syms theta theta_dot theta_Ddot x x_dot x_Ddot

% Known Variables
syms M_t R_t R_tg I_b M_b R_bg R_ti r_b g

% M_t =     %Mass of the track n kg
% R_t =     %Outside radius of the track in m
% R_tg =    %Radius of the centre of mass in m
% M_b =     %Mass of the ball in kg
% r_b =     %Radius of the ball in m
% I_b = 2/5 * M_b * r_b^2 %Moment of inertia in kg m^2
% R_bg =    %Radius from centre of track to centre of  ball in m
% R_ti =    $Inside radius of the track in m
% g = 9.81  %accel due to gravity in m/s^2


%%%% Equations of Motion
%%%% x is left-right, y is up-down
%%%% Second commented out equations are simplified using small angle approx
%%% Track
x_t = (R_t - R_tg) * sin(x/R_t);
x_t_dot = (R_t - R_tg) * x_dot /R_t * cos(x/R_t);
%x_t = (R_t - R_tg) * x/R_t;
%x_t_dot = (R_t - R_tg) * x_dot /R_t;

y_t = R_tg * (1- cos(x/R_t));
y_t_dot = R_tg / R_t * x_dot *sin(x/R_t);
%y_t = R_tg * (1- cos(x/R_t));
%y_t_dot = R_tg / R_t * x_dot * x/R_t;

%%% Ball
x_b = x + R_bg*sin(theta);
x_b_dot = x_dot + R_bg*theta_dot* cos(theta);
%x_b = x + R_bg*theta;
%x_b_dot = x_dot + R_bg*theta_dot;

y_b = R_bg * (1- cos(theta));
y_b_dot = R_bg * theta_dot * sin(theta);
%y_b = R_bg * (1- cos(theta));
%y_b_dot = R_bg * theta_dot * theta;

%%% Kinetic and potential energies
T = 1/2 * M_t * (x_t_dot^2 + y_t_dot^2) + 1/2 * M_b * (x_b_dot^2 + y_b_dot^2) + 1/2 * I_b *(R_ti/r_b * (theta_dot - x_dot/R_t));
V = M_t*g*x_t + M_b*g*x_b;

%eq1 = diff(T, x_dot) + diff(V,x);
%eq2 = diff(T, theta_dot) + diff(V,theta);