%%%%%%%% Runge-Kutta Method
%%%%%%%% Seems to be less efficient than straight iteration at this stage
%%% Calculate lagrange equations
run Project

%%%Initial conditions
x0 = 0;
x_dot0 = 0;

theta0 = 1;
theta_dot0 = 0;

%%% time steps
%%%% NOTES
% dt >= 0.1 essentially useless. Way too much error
% Errors propogate quickly for dt > around 0.01.
% Even for 0.002 they can be seen over a interval 5-10 sec long

h = 0.002;
end_t = 5;                      % Must be a multiple of delta_t
n_entries = end_t/h+1;          % Number of plottable entries

t_plot = [0:n_entries-1]*h;

%%% Empty lists for data to be stored in
disp_x = zeros(1,n_entries);
vel_x = zeros(1,n_entries);
%accel_x = zeros(1,n_entries);

disp_theta = zeros(1,n_entries);
vel_theta = zeros(1,n_entries);
%accel_theta = zeros(1,n_entries);

calc_time_Iteration = zeros(1,n_entries);

% Store initial conditions
disp_x(1) = x0;
vel_x(1) = x_dot0;
disp_theta(1) = theta0;
vel_theta(1) = theta_dot0;


% ==== Have to swap symbolic functions to symbolic variables for solve function ====
syms thetaV theta_dotV theta_DdotV xV x_dotV x_DdotV
eq1q = subs(eq1, [theta, theta_dot, theta_Ddot, x, x_dot, x_Ddot], [thetaV, theta_dotV, theta_DdotV, xV, x_dotV, x_DdotV,]);
eq2q = subs(eq2, [theta, theta_dot, theta_Ddot, x, x_dot, x_Ddot], [thetaV, theta_dotV, theta_DdotV, xV, x_dotV, x_DdotV,]);

% Generic equations for acceleration updates 
new_accel = solve([eq1q == 0,eq2q == 0], [x_DdotV, theta_DdotV]);

%fprintf('calc_time_Setup = %f seconds \n', toc(calc_time_Setup))


for t = [1:n_entries-1];

    t
    
    a_x = vel_x(t);
    a_x_dot = double(subs(new_accel.x_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]));
    a_theta = vel_theta(t);
    a_theta_dot = double(subs(new_accel.theta_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]));

    b_x = vel_x(t) + h/2 * a_x;
    b_x_dot = double(subs(new_accel.x_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h/2 * a_x, vel_x(t) + h/2 * a_x_dot, disp_theta(t)+ h/2 * a_theta, vel_theta(t)+ h/2 * a_theta_dot]));
    b_theta = vel_theta(t) + h/2 * a_theta;
    b_theta_dot = double(subs(new_accel.theta_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h/2 * a_x, vel_x(t) + h/2 * a_x_dot, disp_theta(t)+ h/2 * a_theta, vel_theta(t)+ h/2 * a_theta_dot]));

    c_x = vel_x(t) + h/2 * b_x;
    c_x_dot = double(subs(new_accel.x_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h/2 * b_x, vel_x(t) + h/2 * b_x_dot, disp_theta(t)+ h/2 * b_theta, vel_theta(t)+ h/2 * b_theta_dot]));
    c_theta = vel_theta(t) + h/2 * b_theta;
    c_theta_dot = double(subs(new_accel.theta_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h/2 * b_x, vel_x(t) + h/2 * b_x_dot, disp_theta(t)+ h/2 * b_theta, vel_theta(t)+ h/2 * b_theta_dot]));

    d_x = vel_x(t) + h * c_x;
    d_x_dot = double(subs(new_accel.x_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h * c_x, vel_x(t) + h * c_x_dot, disp_theta(t)+ h * c_theta, vel_theta(t)+ h * c_theta_dot]));
    d_theta = vel_theta(t) + h * c_theta;
    d_theta_dot = double(subs(new_accel.theta_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t) + h * c_x, vel_x(t) + h * c_x_dot, disp_theta(t)+ h * c_theta, vel_theta(t)+ h * c_theta_dot]));

    disp_x(t+1) = disp_x(t) + h/6 * (a_x + 2* b_x + 2*c_x + d_x);
    vel_x(t+1) = vel_x(t) + h/6 * (a_x_dot + 2* b_x_dot + 2*c_x_dot + d_x_dot);

    disp_theta(t+1) = disp_theta(t) + h/6 * (a_theta + 2* b_theta + 2*c_theta + d_theta);
    vel_theta(t+1) = vel_theta(t) + h/6 * (a_theta_dot + 2* b_theta_dot + 2*c_theta_dot + d_theta_dot);

end