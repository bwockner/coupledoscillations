calc_time_Total = tic;
calc_time_Project = tic;
%%% Calculate lagrange equations
run Project

fprintf('calc_time_Project = %f seconds \n', toc(calc_time_Project))

calc_time_Setup = tic;
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

delta_t = 0.001;
end_t = 1;                      % Must be a multiple of delta_t
n_entries = end_t/delta_t+1;    % Number of plottable entries

t_plot = [0:n_entries-1]*delta_t;

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

fprintf('calc_time_Setup = %f seconds \n', toc(calc_time_Setup))

%%%% Loop which solves for accelerations using lagranges eqns based on
%%%% previous time steps displacements and velocities, the subs into uvats
%%%% to calculate new displacements and velocities

for t = [1:n_entries-1];
    tic
    
    %  =====  OLD EQNS  =====
    %%%Solve simultaneous eq1/eq2 for accelerations
    %eq1s = subs(eq1, [x, x_dot, theta, theta_dot], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]) == 0;
    %eq2s = subs(eq2, [x, x_dot, theta, theta_dot], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]) == 0;
    %new_accel = solve([eq1s,eq2s], [x_Ddot, theta_Ddot]);

    %%% Subs values into generic eqns 
    new_x_Ddot = double(subs(new_accel.x_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]));
    new_theta_Ddot = double(subs(new_accel.theta_DdotV, [xV, x_dotV, thetaV, theta_dotV], [disp_x(t), vel_x(t), disp_theta(t), vel_theta(t)]));

    %%% uvats for x and theta
    disp_x(t+1) = disp_x(t) + vel_x(t)*delta_t + 0.5*delta_t^2*new_x_Ddot;
    disp_theta(t+1) = disp_theta(t) + vel_theta(t)*delta_t + 0.5*delta_t^2*new_theta_Ddot;

    vel_x(t+1) = vel_x(t) + delta_t*new_x_Ddot;
    vel_theta(t+1) = vel_theta(t) + delta_t*new_theta_Ddot;
    
    calc_time_Iteration(t) = toc;
end

fprintf('calc_time_Total = %f seconds \n', toc(calc_time_Total))

%%% Export data to text file for visualisation
dlmwrite('DATA.txt', t_plot);
dlmwrite('DATA.txt', disp_x, '-append')
dlmwrite('DATA.txt', disp_theta, '-append')
