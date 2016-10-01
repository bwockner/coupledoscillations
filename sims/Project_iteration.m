run Project

%%%Initial conditions
x0 = 0.5;
x_dot0 = 0;

theta0 = 0.4;
theta_dot0 = 0;

%%% time steps
delta_t = 0.001;
end_t = 5;                      % Must be a multiple of delta_t
n_entries = end_t/delta_t+1;  % Number of plottable entries

t_plot = [0:n_entries-1]*delta_t;

disp_x = zeros(1,n_entries);
vel_x = zeros(1,n_entries);
%accel_x = zeros(1,n_entries);

disp_theta = zeros(1,n_entries);
vel_theta = zeros(1,n_entries);
%accel_theta = zeros(1,n_entries);

disp_x(1) = x0;
vel_x(1) = x_dot0;
disp_theta(1) = theta0;
vel_theta(1) = theta_dot0;



t = 0;
while t<end_t;
    t_pos = round(t/delta_t + 1);
        
    %%%Solve simultaneous eq1/eq2 for accelerations
    eq1s = subs(eq1, [x, x_dot, theta, theta_dot], [disp_x(t_pos), vel_x(t_pos), disp_theta(t_pos), vel_theta(t_pos)]) == 0;
    eq2s = subs(eq2, [x, x_dot, theta, theta_dot], [disp_x(t_pos), vel_x(t_pos), disp_theta(t_pos), vel_theta(t_pos)]) == 0;

    new_accel = solve([eq1s,eq2s], [x_Ddot, theta_Ddot]);
    
    disp_x(t_pos+1) = disp_x(t_pos) + vel_x(t_pos)*delta_t + 0.5*delta_t^2*new_accel.x_Ddot;
    disp_theta(t_pos+1) = disp_theta(t_pos) + vel_theta(t_pos)*delta_t + 0.5*delta_t^2*new_accel.theta_Ddot;
    
    vel_x(t_pos+1) = vel_x(t_pos) + delta_t*new_accel.x_Ddot;
    vel_theta(t_pos+1) = vel_theta(t_pos) + delta_t*new_accel.theta_Ddot;
    
    t = t + delta_t;
end

dlmwrite('DATA.txt', t_plot);
dlmwrite('DATA.txt', disp_x, '-append')
dlmwrite('DATA.txt', disp_theta, '-append')
