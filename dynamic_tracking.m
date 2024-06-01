clc;
clear;
close all;

%% Target trajectray and senario parameters
% time parameters
ndt = 100;
ndt_gps = 10;
ndt_ranging = 10;
dt = 1/ndt;
dt_range = ndt_ranging * dt;
preparation = 8;
duration = 50;

% robot control parameters
gain_v = 3;
gain_omega = 2;
max_v = 1.8;
max_o = 2;

% time series
time = (0:dt:duration)';
pre_time = (0:dt:preparation)';
N_pre = length(pre_time);
N = length([pre_time; time]);

% target trajectory
R = 8;
K_theta = 2*pi/duration;
T = time' * K_theta;
mov_pos = [R*cos(T); sqrt(2) * R* (sin(T).* cos(T))];
pre_pos = repmat(mov_pos(:,1),1, N_pre);
ref_pos = [pre_pos, mov_pos];

%% Robot Dead reckoning, GPS and IMU problems
rob_x_sym = str2sym({'x'; 'y'; 'theta'; 'v'});
rob_u_real_sym = str2sym({'omega'; 'v'}); % for calculating true value
rob_u_imu_sym = str2sym({'omega'; 'a'}); % for prediction kalman filter

% for calculating true value
rob_f_real_sym = str2sym({ 
    'x + v * dt * cos(theta)';...
    'y + v * dt * sin(theta)';...
    'theta + omega * dt';...
    'v';
    });

% for prediction using gyro and accelerator (omega & a)
rob_f_sym = str2sym({ 
    'x + v * dt * cos(theta)';...
    'y + v * dt * sin(theta)';...
    'theta + omega * dt';...
    'v + a*dt';
    });

% for updating x,y using gps
rob_h_gps_sym = str2sym({
    'x';...
    'y';
    });

% for updating theta,v using odometry and magnetometer
rob_h_odma_sym = str2sym({
    'theta';...
    'v';
    });

% linearization
rob_A_sym = jacobian(rob_f_sym, rob_x_sym);
rob_G_sym = jacobian(rob_f_sym, rob_u_imu_sym);
rob_H_gps_sym = jacobian(rob_h_gps_sym, rob_x_sym);
rob_H_odma_sym = jacobian(rob_h_odma_sym, rob_x_sym);

% code generation
rob_f_func = matlabFunction(subs(rob_f_sym, 'dt', dt), 'Vars', {rob_u_imu_sym rob_x_sym});
rob_f_real_func = matlabFunction(subs(rob_f_real_sym, 'dt', dt), 'Vars', {rob_u_real_sym rob_x_sym});
rob_A_func = matlabFunction(subs(rob_A_sym, 'dt', dt), 'Vars', {rob_u_imu_sym rob_x_sym});
rob_G_func = matlabFunction(subs(rob_G_sym, 'dt', dt), 'Vars', {rob_u_imu_sym rob_x_sym});

rob_h_gps_func = matlabFunction(rob_h_gps_sym, 'Vars', {rob_x_sym});
rob_H_gps_func = matlabFunction(rob_H_gps_sym, 'Vars', {rob_x_sym});
rob_h_odma_func = matlabFunction(rob_h_odma_sym, 'Vars', {rob_x_sym});
rob_H_odma_func = matlabFunction(rob_H_odma_sym, 'Vars', {rob_x_sym});

% sensors unvertainty
% odometry IMU
Q_imu =  diag([0.01, 0.01]);
% measurement GPS
R_gps = diag([0.04, 0.04]);
% measurement encoder and magnetometer
R_odma = diag([0.04, 0.04]);

%% Define robots
N_robot = 5;
% initial positions
init_pos = [0 0 -pi/2;
            0 5 -pi/2;
            5 0 -pi/2;
            -5 0 -pi/2;
            0 -5 -pi/2];
% relative positons refer to the target
rel_pos = [ -0.1  1.5;
             1.2  0.8;
             1.2 -0.8 ;
            -1.5  0;
            -0.1 -1.5];

figure;
title('Multi-robot Target Tracking');

robot_state = cell(N_robot);
for r_i = 1:N_robot    
    robot_state{r_i}.x_real = zeros(height(rob_x_sym), N);
    robot_state{r_i}.x_real(:,1) = [init_pos(r_i,1); init_pos(r_i,2); init_pos(r_i,3); 0];
    robot_state{r_i}.acc_real = zeros(1,N);
    
    robot_state{r_i}.x = zeros(height(rob_x_sym), N);
    robot_state{r_i}.x(:,1) = robot_state{r_i}.x_real(:,1);
    robot_state{r_i}.P = zeros(height(rob_x_sym), height(rob_x_sym), N);
    robot_state{r_i}.P(:,:,1) = diag([0, 0, 0, 0]);

    robot_state{r_i}.dist_err = zeros(1, N);
    robot_state{r_i}.orie_err = zeros(1, N);

    robot_state{r_i}.ctrl_v = zeros(1,N);
    robot_state{r_i}.ctrl_omega = zeros(1,N);

    robot_state{r_i}.u = zeros(2,N-1)';
    robot_state{r_i}.z_gps = zeros(2,N-1)';
    robot_state{r_i}.z_odma = zeros(2,N-1)';
end

robot = cell(N_robot);
r_patch = robot_patch;
for i = 1:N_robot
    robot{i}.robot_body = r_patch.vertices;
    robot{i}.robot_handle = patch(...
    'Vertices', robot{i}.robot_body(:, 1:2), ...
    'Faces', r_patch.faces, ...
    'FaceColor', 'flat', ...
    'FaceVertexCData', r_patch.colors, ...
    'EdgeColor','none');
end

%% Target dynamic
A = [1,0,dt_range,0; ... 
    0,1,0,dt_range; ...
    0,0,1,0; ...
    0,0,0,1];

G = [dt_range^2/2,0; ...
    0,dt_range^2/2; ...
    dt_range,0; ...
    0,dt_range];

Q = diag([1, 1]);

sv_sensor = 0.1;
R_sensor = sv_sensor^2 * eye(N_robot);

range_sensor = cell(N_robot,1);
for i = 1:N_robot
    %range_sensor{i} = ['norm([x - ' num2str(init_pos(i,1)) ', y - ' num2str(init_pos(i,2)) '])'];
    range_sensor{i} = ['norm([x - x_' num2str(i) ', y - y_' num2str(i) '])'];
end

anchor_x = cell(2*N_robot, 1);
for i = 1:N_robot
    anchor_x{2 * i - 1} = ['x_' num2str(i)];
    anchor_x{2 * i} = ['y_' num2str(i)];
end

anchor_x_sym = str2sym(anchor_x);

x_sym = str2sym({'x'; 'y' ; 'vx'; 'vy'});
h_sym = str2sym(range_sensor);
H_sym = jacobian(h_sym, x_sym);

h_func = matlabFunction(h_sym, 'Vars', {anchor_x_sym x_sym});
H_func = matlabFunction(H_sym, 'Vars', {anchor_x_sym x_sym});

%% EKF data
X_ekf = zeros(height(x_sym), N);
X_ekf(:,1) = [R, 0, 0, 2*pi*R*ndt/N]; %2*pi*R*ndt/N

P_ekf = zeros(height(x_sym), height(x_sym), N);
P_ekf(:,:,1) = diag([0.1, 0.1, 0.1, 0.1]);

%% UKF data
X_ukf = zeros(height(x_sym), N);
X_ukf(:,1) = X_ekf(:,1);

P_ukf = zeros(height(x_sym), height(x_sym), N);
P_ukf(:,:,1) = diag([0.1, 0.1, 0.1, 0.1]);
D = 4;

%% Plot settings
xlabel('x (m)')
ylabel('y (m)')
grid on
axis equal
xlim([-12 12])
ylim([-10 15])

traj = animatedline('color','b');
est_traj = animatedline('color','r');
ekf_traj = animatedline('color','g');
ukf_traj = animatedline('color','c');
tr = 0.3;
target_patch = rectangle('Curvature', [1 1], 'Position', [-tr, -tr, tr*2, tr*2], ...
        'facecolor', [205 132 85]/255, 'edgecolor', 'none');
    
legend('robot1', 'robot2', 'robot3', 'robot4', 'robot5', ...
    'Truth','MulLa','EKF','UKF', 'Location','north','NumColumns', 3)

target_pos_calculated = zeros(2, N);
target_pos_calculated(:, 1) = [ref_pos(1,1), ref_pos(2,1)];
%% main loop
err_mt = zeros(1, N-1);
err_ekf = zeros(1, N-1);
err_ukf = zeros(1, N-1);
err_mt_x = zeros(1, N-1);
err_ekf_x = zeros(1, N-1);
err_ukf_x = zeros(1, N-1);
err_mt_y = zeros(1, N-1);
err_ekf_y = zeros(1, N-1);
err_ukf_y = zeros(1, N-1);

anchor_vec = zeros(2*N_robot,1);
for k = 1:N-1
    %% update target position
    x = ref_pos(1,k);
    y = ref_pos(2,k);
    
    % plot trajectory and positions
    addpoints(traj, x, y);
    target_patch.Position = [x-tr, y-tr, tr*2, tr*2];
    
    X=[];
    Z=[];
    rho = [];
    for i = 1:N_robot
        %% Robot Control
        % Calculate control error
        err = target_pos_calculated(:,k)' + rel_pos(i, :) - robot_state{i}.x(1:2,k)';
        [rtn_theta, rtn_rho] = cart2pol(err(1), err(2));
        robot_state{i}.dist_err(1,k) = rtn_rho;
        robot_state{i}.orie_err(1,k) = wrapToPi(rtn_theta - robot_state{i}.x(3,k));
        
        % Calculate control output
        if robot_state{i}.dist_err(1,k) < 0.1
            robot_state{i}.ctrl_v(1,k) = 0.0;
            robot_state{i}.ctrl_omega(1,k) = 0.0;
        else
            robot_state{i}.ctrl_v(1,k) = gain_v * robot_state{i}.dist_err(1,k);
            robot_state{i}.ctrl_omega(1,k) = gain_omega * robot_state{i}.orie_err(1,k);
        end
        
        robot_state{i}.ctrl_v(1,k) = min( max(robot_state{i}.ctrl_v(1,k), -max_v), max_v);
        robot_state{i}.ctrl_omega(1,k) = min( max(robot_state{i}.ctrl_omega(1,k), -max_o), max_o);
        
        %% Robot Mesurement
        % Calculate true state
        robot_state{i}.x_real(:,k+1) = rob_f_real_func([robot_state{i}.ctrl_omega(1,k); robot_state{i}.ctrl_v(1,k)], robot_state{i}.x_real(:,k));

        % Calculate mesurement
        % odometry IMU
        robot_state{i}.acc_real(:, k) = (robot_state{i}.x_real(4, k+1) - robot_state{i}.x_real(4, k))/dt;
        robot_state{i}.u(k,:) = [robot_state{i}.ctrl_omega(1,k), robot_state{i}.acc_real(:, k),] + mvnrnd(zeros(1, height(Q_imu)), Q_imu);
        % measurement GPS
        robot_state{i}.z_gps(k,:) = [robot_state{i}.x_real(1,k+1), robot_state{i}.x_real(2,k+1)] + mvnrnd(zeros(1, height(R_gps)), R_gps);
        % measurement encoder and magnetometer
        robot_state{i}.z_odma(k,:) = [robot_state{i}.x_real(3,k+1), robot_state{i}.x_real(4,k+1)] + mvnrnd(zeros(1, height(R_odma)), R_odma);
        
        %% Robot Estimation
         % Prediction
        robot_state{i}.x(:, k + 1) = rob_f_func(robot_state{i}.u(k, :)', robot_state{i}.x(:,k));
        rob_A = rob_A_func(robot_state{i}.u(k, :)', robot_state{i}.x(:, k));
        rob_G = rob_G_func(robot_state{i}.u(k, :)', robot_state{i}.x(:, k));
        robot_state{i}.P(:, :, k+1) = rob_A * robot_state{i}.P(:, :, k) * rob_A' + rob_G * Q_imu * rob_G';

        % Update encoder and magnetometer
        H_odma = rob_H_odma_func(robot_state{i}.x(:, k+1));
        S_odma = H_odma * robot_state{i}.P(:, :, k+1) * H_odma' + R_odma;
        W_odma = robot_state{i}.P(:, :, k+1) * H_odma' / S_odma;
        I_odma = robot_state{i}.z_odma(k, :)' - rob_h_odma_func(robot_state{i}.x(:, k+1));

        robot_state{i}.x(:, k+1) = robot_state{i}.x(:, k+1) + W_odma*I_odma;
        robot_state{i}.P(:, :, k+1) = robot_state{i}.P(:,:,k+1) - W_odma * H_odma * robot_state{i}.P(:, :, k+1);
        
        % Update GPS
        if mod(k, ndt_gps) == 0
            H_gps = rob_H_gps_func(robot_state{i}.x(:, k+1));
            S_gps = H_gps * robot_state{i}.P(:, :, k+1) * H_gps' + R_gps;
            W_gps = robot_state{i}.P(:, :, k+1) * H_gps' / S_gps;
            I_gps = robot_state{i}.z_gps(k, :)' - rob_h_gps_func(robot_state{i}.x(:, k+1));

            robot_state{i}.x(:, k+1) = robot_state{i}.x(:, k+1) + W_gps*I_gps;
            robot_state{i}.P(:, :, k+1) = robot_state{i}.P(:,:,k+1) - W_gps * H_gps * robot_state{i}.P(:, :, k+1);
        end
        
        % for kf calc
        anchor_vec(2*i -1, 1) = robot_state{i}.x(1, k+1);
        anchor_vec(2*i, 1) = robot_state{i}.x(2, k+1);
        
        %% Update robot patch
        r_x = robot_state{i}.x_real(1,k+1);
        r_y = robot_state{i}.x_real(2,k+1);
        r_theta = robot_state{i}.x_real(3,k+1) - pi/2;
        r_rho = r_x^2 + r_y^2;
        rotation_matrix = [
            cos(r_theta), -sin(r_theta), r_x;
            sin(r_theta), cos(r_theta), r_y;
            0 0 1 ];
        t_pos = robot{i}.robot_body*rotation_matrix';
        set(robot{i}.robot_handle, 'Vertices', t_pos(:, 1:2));
        
        [d1, d2] = DIST([r_x, r_y], [x, y]);
        X = [X; r_x, r_y];
        Z = [Z, d1];
        rho = [rho; r_rho];
    end
    Z = Z + mvnrnd(zeros(1, height(R_sensor)), R_sensor);
    
    %% Target Localisition
    if mod(k, ndt_ranging)==0
       %% Multilateration
        H_mt = [];
        b_mt = [];
        for i = 2:N_robot
            H_mt=[H_mt; 2*(X(i,1)-X(1,1)), 2*(X(i,2)-X(1,2))];
            b_mt=[b_mt; Z(1)^2 - Z(i)^2 + rho(i) - rho(1)];
        end
        Est_mt = inv(H_mt'*H_mt)*H_mt'*b_mt;

        err_mt(k) = norm([Est_mt(1)-x, Est_mt(2)-y]);
        err_mt_x(k) = Est_mt(1)-x;
        err_mt_y(k) = Est_mt(2)-y;

        %% EKF
        X_ekf_pred = A * X_ekf(:,k);
        P_ekf_pred = A * P_ekf(:,:,k) * A' + G*Q*G';

        H_ekf = H_func(anchor_vec, X_ekf_pred);
        S_ekf = H_ekf * P_ekf_pred * H_ekf' + R_sensor;
        W_ekf = P_ekf_pred * H_ekf' / S_ekf;
        I = Z' - h_func(anchor_vec, X_ekf_pred);
        X_ekf(:, k+1) = X_ekf_pred + W_ekf * I;
        P_ekf(:,:,k+1) = P_ekf_pred - W_ekf*H_ekf*P_ekf_pred;

        err_ekf(k) = norm([X_ekf(1, k+1)-x, X_ekf(2, k+1)-y]);
        err_ekf_x(k) = X_ekf(1, k+1)-x;
        err_ekf_y(k) = X_ekf(2, k+1)-y;

        %% UKF
        X_ukf_pred = A * X_ukf(:,k);
        P_ukf_pred = A * P_ukf(:,:,k) * A' + G*Q*G';

        A_cho = chol(P_ukf_pred);
        SigmaPoints = zeros(D, 2*D);
        Z_SP = zeros(N_robot, 2*D);
        for i=1:D
            SigmaPoints(:,i) = X_ukf_pred + sqrt(D)*A_cho(i,:)';
            SigmaPoints(:,i+D) = X_ukf_pred - sqrt(D)*A_cho(i,:)';

            Z_SP(:,i) = h_func(anchor_vec, SigmaPoints(:,i));
            Z_SP(:,i+D) = h_func(anchor_vec, SigmaPoints(:,i+D));
        end
        Zpred = mean(Z_SP, 2);

        Pzz = zeros(N_robot, N_robot);
        for i=1:2*D
            Pzz = Pzz + (Z_SP(:,i) - Zpred)*(Z_SP(:,i) - Zpred)'/(2*D);
        end
        Pzz = Pzz + R_sensor;

        Pxz = zeros(D, N_robot);
        for i=1:2*D
            Pxz = Pxz + (SigmaPoints(:,i) - X_ukf_pred)*(Z_SP(:,i) - Zpred)'/(2*D);
        end

        W_ukf = Pxz/Pzz;
        X_ukf(:,k+1) = X_ukf_pred + W_ukf*(Z' - Zpred);
        P_ukf(:,:,k+1) = P_ukf_pred - W_ukf*Pzz*W_ukf';
        err_ukf(k) = norm([X_ukf(1, k+1)-x, X_ukf(2, k+1)-y]);
        err_ukf_x(k) = X_ukf(1, k+1)-x;
        err_ukf_y(k) = X_ukf(2, k+1)-y;
        
        target_pos_calculated(1,k+1) = X_ekf(1, k+1);
        target_pos_calculated(2,k+1) = X_ekf(2, k+1);
    else
        X_ekf(:,k+1) = X_ekf(:,k);
        P_ekf(:,:,k+1) = P_ekf(:,:,k);
        
        X_ukf(:,k+1) = X_ukf(:,k);
        P_ukf(:,:,k+1) = P_ukf(:,:,k);
        
        target_pos_calculated(1,k+1) = target_pos_calculated(1,k);
        target_pos_calculated(2,k+1) = target_pos_calculated(2,k);
    end
    
    %% Plot traj
    if k > N_pre
        addpoints(est_traj, Est_mt(1), Est_mt(2));
        addpoints(ekf_traj, X_ekf(1, k+1), X_ekf(2, k+1));
        addpoints(ukf_traj, X_ukf(1, k+1), X_ukf(2, k+1));
    end
    
    pause(0.001)
end

%% plot results
figure
hold on;box on;
plot(err_mt, '-k.');
plot(err_ekf, '-b.');
plot(err_ukf, '-c.');
legend('MT','EKF','UKF');
xlabel('Time Step')
ylabel('Error Distance')


figure;
hold on; axis equal;
xlim([-10 10])
ylim([-8 12])
title('Robot1 Path')
plot(robot_state{1}.x(1, :),robot_state{1}.x(2, :), 'color', 'k');
plot(robot_state{1}.z_gps(1:10:end,1),robot_state{1}.z_gps(1:10:end,2), 'color', 'r');
plot(ref_pos(1,:) + rel_pos(1, 1), ref_pos(2,:) + rel_pos(1, 2), 'color', 'b');
legend('Estimated','GPS','Ref');


err_mt = err_mt(ndt_ranging: ndt_ranging :end);
err_ekf = err_ekf(ndt_ranging: ndt_ranging :end);
err_ukf = err_ukf(ndt_ranging: ndt_ranging :end);
err_mt_x = err_mt_x(ndt_ranging: ndt_ranging :end);
err_ekf_x = err_ekf_x(ndt_ranging: ndt_ranging :end);
err_ukf_x = err_ukf_x(ndt_ranging: ndt_ranging :end);
err_mt_y = err_mt_y(ndt_ranging: ndt_ranging :end);
err_ekf_y = err_ekf_y(ndt_ranging: ndt_ranging :end);
err_ukf_y = err_ukf_y(ndt_ranging: ndt_ranging :end);

mean_err_mt = mean(err_mt);
mean_err_ekf = mean(err_ekf);
mean_err_ukf = mean(err_ukf);
err_mt_m = mean([err_mt_x; err_mt_y], 2);
err_mt_v = var([err_mt_x; err_mt_y], 0, 2);
err_ekf_m = mean([err_ekf_x; err_ekf_y], 2);
err_ekf_v = var([err_ekf_x; err_ekf_y], 0, 2);
err_ukf_m = mean([err_ukf_x; err_ukf_y], 2);
err_ukf_v = var([err_ukf_x; err_ukf_y], 0, 2);

fprintf('distance err, mt:%f, ekf:%f, ukf:%f\n', mean_err_mt,mean_err_ekf,mean_err_ukf);
fprintf('err_mt_mean, x:%f, y:%f\n', err_mt_m(1),err_mt_m(2));
fprintf('err_mt_var, x:%f, y:%f\n', err_mt_v(1),err_mt_v(2));
fprintf('err_ekf_mean, x:%f, y:%f\n', err_ekf_m(1),err_ekf_m(2));
fprintf('err_ekf_var, x:%f, y:%f\n', err_ekf_v(1),err_ekf_v(2));
fprintf('err_ukf_mean, x:%f, y:%f\n', err_ukf_m(1),err_ukf_m(2));
fprintf('err_ukf_var, x:%f, y:%f\n', err_ukf_v(1),err_ukf_v(2));

fprintf('%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n',...
        mean_err_mt,mean_err_ekf,mean_err_ukf,err_mt_m(1),err_ekf_m(1),err_ukf_m(1),...
        err_mt_v(1),err_ekf_v(1),err_ukf_v(1),err_mt_m(2),err_ekf_m(2),err_ukf_m(2),...
        err_mt_v(2),err_ekf_v(2),err_ukf_v(2));

%% calculate distance
function [dist, dist2]=DIST(A,B)
dist2 = (A(1)-B(1))^2 + (A(2)-B(2))^2;
dist = sqrt(dist2);
end
