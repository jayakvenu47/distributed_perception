clc;
clear;
close all;

%% Target trajectray
R = 8;
ndt = 10;
dt = 1/ndt;

T = linspace(0, 2*pi, 200);
%ref_pos = [R*cos(T); sqrt(2) * R* (sin(T).* cos(T)); T];
ref_pos = [R*cos(T); R* sin(T); T];
N = length(T);
%% Define robots
init_pos = [0 0 -pi;
            10 0 -pi;
            -10 0 -pi;
            0 10 -pi;
            0 -10 -pi];
N_robot = 5;
robot = cell(N_robot);

%% Target dynamic
A = [1,0,dt,0; ... 
    0,1,0,dt; ...
    0,0,1,0; ...
    0,0,0,1];

G = [dt^2/2,0; ...
    0,dt^2/2; ...
    dt,0; ...
    0,dt];

Q = diag([1, 1]);

sv_sensor = 0.4;
R_sensor = sv_sensor^2 * eye(N_robot);

range_sensor = cell(N_robot,1);
for i = 1:N_robot
    range_sensor{i} = ['norm([x - ' num2str(init_pos(i,1)) ', y - ' num2str(init_pos(i,2)) '])'];
end

x_sym = str2sym({'x'; 'y' ; 'vx'; 'vy'});
h_sym = str2sym(range_sensor);
H_sym = jacobian(h_sym, x_sym);

h_func = matlabFunction(h_sym, 'Vars', {x_sym});
H_func = matlabFunction(H_sym, 'Vars', {x_sym});

%% EKF data
X_ekf = zeros(height(x_sym), N);
X_ekf(:,1) = [R, 0, 0, 0]; %2*pi*R*ndt/N

P_ekf = zeros(height(x_sym), height(x_sym), N);
P_ekf(:,:,1) = diag([0.1, 0.1, 0.1, 0.1]);

%% UKF data
X_ukf = zeros(height(x_sym), N);
X_ukf(:,1) = X_ekf(:,1);

P_ukf = zeros(height(x_sym), height(x_sym), N);
P_ukf(:,:,1) = diag([0.1, 0.1, 0.1, 0.1]);
D = 4;

%% Plot settings
f3=figure;
xlabel('x (m)')
ylabel('y (m)')
grid on
axis equal
xlim([-12 12])
ylim([-12 12])

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

traj = animatedline('color','b');
est_traj = animatedline('color','r');
ekf_traj = animatedline('color','g');
ukf_traj = animatedline('color','c');
tr = 0.3;
target_patch = rectangle('Curvature', [1 1], 'Position', [-tr, -tr, tr*2, tr*2], ...
        'facecolor', [205 132 85]/255, 'edgecolor', 'none');

pos_mesurement = cell(N_robot);
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
T_ekf = [];
T_ukf = [];
T_mt = [];
for k = 1:N-1
    x = ref_pos(1,k);
    y = ref_pos(2,k);
    theta = ref_pos(3,k);
    
    % plot trajectory and positions
    addpoints(traj, x, y);
    target_patch.Position = [x-tr, y-tr, tr*2, tr*2];
    
    X=[];
    Z=[];
    rho = [];
    for i = 1:N_robot
        r_x = init_pos(i,1);
        r_y = init_pos(i,2);
        r_theta = init_pos(i,3) - pi/2;
        r_rho = r_x^2 + r_y^2;
        rotation_matrix = [
            cos(r_theta), -sin(r_theta), r_x;
            sin(r_theta), cos(r_theta), r_y;
            0 0 1 ];
        t_pos = robot{i}.robot_body*rotation_matrix';
        set(robot{i}.robot_handle, 'Vertices', t_pos(:, 1:2));
        

        [d1, d2] = DIST([r_x, r_y], [x, y]);
        %d1 = d1 + sv_sensor * randn;
        X = [X; r_x, r_y];
        Z = [Z, d1];
        rho = [rho; r_rho];
    end
    Z = Z + mvnrnd(zeros(1, height(R_sensor)), R_sensor);
    
    %% Multilateration
    tic;
    H_mt = [];
    b_mt = [];
    for i = 2:N_robot
        H_mt=[H_mt; 2*(X(i,1)-X(1,1)), 2*(X(i,2)-X(1,2))];
        b_mt=[b_mt; Z(1)^2 - Z(i)^2 + rho(i) - rho(1)];
    end
    Est_mt = inv(H_mt'*H_mt)*H_mt'*b_mt;
    toc;
    T_mt = [T_mt, toc];
    
    addpoints(est_traj, Est_mt(1), Est_mt(2));
    err_mt(k) = norm([Est_mt(1)-x, Est_mt(2)-y]);
    err_mt_x(k) = Est_mt(1)-x;
    err_mt_y(k) = Est_mt(2)-y;

    
    %% EKF
    tic;
    X_ekf_pred = A * X_ekf(:,k);
    P_ekf_pred = A * P_ekf(:,:,k) * A' + G*Q*G';
    
    H_ekf = H_func(X_ekf_pred);
    S_ekf = H_ekf * P_ekf_pred * H_ekf' + R_sensor;
    W_ekf = P_ekf_pred * H_ekf' / S_ekf;
    I = Z' - h_func(X_ekf_pred);
    X_ekf(:, k+1) = X_ekf_pred + W_ekf * I;
    P_ekf(:,:,k+1) = P_ekf_pred - W_ekf*H_ekf*P_ekf_pred;
    
    toc;
    T_ekf = [T_ekf, toc];
    
    addpoints(ekf_traj, X_ekf(1, k+1), X_ekf(2, k+1));
    err_ekf(k) = norm([X_ekf(1, k+1)-x, X_ekf(2, k+1)-y]);
    err_ekf_x(k) = X_ekf(1, k+1)-x;
    err_ekf_y(k) = X_ekf(2, k+1)-y;

    
    %% UKF
    tic;
    X_ukf_pred = A * X_ukf(:,k);
    P_ukf_pred = A * P_ukf(:,:,k) * A' + G*Q*G';
    
    A_cho = chol(P_ukf_pred);
    SigmaPoints = zeros(D, 2*D);
    Z_SP = zeros(N_robot, 2*D);
    for i=1:D
        SigmaPoints(:,i) = X_ukf_pred + sqrt(D)*A_cho(i,:)';
        SigmaPoints(:,i+D) = X_ukf_pred - sqrt(D)*A_cho(i,:)';
        
        Z_SP(:,i) = h_func(SigmaPoints(:,i));
        Z_SP(:,i+D) = h_func(SigmaPoints(:,i+D));
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
    
    toc;
    T_ukf = [T_ukf, toc];
    
    addpoints(ukf_traj, X_ukf(1, k+1), X_ukf(2, k+1));
    err_ukf(k) = norm([X_ukf(1, k+1)-x, X_ukf(2, k+1)-y]);
    err_ukf_x(k) = X_ukf(1, k+1)-x;
    err_ukf_y(k) = X_ukf(2, k+1)-y;

    
    pause(0.001)
end

figure
hold on;box on;
xlabel('Time Step')
ylabel('Error Distance')
plot(err_mt, '-k.');
plot(err_ekf, '-b.');
plot(err_ukf, '-c.');
legend('Mt','EKF','UKF');

mean_err_mt = mean(err_mt);
mean_err_ekf = mean(err_ekf);
mean_err_ukf = mean(err_ukf);
err_mt_m = mean([err_mt_x; err_mt_y], 2);
err_mt_v = var([err_mt_x; err_mt_y], 0, 2);
err_ekf_m = mean([err_ekf_x; err_ekf_y], 2);
err_ekf_v = var([err_ekf_x; err_ekf_y], 0, 2);
err_ukf_m = mean([err_ukf_x; err_ukf_y], 2);
err_ukf_v = var([err_ukf_x; err_ukf_y], 0, 2);
T_mt_m = mean(T_mt);
T_ekf_m = mean(T_ekf);
T_ukf_m = mean(T_ukf);

fprintf('distance err, mt:%f, ekf:%f, ukf:%f\n', mean_err_mt,mean_err_ekf,mean_err_ukf);
fprintf('err_mt_mean, x:%f, y:%f\n', err_mt_m(1),err_mt_m(2));
fprintf('err_mt_var, x:%f, y:%f\n', err_mt_v(1),err_mt_v(2));
fprintf('err_ekf_mean, x:%f, y:%f\n', err_ekf_m(1),err_ekf_m(2));
fprintf('err_ekf_var, x:%f, y:%f\n', err_ekf_v(1),err_ekf_v(2));
fprintf('err_ukf_mean, x:%f, y:%f\n', err_ukf_m(1),err_ukf_m(2));
fprintf('err_ukf_var, x:%f, y:%f\n', err_ukf_v(1),err_ukf_v(2));
fprintf('time, mt:%f, ekf:%f, ukf:%f\n', T_mt_m, T_ekf_m,T_ukf_m);

fprintf('%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n',...
        mean_err_mt,mean_err_ekf,mean_err_ukf,err_mt_m(1),err_ekf_m(1),err_ukf_m(1),...
        err_mt_v(1),err_ekf_v(1),err_ukf_v(1),err_mt_m(2),err_ekf_m(2),err_ukf_m(2),...
        err_mt_v(2),err_ekf_v(2),err_ukf_v(2), T_mt_m, T_ekf_m,T_ukf_m);


function [dist, dist2]=DIST(A,B)
dist2 = (A(1)-B(1))^2 + (A(2)-B(2))^2;
dist = sqrt(dist2);
end
