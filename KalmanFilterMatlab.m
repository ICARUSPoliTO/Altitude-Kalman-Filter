%Both for Altitude and Vertical speed
clear
close all
AltData = importdata("/Users/denizakyuz/Desktop/ICARUS/AltimeterFilter/ALT.csv");
AltData = AltData.data(:,2);
p0 = 1013.25;
AltData = 44330 * (1 - (AltData ./ p0).^0.190263) - 296;
%AltDataX = linspace(0,3*2*pi,10000);
%AltData = sin(AltDataX)*50 + randn(1,length(AltDataX))*2;


% ==========================================================
% Simple 2-state Kalman filter for MS5607 altitude signal
% States: x = [altitude; vertical_velocity]
% No MATLAB toolboxes required, easy to port to C
% ==========================================================
dt = 1/89;                     % sampling period [s] (adjust to your data)

%Parameters to be adjusted
R = 0.4;                       % measurement noise variance (m^2) defult 0.4
q_pos = 1e-5;                  % process noise (altitude) default 1e-5
q_vel = 1e-3;                  % process noise (velocity) default 1e-3

% --- LOAD DATA ---
z = AltData;
N = length(z);

% --- INITIALIZATION ---
x = zeros(2, N);                % state history [alt; vel]
P = eye(2);                     % initial covariance
Q = diag([q_pos, q_vel]);       % process noise matrix
H = [1 0];                      % measurement matrix
I = eye(2);                     % identity

% --- MAIN KF LOOP ---
for k = 2:N
    % Predict step
    A = [1 dt; 0 1];
    x_pred = A * x(:,k-1);
    P_pred = A * P * A' + Q;

    % Update step
    z_meas = z(k);
    y = z_meas - H * x_pred;               % innovation
    S = H * P_pred * H' + R;               % innovation covariance
    K = (P_pred * H') / S;                 % Kalman gain (2x1)
    x(:,k) = x_pred + K * y;               % updated state
    P = (I - K * H) * P_pred;              % updated covariance
end

% --- RESULTS ---
alt_est = x(1,:);
vel_est = x(2,:);

% --- PLOTS ---
t = (0:N-1)*dt;
figure;
yyaxis left
plot(t(1000:end-1000), z(1000:end-1000), 'k:', t(1000:end-1000), alt_est(1000:end-1000), 'b-',LineWidth=2);
%plot(t(1000:end-1000), alt_est(1000:end-1000), 'b-');
xlabel('Time [s]'); ylabel('Altitude [m]');
legend('Measured','KF estimate'); grid on;

%yyaxis right
hold on
plot(t(1000:end-1000), vel_est(1000:end-1000), 'r',LineWidth=2);
xlabel('Time [s]'); ylabel('Vertical velocity [m/s]');
grid minor;

% AltData = importdata("GPS.csv");
% AltData = AltData.data(2000:end,12) - 296;
% figure(2);
% plot((AltData(1:end-1)-AltData(2:end))*25);

alt_est = alt_est(25000:end-20000);
vel_est = vel_est(25000:end-20000);

for i=1:length(alt_est/10);
    disp("Altitude: " + round(alt_est(i*10),1) + "m      Vertical speed: " + round(vel_est(i),2) + "m/s")
    %pause(1/90*10);
end
