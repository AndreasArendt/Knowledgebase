clear;

std_accel = 1;
std_vel = 0.01;
s = generateTrajectory('tEnd', 30, 'std_accel', std_accel, 'std_vel', std_vel);

%% Factor Graph
% X := [pos, vel]
X0 = [0; 0];

nStates = numel(X0);
k = numel(s.timestamp);

F = [1  s.dt  ; ...
     0  1];

B = [0.5*s.dt^2; ...
     s.dt];

H = [0 1];

P0 = diag([1, 1]);
Q  = (std_accel^2) * (B * B');
R  = std_vel^2;

fac = estimation.factorgraph(X0, P0);

u = s.imu.a__mDs2;
z = s.mes.v__mDs;
for ii = 1:numel(s.timestamp)    
    fac.addStep(F, B, u(ii), Q, H=H,R=diag(R),z=z(ii));    
end

[x_est] = fac.solve();
% [x_est, sigma] = fac.solve();
% [x_est, sigma, res] = fac.solve();

%%
af = afigure;
subplot(3,1,1);
hold on; grid on;
plot(s.timestamp, s.imu.a__mDs2, 'DisplayName', 'accel [meas]')
plot(s.timestamp, s.ref.a__mDs2, 'DisplayName', 'accel [ref]')

subplot(3,1,2)
hold on; grid on;
plot(s.timestamp, s.mes.v__mDs, 'DisplayName', 'vel [meas]')
plot(s.timestamp, x_est(1:end-1,2), 'DisplayName', 'vel [est]');
plot(s.timestamp, s.ref.v__mDs, 'DisplayName', 'vel [ref]')
plot(s.timestamp, cumsum(s.imu.a__mDs2 .* s.dt), 'DisplayName', 'vel [int]')

subplot(3,1,3)
hold on; grid on;
plot(s.timestamp, s.ref.x__m, 'DisplayName', 'pos [ref]')
plot(s.timestamp, x_est(1:end-1,1), 'DisplayName', 'pos [est]');
plot(s.timestamp, cumsum(s.mes.v__mDs .* s.dt), 'DisplayName', 'pos [int]')
plot(s.timestamp, cumsum(cumsum(s.imu.a__mDs2 .* s.dt) .* s.dt), 'DisplayName', 'pos [int2]')

af.showlegend()