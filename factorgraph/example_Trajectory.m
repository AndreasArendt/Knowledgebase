clear;

std_accel = 1;
std_vel = 0.01;
s = generateTrajectory('tEnd', 60, 'std_accel', std_accel, 'std_vel', std_vel);

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

% z := [X0; B*u; y]  % prior, dynamics, measurements
Bu = reshape(B * s.imu.a__mDs2(1:end-1), [], 1); % B * u - control input
z = s.mes.v__mDs;                                % y - Measurement

y = sparse(nan(nStates + numel(Bu) + numel(z), 1));
y(1:nStates)                   = X0; % initial state
y(nStates+1:nStates+numel(Bu)) = Bu;
y(nStates+numel(Bu)+1:end)     = z;

% Final system matrix
A_prior = sparse(zeros(nStates,nStates*k));
A_prior(1:nStates,1:nStates) = eye(nStates);

A_dyn = sparse(zeros((k-1)*nStates,k*nStates));
for r = 1:nStates:(k-1)*nStates
    idxR = r:r+nStates-1;
    idxC = r:r+nStates-1;
    A_dyn(idxR, idxC)         = -F;
    A_dyn(idxR, idxC+nStates) = eye(nStates);
end

A_meas = sparse(zeros(numel(z),k*nStates));
r_ctr = 1;
for r = 1:nStates:numel(z)*nStates
    A_meas(r_ctr, r:r+nStates-1) = H;
    r_ctr = r_ctr+1;
end

A = [ A_prior;
      A_dyn;
      A_meas];

% Covariance
P0 = diag([10, 1]);
Q  = (std_accel^2) * (B * B');
R  = std_vel^2;

P0P0 = diag(P0);
QQ   = repmat(diag(Q), k-1, 1);
RR   = repmat(diag(R), numel(z), 1);

W = sparse(diag([P0P0; QQ; RR]));

% solving MAP with linear least squares
[x_hat, Sigma_x, residuals] = estimation.lls(A, y, R=W);

X_est = reshape(full(x_hat), 2, [])';
Sigma = reshape(diag(Sigma_x), nStates, [])';
res = full(residuals(end-numel(z)+1:end));

%%
af = afigure;
subplot(3,1,1);
hold on; grid on;
plot(s.timestamp, s.imu.a__mDs2, 'DisplayName', 'accel [meas]')
plot(s.timestamp, s.ref.a__mDs2, 'DisplayName', 'accel [ref]')

subplot(3,1,2)
hold on; grid on;
plot(s.timestamp, s.mes.v__mDs, 'DisplayName', 'vel [meas]')
plot(s.timestamp, X_est(:,2), 'DisplayName', 'vel [est]');
plot(s.timestamp, s.ref.v__mDs, 'DisplayName', 'vel [ref]')
plot(s.timestamp, cumsum(s.imu.a__mDs2 .* s.dt), 'DisplayName', 'vel [int]')

subplot(3,1,3)
hold on; grid on;
plot(s.timestamp, s.ref.x__m, 'DisplayName', 'pos [ref]')
plot(s.timestamp, X_est(:,1), 'DisplayName', 'pos [est]');
plot(s.timestamp, cumsum(s.mes.v__mDs .* s.dt), 'DisplayName', 'pos [int]')
plot(s.timestamp, cumsum(cumsum(s.imu.a__mDs2 .* s.dt) .* s.dt), 'DisplayName', 'pos [int2]')

af.showlegend()