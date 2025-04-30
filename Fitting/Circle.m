%% Setup
rng(2603);

r = 10;  % radius
N = 100; % number samples
center = [2 0];

sigma = diag([0.01, 0.05]);
noise = randn(N,2) * sigma;

angle = deg2rad(randn(N,1)*30-15);

x = r.*cos(angle) + center(1) + noise(:,1);
y = r.*sin(angle) + center(2) + noise(:,2);

%% Fitting
circle_residual = @(Z, THETA) sqrt((Z(:,1) - THETA(1)).^2 + (Z(:,2) - THETA(2)).^2) - THETA(3);

% Theta0 = [0; 0; mean(sqrt(x.^2 + y.^2))];
Theta0 = [0; 0; 4];

% p_est = lsqnonlin(circle_residual, p0, [], [], opts);

Jx = @(Z, X)[(-Z(:,1)-X(1)) ./ vecnorm([Z(:,1)-X(1), Z(:,2)-X(2)],2,2), ...
             (-Z(:,2)-X(2)) ./ vecnorm([Z(:,1)-X(1), Z(:,2)-X(2)],2,2), ...
              -X(3).*ones(size(Z,1),1)];
Z = [x,y];

[Theta, res] = estimation.LevenbergMarquardt(Jx, circle_residual, Z, Theta0);
circle_residual(Z, Theta);

%% Plottnig
af = afigure;
subplot(2,1,1);
hold on; grid on; axis equal;
plot(x,y, '.', 'DisplayName', 'Data');
plot(center(1), center(2), 'x', 'Color', Color.RED, 'DisplayName', 'Center');
plot(center(1) + r.*sin(-pi:0.01:pi), center(2) + r.*cos(-pi:0.01:pi), '--', 'DisplayName', 'Perfect Circle');
plot(Theta(1), Theta(2), 'o', 'DisplayName', 'Center (est)');
plot(Theta(1) + Theta(3).*sin(-pi:0.01:pi), Theta(2) + Theta(3).*cos(-pi:0.01:pi), 'DisplayName', 'Circle (est)');

subplot(2,1,2);
hold on; grid on;
plot(res);

af.showlegend();


