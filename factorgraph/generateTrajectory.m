function s = generateTrajectory(options)

arguments
    options.seed = 1208
    options.dt = 0.01;
    options.tEnd = 100;
    options.std_accel = 1;
    options.std_vel = 0.2;
end

rng(options.seed);

f = 0.05; % [Hz]

s = struct();
s.timestamp = 0:options.dt:options.tEnd-options.dt;
s.dt = options.dt;

s.ref.a__mDs2 = sin(2*pi*s.timestamp*f) + cos(2*pi*s.timestamp*3*f);
s.ref.v__mDs  = cumsum(s.ref.a__mDs2 * s.dt);
s.ref.x__m    = cumsum(0.5 .* s.ref.v__mDs .* s.dt);

epsilon_acc = options.std_accel .* randn(1, numel(s.timestamp)); % noise
s.imu.a__mDs2 = s.ref.a__mDs2 + epsilon_acc;

epsilon_mes = options.std_vel .* randn(1, numel(s.timestamp)); % noise
s.mes.v__mDs  = s.ref.v__mDs + epsilon_mes;