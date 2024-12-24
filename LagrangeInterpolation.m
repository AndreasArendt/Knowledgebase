clear; clf;

% given points
x = [0 1 2 3 4 5 6 7];
y = [0 -5 -10 2 4  9 -4 4];

p = [];
X = -10:0.1:10;
for x_eval = X
    p(end+1) = lagrange_interp(x, y, x_eval);
end

figure;
hold on; grid on;
plot(x,y, 'x')
plot(X, p);

