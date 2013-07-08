load('data.mat')

xf0 = [10;-1];

f = @(x,u,w) x + ([x(2); u])*dt;
h = @(x,u,w) x(2)/x(1) + w;
[xf, P] = ukf_sqrt(optic_flow', xf0, f, h, Q, R, control');

subplot(1,3,1);
plot(time, position, 'blue', time, xf(1,:), 'red');
subplot(1,3,2);
plot(time, velocity, 'blue', time, xf(2,:), 'red');
subplot(1,3,3);

l = size(P);
Pp = reshape(P(1,1,:),l(3),1);

error = position - xf(1,:)';
plot(time, 3*sqrt(Pp), 'black', time, -3*sqrt(Pp), 'black', time, error, 'red');

savefile = 'ukf_matlab_output.mat'
save(savefile, 'time', 'xf');
