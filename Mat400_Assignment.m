% MAT400 Assignment

disp('Mat400 Assignment');

% Task 1:
% Solve Eulers equations using fourth order runge-kutta on a cone with the
% following properties:
% - M = 10kg
% - r = 1m
% - h = 4m
% - ? = (3,1,2) rads^-1
% - t 0 -> 20 s

% set values
M = 10;
r = 1;
h = 4;
w_init = [3,1,2];

t_max = 20;
dt = 0.1;

% solve for angular velocity
[x y z t] = SolveRK4(M, r, h, w_init, t_max, dt);

% plot results
plot(t, x, 'r');
hold on
plot(t, y, 'b');
hold on
plot(t, z, 'g');
hold off
title('Task 1 & 2')
xlabel('time (s)')
ylabel('angular velocity (rads-1)')
legend('x','y','z')

% Task 2:
% Solve the Semi-Implcit Euler equations of the cone  to obtrain its 
% trajctory over the same time period using the following properties:
% - u = (0,0,200)ms-1
% - a = (0,0,-9.8)ms-2
% - t 0 -> 20 s

