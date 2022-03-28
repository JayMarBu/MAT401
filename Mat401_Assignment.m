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

% Task 2 :
% plot results of Task 1
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

print -dpdf -r500 -painters task_1_and_2.pdf

% Task 3:
% Solve the Semi-Implcit Euler equations of the cone  to obtrain its 
% trajctory over the same time period using the following properties:
% - u = (0,0,200)ms-1
% - a = (0,0,-9.8)ms-2
% - t 0 -> 20 s

[x2 y2 z2 vx2 vy2 vz2 t2] = SolveSemiImplicitEuler([0,0,200],[0,0,9.8],t_max,dt);

% Task 4:
% plot results of Task 3

% displacment graph
plot(t2, x2, 'r');
hold on
plot(t2, y2, 'b');
hold on
plot(t2, z2, 'g');
hold off
title('Task 3 & 4 (Displacement graph)')
xlabel('time (s)')
ylabel('displacement (m)')
legend('x','y','z')
print -dpdf -r500 -painters task_3_and_4_displacement.pdf

% Velocity graph
plot(t2, vx2, 'r');
hold on
plot(t2, vy2, 'b');
hold on
plot(t2, vz2, 'g');
hold off
title('Task 3 & 4 (Velocity graph)')
xlabel('time (s)')
ylabel('velocity (ms-1)')
legend('x','y','z')
print -dpdf -r500 -painters task_3_and_4_velocity.pdf