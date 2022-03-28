% MAT401 Assignment

disp('Mat401 Assignment');

% Task 1:
% Solve Eulers equations using fourth order runge-kutta on a cone with the
% following properties:
% - M = 10kg
% - r = 1m
% - h = 4m
% - w = (3,1,2) rads^-1
% - t 0 -> 20 s

% set values
M = 10;
r = 1;
h = 4;
w_init = [3,1,2];

t_max = 20;
dt = 0.001;

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

% Task 5:
% plot the complete motion in the x-y plane of a point p on the cone
% - p = (0, 3r/4, 0)

p = [0, (3/4)*r, 0];

[x3, y3, z3] = ComputeTrajectoryOfPoint(p, x,y,z, x2,y2,z2, t_max, dt);

% displacment time graph
plot(t2, x3, 'r');
hold on
plot(t2, y3, 'b');
grid on
hold off
title('Task 5')
xlabel('time (s)')
ylabel('displacement (m)')
legend('x','y')
print -dpdf -r500 -painters task_5_displacment_over_time.pdf

% x y graph
plot(x3, y3, 'r');
grid on
hold off
title('Task 5')
xlabel('x (m)')
ylabel('y (m)')
print -dpdf -r500 -painters task_5_x_by_y.pdf
