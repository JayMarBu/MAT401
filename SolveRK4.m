function [x y z, t] = SolveRK4( M, r, h, w, t_max, dt)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % Define Principle moments of Inertia (diagonal values of inertia tensor)
    I1 = (3/20)*M*(r.^2 + (1/4)*h.^2);
    I2 = (3/20)*M*(r.^2 + (1/4)*h.^2);
    I3 = (3/20)*M*(2*r.^2);

    % Obtain values for gamma
    g1 = (I3 - I2)/I1;
    g2 = (I1 - I3)/I2;
    g3 = (I2 - I1)/I3;

    % calculate number of steps to run
    nmax = t_max / dt;

    % initialise outputs
    x = zeros(1,nmax);
    y = zeros(1,nmax);
    z = zeros(1,nmax);

    t = zeros(1,nmax);

    % set starting angular velocity
    x(1) = w(1);
    y(1) = w(2);
    z(1) = w(3);

    t(1) = 0;


    for n = 1: nmax-1

        % h = dt
        % x(n) = n
        % y(n) = output value at index n

        % k1 = hf(x(n), y(n)) 
        kx1 = -dt*g1*y(n)*z(n);
        ky1 = -dt*g2*x(n)*z(n);
        kz1 = -dt*g3*x(n)*y(n);

        % k2 = hf(x(n) + h/2, y(n) + k1/2)
        kx2 = -dt*g1*(y(n) + (ky1 * 0.5))* (z(n) + (kz1 * 0.5));
        ky2 = -dt*g2*(x(n) + (kx1 * 0.5))* (z(n) + (kz1 * 0.5));
        kz2 = -dt*g3*(x(n) + (kx1 * 0.5))* (y(n) + (ky1 * 0.5));

        % k3 = hf(x(n) + h/2, y(n) + k1/2)
        kx3 = -dt*g1*(y(n) + (ky2 * 0.5))* (z(n) + (kz2 * 0.5));
        ky3 = -dt*g2*(x(n) + (kx2 * 0.5))* (z(n) + (kz2 * 0.5));
        kz3 = -dt*g3*(x(n) + (kx2 * 0.5))* (y(n) + (ky2 * 0.5));

        % k4 = hf(x(n) + h, y(n) + k3)
        kx4 = -dt*g1*(y(n) + ky3)*(z(n) + kz3);
        ky4 = -dt*g2*(x(n) + kx3)*(z(n) + kz3);
        kz4 = -dt*g3*(x(n) + kx3)*(y(n) + ky3);

        % y(n+1) = y(n) + k1/6 + k2/3 + k3/3 + k4/6
        x(n+1) = x(n) + kx1/6 + kx2/3 + kx3/3 + kx4/6; 
        y(n+1) = y(n) + ky1/6 + ky2/3 + ky3/3 + ky4/6; 
        z(n+1) = z(n) + kz1/6 + kz2/3 + kz3/3 + kz4/6; 

        t(n+1) = t(n) + dt;
    end

end

