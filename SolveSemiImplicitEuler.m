function [x y z vx vy vz t] = SolveSemiImplicitEuler( u, a, t_max, dt )
% - u = initial velocity
% - a = acceleration
% - t_max = max time step
% - dt = time in seconds of one step

    nmax=t_max/dt;
    
    % initialise outputs
    x = zeros(1,nmax);
    y = zeros(1,nmax);
    z = zeros(1,nmax);

    t = zeros(1,nmax);
    
    vx = zeros(1,nmax);
    vy = zeros(1,nmax);
    vz = zeros(1,nmax);

    % set starting position & velocity
    x(1) = 0;
    y(1) = 0;
    z(1) = 0;

    t(1) = 0;
    
    vx(1) = u(1);
    vy(1) = u(2);
    vz(1) = u(3);

    for n = 1:nmax-1;
        % calculate velocity
        vx(n+1) = vx(n) - dt*a(1);
        vy(n+1) = vy(n) - dt*a(2);
        vz(n+1) = vz(n) - dt*a(3);
        
        % calculte displacement
        x(n+1) = x(n) + dt*vx(n+1);
        y(n+1) = y(n) + dt*vy(n+1);
        z(n+1) = z(n) + dt*vz(n+1);
        
        t(n+1) = t(n) + dt;
    end

end

