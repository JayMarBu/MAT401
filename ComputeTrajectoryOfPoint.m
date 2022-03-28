function [ x y z ] = ComputeTrajectoryOfPoint( p, wx, wy, wz, dx, dy, dz, t_max, dt )
% - p = point on cone
% - [wx(t), wy(t), wz(t)] = angular velocity at time t
% - [dx(t), dy(t), dz(t)] = position of the centre of mass at time t
% - t_max = max time step
% - dt = time in seconds of one step

    nmax = t_max / dt;
    
    % initialise outputs
    x = zeros(1,nmax);
    y = zeros(1,nmax);
    z = zeros(1,nmax);
    
    rx = zeros(1,nmax);
    ry = zeros(1,nmax);
    rz = zeros(1,nmax);

    % set starting position
    x(1) = p(1);
    y(1) = p(2);
    z(1) = p(3);
    
    rx(1) = p(1);
    ry(1) = p(2);
    rz(1) = p(3);
    
    for n = 1:nmax-1;
        % build rotation matrix
        wPrime = [wx(n), wy(n), wz(n)] / sqrt(wx(n).^2 + wy(n).^2 + wz(n).^2);

        a = wPrime(1);
        b = wPrime(2);
        c = wPrime(3);

        theta = sqrt(wx(n).^2 + wy(n).^2 + wz(n).^2)*dt;

        cT = cos(theta);
        sT = sin(theta);

        rotMat = [ a.^2*(1-cT)+cT,  a*b*(1-cT)-c*sT,    a*c*(1-cT)+b*sT ;
                   a*b*(1-cT)+c*sT, b.^2*(1-cT)+cT,     b*c*(1-cT)-a*sT ;
                   a*c*(1-cT)-b*sT, b*c*(1-cT)+a*sT,    c.^2*(1-cT)+cT ];
        
        
       % transform point using the matrix
       rn = [rx(n), ry(n), rz(n)];
       rn = rn(:);
       rnplus1 = rotMat*rn;
       
       rx(n+1) = rnplus1(1);
       ry(n+1) = rnplus1(2);
       rz(n+1) = rnplus1(3);
       
       % calculate global position at t = n+1
       x(n+1) = dx(n+1) + rx(n+1);
       y(n+1) = dy(n+1) + ry(n+1);
       z(n+1) = dz(n+1) + rz(n+1);   
    end      
end

