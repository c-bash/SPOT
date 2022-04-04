function nl = theNonlinearConstraints0(x,N,Ad,Bd,xi,Sj,pobj1hz,pobj2hz,ct,nc1,nc2,pdock,withEC,rhold)
    % The x vector is to be optimized. If there are 3 degrees of freedom (i.e.
    % we are considering x y theta, and their velocities, xdot ydot thetadot,
    % and their control forces, ux uy tau), then the length of x is 9*N
    % and has the form:
    %
    % x = [ux0 uy0 tau0 x1 y1 theta1 xdot1 ydot1 thetadot1
    %   ...ux1 uy1 tau1 x2 y2 theta2 xdot2 ydot2 thetadot2 ...
    %   ...ux2 uy2 tau2 x3 y3 theta3 xdot3 ydot3 thetadot3...
    %   ...ux3 uy3 tau3 x4 y4 theta4 xdot4 ydot4 thetadot4...
    %   ... 
    %   ...uxN-1 uyN-1 tauN-1 xN yN thetaN xdotN ydotN thetadotN]

    [n,~] = size(xi); % the length of the initial conditions
    Nu = 3; % the number of control forces
    Nxu = n + Nu; % the combined size of the x and u parameters for a single horizon step
    
    nl = zeros(n,1);
    
    % Constraints for the dynamics of the problem, dotx-Ax-Bu
    nl(1:n) = x(4:9)-Ad*xi-Bd*x(1:3);
    for i = 0:N-2    
        nl(end+1:end+n) = x((i+1)*9+4:(i+2)*9)-Ad*x(i*9+4:(i+1)*9)-Bd*x(i*9+10:i*9+12);
    end
    
    % Constraints for the control forces, u, i.e. Fx and Fy
    Ni2 = (0:N-1).';
    nl(end+1:end+N) = x(1+9*Ni2); 
    nl(end+1:end+N) = x(2+9*Ni2);  
    
    % Constraints for the control forces, tau
    nl(end+1:end+N) = x(3+9*Ni2);
    
    % Constraints for the obstacle avoidance, obj1 and obj2, moving over
    % the horizon
    for i = 0:N-1
        pobj1 = pobj1hz(:,i+2);
        %pobj1 = pobj1hz; % for constant object position
        nl(end+1) = 1 - (x(i*9+4:i*9+5) - pobj1).' * Sj * (x(i*9+4:i*9+5) - pobj1);
    end
    for i = 0:N-1
        pobj2 = pobj2hz(:,i+2);
        %pobj2 = pobj2hz; % for constant object position
        nl(end+1) = 1 - (x(i*9+4:i*9+5) - pobj2).' * Sj * (x(i*9+4:i*9+5) - pobj2);
    end
    
    % Constraints for the collision avoidance with target, rhold
    for i = 0:N-1
        nl(end+1) = rhold*rhold - (x(i*9+4:i*9+5) - ct).' * (x(i*9+4:i*9+5) - ct);
    end
    
    % Constraints for entry cone
    if withEC
        for i = 0:N-1
           %nc1 = [m1(i+2), -1];
           nc1i = nc1(i+2,:);
           pdocki = pdock(:,i+2);
           %nc1 = [m1, -1];
           %pdocki = pdock;
           nl(end+1) = -nc1i * x(i*9+4:i*9+5) + nc1i * pdocki;
        end
        for i = 0:N-1
           %nc2 = [m2(i+2), -1];
           nc2i = nc2(i+2,:);
           pdocki = pdock(:,i+2);
           %nc2 = [m2, -1];
           %pdocki = pdock;
           nl(end+1) = nc2i * x(i*9+4:i*9+5) - nc2i * pdocki;
        end
    end


end