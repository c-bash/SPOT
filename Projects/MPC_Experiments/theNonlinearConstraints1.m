function nl = theNonlinearConstraints1(x,N,Ad,Cmat,xi,pobj1hz,pobj2hz,ct,nc1m,nc2m,pdock,withEC,KOZobj,rhold)
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
    
    if withEC
        nl = zeros(6*N+3*N+3*N+2*N,1);
    else
        nl = zeros(6*N+3*N+3*N,1);
    end
        
    % Constraints for the dynamics of the problem, dotx-Ax-Bu
    bvec = zeros(6*N,1);
    bvec(1:6) = Ad*xi;
    
    nl(1:6*N) = Cmat * x - bvec;
        
    % Constraints for the control forces, u, i.e. Fx and Fy
    Ni2 = (0:N-1).';
    nl(6*N+1:7*N) = x(1+9*Ni2); 
    nl(7*N+1:8*N) = x(2+9*Ni2);  
    
    % Constraints for the control forces, tau
    nl(8*N+1:9*N) = x(3+9*Ni2);
       
    % Make the x array a 2xN matrix of position values    
    xreshape = reshape(x,9,N);
    px = xreshape(4:5,:);   % 2xN matrix of chaser position values
    
    % Constraints for the obstacle avoidance, obj1 and obj2, moving over
    % the horizon
    zbar1 = px - pobj1hz(:,2:N+1); % 2xN
    zbar2 = px - pobj2hz(:,2:N+1); % 2xN
    
    nl(9*N+1:10*N) = ones(N,1) - 1/(KOZobj*KOZobj) * diag(zbar1.' * zbar1);
    nl(10*N+1:11*N) = ones(N,1) - 1/(KOZobj*KOZobj) * diag(zbar2.' * zbar2);
        
    % Constraints for the collision avoidance with target, rt and rhold
    ptar = ct;
    ptar = repmat(ptar,1,N); % 2xN
    
    zbartar = px - ptar; % 2xN
    
    nl(11*N+1:12*N) = ones(N,1) - 1/(rhold*rhold) * diag(zbartar.' * zbartar);
    
    % Constraints for entry cone
    if withEC
        
        pxvec = reshape(px,2*N,1); % 2Nx1 vector of the chaser position
        
        % pdock is 2xN+1, and we want to convert it to 2Nx1, remember to
        % remove the first column which coincides with the current pos
        pdockvec = reshape(pdock(:,2:N+1),2*N,1);
               
        % Entry cone constraint equations
        nl(12*N+1:13*N) = -nc1m * pxvec + nc1m * pdockvec;
        nl(13*N+1:14*N) = nc2m * pxvec - nc2m * pdockvec;
                
    end


end