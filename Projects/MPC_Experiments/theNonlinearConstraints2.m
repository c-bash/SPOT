function nl = theNonlinearConstraints2(x,N,bvec,Cmat,con0,pobj1hz,pobj2hz,nc1m,nc2m,pdock,withEC,ro_1,ro_2,ro_t,ptar,KOZobj,rhold)
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
    %
    % This function defines the constraint equations for the problem,
    % including:
    %     State equation constraints (xdot = Ax+Bu): n*N
    %     Thrust constraints (umax): 2*N
    %     Thrust constraints (taumax): N
    %     Obstacle avoidance w/ 2 moving obstacles: 2*N
    %     Obstacle avoidance w/ rotating target: N
    %     [POSITION-DEPENDENT] Entry cone constraints: 2*N
    %
    % Inputs:
    %   x: the optimization state and control vector; (9N x 1)
    %   N: the horizon length
    %   bvec: in regards to the Cx = b dynamics constraint; (6N x 1)
    %   Cmat: in regards to the Cx = b dynamics constraint; (6N x 9N)
    %   con0: an array of the number of constraints in each of the
    %         following categories (in order): state equations, control 
    %         forces, control torques, obstacle constraints
    %   pobjxhz/ptar: the position of obstacle x/target; (2 x N+1)
    %   ncXm: the array for the entry cone hyperplane constraint; (N x 2N)
    %   pdock: the position of the docking apex; (2 x N+1)
    %   withEC: Boolean of whether entry cone constraint is in effect
    %   ro_X: expansion points on the KOZ of the obstacle; (2 x N)
    %   KOZobj: the KOZ radius of the obstacles
    %   rhold: the holding radius around the target
    % Outputs:
    %   nl: the constraint equations
    
    if withEC
        nl = zeros(sum(con0)+2*N,1);
    else
        nl = zeros(sum(con0),1);
    end
        
    % Constraints for the dynamics of the problem, dotx-Ax-Bu  
    nl(1:con0(1)) = Cmat * x - bvec;
        
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
    nl(9*N+1:10*N) = ones(N,1) - 1/(KOZobj*KOZobj) * diag((ro_1-pobj1hz(:,2:N+1)).' * (2*px - ro_1 - pobj1hz(:,2:N+1)));
    nl(10*N+1:11*N) = ones(N,1) - 1/(KOZobj*KOZobj) * diag((ro_2-pobj2hz(:,2:N+1)).' * (2*px - ro_2 - pobj2hz(:,2:N+1)));
        
    % Constraints for the collision avoidance with target, rhold
    nl(11*N+1:12*N) = ones(N,1) - 1/(rhold*rhold) * diag((ro_t-ptar(:,2:N+1)).' * (2*px - ro_t - ptar(:,2:N+1)));

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