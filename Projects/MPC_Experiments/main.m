%----------------- main --------------------%

% This file will run the solver for the problem of getting the chaser to
% the target in an optimized fashion. 
%
% It outputs many arguments of interest, particularly:
%   sol: 1 x (9k+6) vector of state and control variables
%   soliter: (9N+6) x k vector of state and control variable solutions over
%              k iterations
%   solthetat: (k+1) x 1 vector of target spacecraft theta angles
%   solpobj: (k+1) x 4 vector of the x,y positions of the 2 obstacles\
%   solptar: (k+1) x 2 vector of the x,y positions of the target
%   solmx: (k+1) x 1 vector of the entry cone hyperplane slopes
%   solJ: 1 x k vector of the cost for each iteration
%   solCPUtime: 1 x k vector of the solver time each iteration


% ================================
% ========== VARIABLES ===========
% ================================

% Key variables of the optimization problem
Td = 2;  % s
N = 15;  % number of horizons

% Run the file that defines some generic problem variables
theInitializer;

% Initial conditions of the chaser
cc = [3.5; 3.5]; % chaser center (initial)
thetac0 = pi; % rad, initial chaser angle
xi = [cc(1),cc(2),thetac0,0,0,0].'; % chaser current state, this will be overwritten in loop

% Initial conditions of the target
ct = [0.3; 0.3]; % target center
vtar = [0.01;0.01]; % target COM velocity, m/s
thetat0 = pi/4; % rad, initial target angle
thetat_dot = 2.0 * d2r; % rad/s, rotation rate of target spacecraft

% Current docking condition, the state of the COM of the chaser spacecraft  
cdock = [ct(1)+rdock*cos(thetat0),...
         ct(2)+rdock*sin(thetat0),...
         thetat0+pi,...
         -rdock*thetat_dot*sin(thetat0),...
         rdock*thetat_dot*cos(thetat0),...
         thetat_dot].';     

% Position and velocity of obstacles
testcase = 'A'; % 'A' or 'B' or 'C' ...

vobj1 = [0.0;0.0];             % velocity of the moving obstacle 1 (0.01,0.01)
vobj2 = [0.0;0.0];            % velocity of the moving obstacle 2 (-0.02,0.02)
if testcase == 'A'
    pobj1 = [2.4; 2.5];          % position of obstacle 1, Test Case A
    pobj2 = [2.5; 1.5];          % position of obstacle 2, Test Case A
elseif testcase == 'B'
    pobj1 = [2.0; 2.5];          % position of obstacle 1, Test Case B
    pobj2 = [3.0; 2.5];          % position of obstacle 2, Test Case B
elseif testcase == 'C'
    pobj1 = [2.5; 2.5];          % position of obstacle 1, Test Case C
    pobj2 = [1.5; 2.0];          % position of obstacle 2, Test Case C
elseif testcase == 'D'
    pobj1 = [1.0; 2.5];          % position of obstacle 1, Test Case D
    pobj2 = [3.0; 2.0];          % position of obstacle 2, Test Case D
elseif testcase == 'E'
    pobj1 = [0.0; 0.0];          % position of obstacle 1, Test Case E
    pobj2 = [0.0; 0.0];          % position of obstacle 2, Test Case E
elseif testcase == 'F'
    pobj1 = [4.0; 4.0];          % position of obstacle 1, Test Case F
    pobj2 = [4.0; 4.0];          % position of obstacle 2, Test Case F
elseif testcase == 'G'
    pobj1 = [3.0; 1.0];          % position of obstacle 1, Test Case G
    pobj2 = [2.0; 0.5];          % position of obstacle 2, Test Case G
end

% Filename for output files, animations, figures, etc.
PRENAME = strcat('12-Park2017-TestCase',testcase,'-Linearized-');

% Setting the upper and lower limits of the constraint equations
% There are:
%     State equation constraints (xdot = Ax+Bu): n*N
%     Thrust constraints (umax): 2*N
%     Thrust constraints (taumax): N
%     Obstacle avoidance w/ 2 moving obstacles: 2*N
%     Obstacle avoidance w/ rotating target: N
[n,~] = size(xi);
con0 = [n*N, 2*N, N, 2*N, N]; % The initial number of constraints
n_constraints0 = sum(con0);

[cl,cu] = theConstraintLimits(con0,umax,taumax);

% Initial Starting Guess (Required for Nonlinear Problems)
guess = [-0.15, -0.15, 0].';
x0 = theInitialGuesser(xi,guess,N,Td,m,Iz);

% ================================
% ===== SOLVING THE PROBLEM ======
% ================================

% Setup Options
opts = optiset('solver','ipopt','display','final');

% Initializing empty arrays to save solution data into
sol = xi.';
soliter = [];
solpobj = [];
solptar = [];
solm1 = [];
solm2 = [];
solthetat = [];
solpdock = [];
solCPUtime = [];
solJ = [];

% Values needed for the loop. The tolerance sets whether or not the chaser
% has 'docked'. The counting variable is used for 'saving data' purposes.
tol = 0.0001; % m^2
count = 1;

% This loop runs until the current position of the chaser matches the desired
% position of the COM as set by cdock (to within the set tolerance)
% COMMENT: should this loop also include the angular difference between
% chaser and target?
isdocked = 0;
while ~isdocked
    
    % The moving target and obstacle positions for the current sampling instant, k
    [pobj1hz, pobj2hz, ptar] = theMovingObstacleHorizonParameters(Td,N,ct,pobj1,pobj2,vtar,vobj1,vobj2);
    
    % The rotating target parameters for the current sampling instant, k
    [thetat,m1,m2,nc1,nc2,pdock,xt,cdock] = theRotatingTargetHorizonParameters(Td,N,ptar,rt,rhold,rdock,thetat_dot,thetah,thetat0);
    
    % Reshape the xt variable for use in theObjectiveFunction, ignore the
    % current state and reform into [03 xt(2) 03 xt(2) 03 ... 03 xt(N+1)]
    xt9 = zeros(9,N);
    xt9(4:9,:) = xt(:,2:N+1);
    xt9 = reshape(xt9, [9*N,1]);
    
    % Define arrays for the hyperplane normal arrays to be used for the entry cone constraints
    nc1m = zeros(N,2*N);
    nc2m = zeros(N,2*N);
    
    % The reference trajectory of the chaser over the horizon of instant k
    pxref = theReferenceTrajectory(x0,N); % 2xN
        
    % The points ro on the KOZ boundaries of the obstacles over the horizon
    ro_1 = theObstacleLinearizer(pxref,pobj1hz(:,2:N+1),KOZ_actual);
    ro_2 = theObstacleLinearizer(pxref,pobj2hz(:,2:N+1),KOZ_actual);
    ro_t = theObstacleLinearizer(pxref,ptar(:,2:N+1),rhold);
        
    % Determine if the entry cone constraints should be on
    if (xi(1:2) - cdock(1:2,1)).' * (xi(1:2) - cdock(1:2,1)) < rcone*rcone % if chaser is within rcone range
        withEntryCone = 1; % turn on
        
        % Add entry cone constraint limits
        if length(cu) == n_constraints0
            cl(end+1:end+2*N) = -inf;
            cu(end+1:end+2*N) = 0;
        end
        
        % Create Nx2N matrix of the entry cone normal to each hyperplane
        % Should be in the form
        % ncm = [nx1 ny1  0   0   0   0  ... 0   0   0
        %         0   0  nx2 ny2  0   0  ... 0   0   0
        %                           ...
        %         0   0   0   0   0   0  ... 0  nxN nyN ]
        for i=1:N
            nc1m(i,2*i-1) = nc1(i+1,1);
            nc1m(i,2*i) = nc1(i+1,2);
            
            nc2m(i,2*i-1) = nc2(i+1,1);
            nc2m(i,2*i) = nc2(i+1,2);
        end
        
    else
        withEntryCone = 0; % turn off
        
        % Remove entry cone constraint limits
        if length(cu) > n_constraints0
            cl(end-2*N+1:end) = [];
            cu(end-2*N+1:end) = [];
        end
        
        % FOR PLOTTING PURPOSES
        m1 = NaN;
        m2 = NaN;
    end

    % Define the b vector for Cx = b dynamics constraint
    bvec = zeros(6*N,1);
    bvec(1:6) = Ad*xi;
    
    % Objective Function
    ofun = @(x)theObjectiveFunctionCompact(x,xt9,H);
    
    % Nonlinear Constraints (cl <= nlcon(x) <= cu)
    nlcon = @(x) theNonlinearConstraintsMovingObstacleCompact(x,N,bvec,...
        Cmat,con0,pobj1hz,pobj2hz,nc1m,nc2m,pdock,withEntryCone,ro_1,ro_2,...
        ro_t,ptar,KOZ_actual,rhold);
        
    % Computing the gradient and jacobian -- Complex Step Differentiation
    grad = @(x) cstepJac(ofun,x);
    jac = @(x) cstepJac(nlcon,x);
    
    % Build OPTI Problem
    Opt = opti('fun',ofun,'grad',grad,'nl',nlcon,cl,cu,'jac',jac,'x0',x0,'options',opts);

    % Solve!
    tic
    [x,fval,ef,info] = solve(Opt,x0);
    solCPUtime(end+1) = toc;

    % Save data
    sol(end+1:end+9) = x(1:9);
    soliter(1:n,count) = xi;
    soliter(n+1:n+9*N,count) = x;
    solpobj(count,1) = pobj1hz(1,1); 
    solpobj(count,2) = pobj1hz(2,1);
    solpobj(count,3) = pobj2hz(1,1);
    solpobj(count,4) = pobj2hz(2,1);
    solptar(count,1) = ptar(1,1);
    solptar(count,2) = ptar(2,1);
    solm1(count,1) = m1(1);
    solm2(count,1) = m2(1);
    solthetat(count,1) = thetat(1);
    solpdock(count,1:2) = pdock(:,1);
    solJ(end+1) = theObjectiveFunctionCompact(x,xt9,H);
    count = count + 1;
    
    % Determine if the chaser is docked using the 'current' positions of
    % the chaser and docking condition; done using only position
    if (xi(1:2) - cdock(1:2,1)).' * (xi(1:2) - cdock(1:2,1)) < tol
        isdocked = 1;
    end
    
    % Set the next guess as the solution from the previous iteration
    xi = x(4:9); % the next 'current' state of the chaser spacecraft
    x0 = x;
    x0(1:9) = [];
    x0(end+1:end+3) = [0,0,0].';
    
    % THIS IS NEW (16 FEB 2022, NOT IN ARCHIVED FILES)
    x0(end+1:end+n) = Ad*x0(4+9*(N-2):9+9*(N-2)) + Bd*x0(1+9*(N-1):3+9*(N-1));
    
    % Update rhold if:
    %  - the position of the chaser matches the holding position, within tolerance
    %  - the attitude of the chaser matches the docking condition, within tolerance
    %  - the holding radius is larger than its minimum allowable size, rdock
    if abs((xi(1:2) - xt(1:2,1)).' * (xi(1:2) - xt(1:2,1))) < 0.5 % if position error is below tolerance
        if abs(xi(3) - cdock(3,1)) < 10*d2r % if attitude error is below tolerance
            if rhold > rdock
                rhold = rhold * gamma; % decrease the holding radius
                if rhold < rdock
                    rhold = rdock; % set to the minimum holding radius, rdock
                end
            end
        end
    end
    
    % Set the next 'current' target and obstacle positions
    ct = ptar(:,2);
    pobj1 = pobj1hz(:,2); 
    pobj2 = pobj2hz(:,2);
        
    % Set the new target theta value for the next iteration
    thetat0 = thetat(2);
    
end

% For last iteration ...
solpobj(count,1) = pobj1hz(1,2);
solpobj(count,2) = pobj1hz(2,2);
solpobj(count,3) = pobj2hz(1,2);
solpobj(count,4) = pobj2hz(2,2);
solptar(count,1) = ptar(1,2);
solptar(count,2) = ptar(2,2);
solm1(count,1) = m1(2); % NaN; %
solm2(count,1) =  m2(2); % NaN; %
solthetat(count,1) = thetat(2);
solpdock(count,1:2) = pdock(:,2);

% Save the data in the workspace to a .mat file
save(strcat(PRENAME,'dat.mat'))

