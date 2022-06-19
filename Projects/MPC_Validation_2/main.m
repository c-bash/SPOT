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
% ========== ALGORITHM ===========
% ================================

% Choose which algorithm configuration to use to solve the problem
%    Linearized Compact MPC (Quadprog): 5
%    MPC (ISNM): 6
%    Fast MPC (ISNM): 7
ALGORITHM = 7;
isCONTROL = 1;

LASTCOUNT = 44;

if ALGORITHM == 5
    % MPC
    % Quadprog
    % Objective function = Split into sections, 2(0.5*x.'*H*x - xt.'*H*x) + xt.'*H*xt
    % Linear constraints = big matrix multiply, split into = and <= constraints
    algorithmname = 'Quadprog';
elseif ALGORITHM == 6
    % MPC
    % Infeasible Newton Step
    algorithmname = 'ISNM';
elseif ALGORITHM == 7
    % Fast MPC
    % Infeasible Newton Step
    algorithmname = 'FastMPC';
end

% ================================
% ========== VARIABLES ===========
% ================================

% Key variables of the optimization problem
Td = 2;  % s
N = 15;  % number of horizons

% Run the file that defines some generic problem variables
theInitializer;

% Position and velocity of obstacles
testcase = 'C'; % 'A' or 'B' or 'C' ...

if testcase == 'A'
    cc = [3.0; 2.0]; % chaser center (initial)
    thetac0 = pi; % rad, initial chaser angle

    ct = [0.5; 0.5]; % target center
    vtar = [0.01;0.01]; % target COM velocity, m/s
    thetat0 = -pi/4; % rad, initial target angle
    thetat_dot = 1.0 * d2r; % rad/s, rotation rate of target spacecraft
elseif testcase == 'B'
    cc = [3.3; 1.5]; % chaser center (initial)
    thetac0 = pi; % rad, initial chaser angle

    ct = [2.0; 1.9]; % target center
    vtar = [-0.01;-0.005]; % target COM velocity, m/s
    thetat0 = -pi/2; % rad, initial target angle
    thetat_dot = -1.25 * d2r; % rad/s, rotation rate of target spacecraft

elseif testcase == 'C'
    cc = [1.25; 1.25]; % chaser center (initial)
    thetac0 = pi; % rad, initial chaser angle

    ct = [2.5; 1.25]; % target center
    vtar = [-0.008;0.008]; % target COM velocity, m/s
    thetat0 = -pi/2; % rad, initial target angle
    thetat_dot = -2.0 * d2r; % rad/s, rotation rate of target spacecraft
end

if ~isCONTROL
    if testcase == 'A'
        pobj1 = [1.5; 1.25];          % position of obstacle 1, Test Case A
        vobj1 = [0.01;0.01];          % velocity of the moving obstacle 1
    elseif testcase == 'B'
        pobj1 = [2.5; 1.2];          % position of obstacle 1, Test Case B
        vobj1 = [-0.01;-0.01];        % velocity of the moving obstacle 1
    elseif testcase == 'C'
        pobj1 = [2.0; 0.75];          % position of obstacle 1, Test Case C
        vobj1 = [-0.005;0.005];        % velocity of the moving obstacle 1
    end
else
    pobj1 = [10.; 10.];          
    vobj1 = [0.00;0.00];
end

% Initial conditions of the chaser
xi = [cc(1),cc(2),thetac0,0,0,0].'; % chaser current state, this will be overwritten in loop

% Current docking condition, the state of the COM of the chaser spacecraft  
cdock = [ct(1)+rdock*cos(thetat0+cdockoa),...
         ct(2)+rdock*sin(thetat0+cdockoa),...
         thetat0+pi,...
         vtar(1)-rdock*thetat_dot*sin(thetat0+cdockoa),...
         vtar(2)+rdock*thetat_dot*cos(thetat0+cdockoa),...
         thetat_dot].';   


% Filename for output files, animations, figures, etc.
PRENAME = strcat('2-TestCase',testcase,'-',algorithmname,'-');

if isCONTROL
    PRENAME = strcat(PRENAME, 'CTRL-');
end

% Setting the upper and lower limits of the constraint equations
% There are:
%     State equation constraints (xdot = Ax+Bu): n*N
%     Thrust constraints (umax): 2*N
%     Thrust constraints (taumax): N
%     Obstacle avoidance w/ 1 moving obstacles: N
%     Obstacle avoidance w/ rotating target: N
[n,~] = size(xi);
con0 = [n*N, 2*N, N, N, N]; % The initial number of constraints
n_constraints0 = sum(con0);

% Initial Starting Guess (Required for Nonlinear Problems)
guess = [-0.14, 0, 0].';
x0 = theInitialGuesser(xi,guess,N,Td,m,Iz);
maxiterations = 10;

% ================================
% ===== SOLVING THE PROBLEM ======
% ================================

% Setup Options
if ALGORITHM == 5
    options = optimoptions('quadprog','MaxIterations',maxiterations);
elseif ALGORITHM >= 6
    nu = zeros(6*N,1);
end

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

% TS_solcost = []; % DELETE ME
% solINEQ = [];
% TS_P = [];
% TS_h = [];
% TS_zt = [];
% TS_dz = [];

% Values needed for the loop. The tolerance sets whether or not the chaser
% has 'docked'. The counting variable is used for 'saving data' purposes.
tol = 0.00015; % m^2
count = 1;

% This loop runs until the current position of the chaser matches the desired
% position of the COM as set by cdock (to within the set tolerance)
% COMMENT: should this loop also include the angular difference between
% chaser and target?
isdocked = 0;
while ~isdocked
    
    % The moving target and obstacle positions for the current sampling instant, k
    [pobj1hz, ptar] = theMovingObstacleHorizonParameters(Td,N,ct,pobj1,vtar,vobj1);
    
    % The rotating target parameters for the current sampling instant, k
    [thetat,m1,m2,nc1,nc2,pdock,xt,cdock] = theRotatingTargetHorizonParameters(Td,N,ptar,vtar,pdockr,rhold,rdock,cdockoa,pdockoa,ecoa,thetat_dot,thetah,thetat0,cdockdeltax);
    
    % Reshape the xt variable for use in theObjectiveFunction, ignore the
    % current state and reform into [03 xt(2) 03 xt(2) 03 ... 03 xt(N+1)]
    xt9 = zeros(9,N);
    xt9(4:9,:) = xt(:,2:N+1);
    xt9 = reshape(xt9, [9*N,1]);
    
    %xt9 = zeros(9*N,1);
    
%         disp(x0(1:20))
    % The reference trajectory of the chaser over the horizon of instant k
    pxref = theReferenceTrajectory(x0,N); % 2xN

    % The points ro on the KOZ boundaries of the obstacles over the horizon
    ro_1 = theObstacleLinearizer(pxref,pobj1hz(:,2:N+1),KOZ_actual);
    ro_t = theObstacleLinearizer(pxref,ptar(:,2:N+1),rhold);
    
    % Reshape pdock for use in the linearized matrix constraints
    pdock9 = zeros(9,N);
    pdock9(4:5,:) = pdock(:,2:N+1);
    pdock9 = reshape(pdock9, [9*N,1]);

    % Make the ro_X and robj9 vectors for obstacle constraints
    ro_1_9 = zeros(9,N);
    ro_1_9(4:5,:) = ro_1;
    ro_1_9 = reshape(ro_1_9,9*N,1);

    ro_t_9 = zeros(9,N);
    ro_t_9(4:5,:) = ro_t;
    ro_t_9 = reshape(ro_t_9,9*N,1);

    robj1_9 = zeros(9,N);
    robj1_9(4:5,:) = pobj1hz(:,2:N+1);
    robj1_9 = reshape(robj1_9,9*N,1);

    rtar_9 = zeros(9,N);
    rtar_9(4:5,:) = ptar(:,2:N+1);
    rtar_9 = reshape(rtar_9,9*N,1);

    % Make the matrix for obstacle constraints
    S1 = Sj;
    St = [1/(rhold*rhold), 0; 0, 1/(rhold*rhold)];

    E1 = (ro_1 - pobj1hz(:,2:N+1)).' * S1;
    Et = (ro_t - ptar(:,2:N+1)).' * St;

    D1 = zeros(N,(9*N));
    Dt = zeros(N,(9*N));
    for row = 1:N
        D1(row,4+(row-1)*9:5+(row-1)*9) = E1(row,:);
        Dt(row,4+(row-1)*9:5+(row-1)*9) = Et(row,:);
    end  
    
    % Define the b vector for Cx = b dynamics constraint
    bvec = zeros(6*N,1);
    bvec(1:6) = Ad*xi;
    
    % Define arrays for the hyperplane normal arrays to be used for the entry cone constraints

    nc1m = zeros(N,9*N);
    nc2m = zeros(N,9*N);
        
    % Determine if the entry cone constraints should be on
    if ((xi(1:2) - cdock(1:2,1)).' * (xi(1:2) - cdock(1:2,1)) < rcone*rcone) % if chaser is within rcone range
        
        temp1 = -nc1(2,:) * x0(4:5) + nc1(2,:) * pdock9(4:5);
        temp2 = nc2(2,:) * x0(4:5) - nc2(2,:) * pdock9(4:5);
        
        if temp1 <= 0 && temp2 <= 0 % only turn ec constraints on if it won't cause infeasibilities at the next step
            
        
            withEntryCone = 1; % turn on
            con = [n*N, 2*N, N, N, N, N, N];

            % Create Nx2N matrix of the entry cone normal to each hyperplane
            % Should be in the form
            % ncm = [ 0_1x3 n1_1x2 0_1x4 0_1x3  0_1x2 0_1x4  ... 0_1x3 0_1x2 0_1x4
            %         0_1x3  0_1x2 0_1x4 0_1x3 n2_1x2 0_1x4  ... 0_1x3 0_1x2 0_1x4
            %                           ...
            %         0_1x3  0_1x2 0_1x4 0_1x3  0_1x2 0_1x4  ... 0_1x3 nN+1_1x2 0_1x4
            for i=1:N
                nc1m(i,4+(i-1)*9:5+(i-1)*9) = nc1(i+1,:);
                nc2m(i,4+(i-1)*9:5+(i-1)*9) = nc2(i+1,:);
            end
        else
            
            withEntryCone = 0; % turn off
            con = con0;

            % FOR PLOTTING PURPOSES
            m1 = NaN;
            m2 = NaN;
        end
                
    else
        withEntryCone = 0; % turn off
        con = con0;
        
        % FOR PLOTTING PURPOSES
        m1 = NaN;
        m2 = NaN;
    end
    
    if ALGORITHM == 5     
        % Constraint Matrix
        if ~withEntryCone
            Wmat = [-2*D1;-2*Dt];
        else
            Wmat = [-2*D1;-2*Dt;-nc1m;nc2m];
        end
    % THIS NEEDS TO BE CLEANED UP
    elseif ALGORITHM == 6 || ALGORITHM == 7
        if withEntryCone
            Wmat = zeros(10*N,9*N);
            hvec = zeros(10*N,1);
        else
            Wmat = zeros(8*N,9*N);
            hvec = zeros(8*N,1);
        end

        for index = 0:N-1
            blockc = 9*index;
            if withEntryCone
                blockr = 10*index;
                Wmat(9+blockr,:) = -nc1m(index+1,:);
                Wmat(10+blockr,:) = nc2m(index+1,:);
            else
                blockr = 8*index;
            end
            Wmat(1+blockr:3+blockr,1+blockc:3+blockc) = eye(3);
            Wmat(4+blockr:6+blockr,1+blockc:3+blockc) = -eye(3);
            Wmat(7+blockr,:) = -2*D1(index+1,:);
            Wmat(8+blockr,:) = -2*Dt(index+1,:);
        end

        f_D1 = -ones(N,1) - D1*(ro_1_9 + robj1_9);
        f_Dt = -ones(N,1) - Dt*(ro_t_9 + rtar_9);   
        if withEntryCone % if entry cone constraints turned on
            f_ec1 = -nc1m * pdock9;
            f_ec2 = nc2m * pdock9;
        end

        for index = 0:N-1
            if withEntryCone % if entry cone constraints turned on
                blockr = 10*index;
                hvec(9+blockr) = f_ec1(index+1);
                hvec(10+blockr) = f_ec2(index+1);
            else
                blockr = 8*index;
            end

            hvec(1+blockr:3+blockr) = [umax;umax;taumax];
            hvec(4+blockr:6+blockr) = [umax;umax;taumax];
            hvec(7+blockr) = f_D1(index+1);
            hvec(8+blockr) = f_Dt(index+1);
        end

    end
    
    
    
    % Constraint upper and lower limits
    if ALGORITHM == 5
        [cl,cu,bineqvec] = theConstraintLimits4(Umat,Fmat,Wmat,con0,con,...
            D1,Dt,ro_1_9,ro_t_9,robj1_9,rtar_9,nc1m,...
            nc2m,pdock9);
    end
    
    if ALGORITHM == 6 || ALGORITHM == 7
        if any((hvec - Wmat*x0) <= 0)
            disp('Started with infeasible solution (outside).');
            infeas0 = find((hvec - Wmat*x0) <= 0);
            disp(infeas0)

            % Find a new feasible z-val for the violated constraint(s)
            for ind = infeas0.'
                if size(con) == size(con0) % no ec constraints yet
                    ind_startrange = 7 + 8*idivide(ind,int16(8+1));
                    ind_endrange = ind_startrange + 1;
                else
                    ind_startrange = 7 + 10*idivide(ind,int16(10+1));
                    ind_endrange = ind_startrange + 3;
                end
%                 disp(ind_startrange)
%                 disp(ind_endrange)
                    
                z_infeas_ind = find(Wmat(ind,:)); % optimization variable(s) to vary
    %             disp('z(infeas)=')
    %             disp(z_infeas_ind)

                newind = z_infeas_ind - 9;
                
                search_range = 1;
                search_stepnum = 1001;
                
                while (1)
                    if newind < 0
                        disp("Sampling...")

                        % make a grid of search points
                        try_x = linspace(x0(z_infeas_ind(1))-search_range,x0(z_infeas_ind(1))+search_range,search_stepnum);
                        try_y = linspace(x0(z_infeas_ind(2))-search_range,x0(z_infeas_ind(2))+search_range,search_stepnum);
                        [xx,yy] = meshgrid(try_x,try_y);

                        % calcuate the constraint at each grid point
%                         temp = zeros(search_stepnum,search_stepnum);
                        temp2 = zeros(search_stepnum*search_stepnum,3);
                        for i1=1:search_stepnum
                            for i2 = 1:search_stepnum
                                temp = hvec(ind_startrange:ind_endrange) - Wmat(ind_startrange:ind_endrange,z_infeas_ind)*[xx(i1,i2);yy(i1,i2)];

                                % take note of grid points that satisfy the constraint
                                if ~any(temp <= 0)
                                    temp2((i1-1)*search_stepnum+i2,1) = min(temp);
                                    temp2((i1-1)*search_stepnum+i2,2) = i1;
                                    temp2((i1-1)*search_stepnum+i2,3) = i2;
                                end
                            end
                        end

                        % get rid of data that does not satisfy the constraints
                        temp2( ~any(temp2,2), : ) = [];  %rows
                        
                        if isempty(temp2)
                            search_range = 3;
                            search_stepnum = search_stepnum * 2;
                            continue
                        end

                        % choose one of the grid points that satisfy the constraint
                        % currently finds the first minimum value... could be further improved
                        ind_of_min = find(temp2(:,1) == min(temp2(:,1)));
                        ind_chosen = ind_of_min(1);

                        % update the x0 with the new grid point
                        [i1i2_chosen] = temp2(ind_chosen,[2,3]);
                        i1_chosen = i1i2_chosen(1);
                        i2_chosen = i1i2_chosen(2);
                        x0(z_infeas_ind) = [xx(i1_chosen,i2_chosen);yy(i1_chosen,i2_chosen)];                     
                        
                    else
                        x0(z_infeas_ind) = x0(newind); % re-tracing previous sol'n
                    end

                    if (hvec(ind)- Wmat(ind,:)*x0) <= 0
                        if newind < 0
                            disp("Solution not in sample range.")
                        else
                            newind = newind - 9;
                        end
                    else
                        disp("Resolved.")
                        break
                        
                    end
                end
            end
        end
       
    end


    
    % Solve!
    if ALGORITHM == 5
        tic
        [x,fval] = quadprog(H,-H.'*xt9,Wmat,bineqvec,Cmat,bvec,cl,cu,x0,options);
        solCPUtime(end+1) = toc;
    elseif ALGORITHM == 6
        disp('----------------------------------------------')
        disp(string(count))
        
        tic
        [x,nu] = theFastMPCSolver2(x0,nu,xt9,H,Cmat,bvec,Wmat,hvec,1e-9,0.01,0.95,maxiterations);
        solCPUtime(end+1) = toc;
    elseif ALGORITHM == 7
        disp('----------------------------------------------')
        disp(string(count))
        tic
        [x,nu] = theFastMPCSolver3(x0,nu,xt9,H,Cmat,bvec,Wmat,hvec,1e-9,0.01,0.95,maxiterations);
        solCPUtime(end+1) = toc;
    end
   

    % Save data
    sol(end+1:end+9) = x(1:9);
    soliter(1:n,count) = xi;
    soliter(n+1:n+9*N,count) = x;
    solpobj(count,1) = pobj1hz(1,1); 
    solpobj(count,2) = pobj1hz(2,1);
%    solpobj(count,3) = pobj2hz(1,1);
%    solpobj(count,4) = pobj2hz(2,1);
    solptar(count,1) = ptar(1,1);
    solptar(count,2) = ptar(2,1);
    solm1(count,1) = m1(1);
    solm2(count,1) = m2(1);
    solthetat(count,1) = thetat(1);
    solpdock(count,1:2) = pdock(:,1);
    solJ(end+1) = theObjectiveFunction123(x,xt9,H);
    
    
%     solINEQ(1:7*N,count) = Wmat*x-hvec;
% THIS IS ACTUALLY USEFUL, UNCOMMENT LATER FOR FASTMPC
%     P_s = length(TS_P);
%     P_e = size(Wmat,1);
%     TS_P(P_s+1:P_s+P_e,1:9*N) = Wmat;
%     TS_h(1:length(hvec),count) = hvec;
%     TS_zt(1:9*N,count) = xt9;
    count = count + 1;
%     
%     % Troubleshooting data % DELETE ME
%     cstcnt_s = length(TS_solcost);
%     cstcnt_e = 1; %length(TS_cost);
%     TS_solcost(cstcnt_s+1:cstcnt_s+cstcnt_e,1) = TS_cost(:,1);
%     TS_solcost(cstcnt_s+1:cstcnt_s+cstcnt_e,2) = TS_cost(:,2);
%     TS_solcost(cstcnt_s+1:cstcnt_s+cstcnt_e,3) = TS_cost(:,3);
%     TS_solcost(cstcnt_s+1:cstcnt_s+cstcnt_e,4) = TS_cost(:,4);
%     TS_solcost(cstcnt_s+1:cstcnt_s+cstcnt_e,5) = TS_cost(:,5);
%     
%     dz_s = size(TS_dz,2);
%     dz_e = size(TS_deltaz,2);
%     TS_dz(:,dz_s+1:dz_s+dz_e) = TS_deltaz;
    

    
    
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
    %x0(end+1:end+n) = Ad*x0(4+9*(N-2):9+9*(N-2)) + Bd*x0(1+9*(N-1):3+9*(N-1));
    % TRYING TO GET RID OF INFEASIBLE START
    x0(end+1:end+n) = x0(4+9*(N-2):9+9*(N-2));
    
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
%     pobj2 = pobj2hz(:,2);
        
    % Set the new target theta value for the next iteration
    thetat0 = thetat(2);
    
    if count == LASTCOUNT
        break
    end
    
end

% For last iteration ...
solpobj(count,1) = pobj1hz(1,2);
solpobj(count,2) = pobj1hz(2,2);
% solpobj(count,3) = pobj2hz(1,2);
% solpobj(count,4) = pobj2hz(2,2);
solptar(count,1) = ptar(1,2);
solptar(count,2) = ptar(2,2);
solm1(count,1) = m1(2); % NaN; %
solm2(count,1) =  m2(2); % NaN; %
solthetat(count,1) = thetat(2);
solpdock(count,1:2) = pdock(:,2);


% Save the data in the workspace to a .mat file
save(strcat(PRENAME,'dat.mat'))

