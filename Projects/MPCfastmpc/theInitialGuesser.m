function initialguess = theInitialGuesser(xi,guess,N,Td,m,Iz)
    % Make an initial guess based on ODE propagation using RK45. This is
    % used in main. We take the initial condition and the initial guess,
    % and propagate this for one iteration-equivalent step. Then, the
    % control force is set to zero and the remainder of the horizon is
    % propagated such that the chaser is just coasting with no control effort.
    % 
    % Inputs:
    %   xi: the initial condition (6 x 1)
    %   guess: the initial guess after the initial condition (3 x 1)
    %   N: the horizon length
    %   Td: the discretization time, s
    %   m: the chaser mass, kg
    %   Iz: the chaser moment of inertia, kg m^2
    % Outputs:
    %   initialguess: the starting initial guess, x0, which is (9*N x 1)

    % Set the timestep for the propagator
    stepsize = Td/10; % arbitrary timestep, doesn't matter since we take the final value at Td
    tspan = 0:stepsize:Td;
    
    % Set the initial condition
    yi = xi;
    
    % Initialize the solution array. The initial control force is non-zero, 
    % and all control forces afterwards are zero.
    initialguess = zeros(N,9);
    initialguess(1,1:3) = guess;

    % Propagate the dynamics over the horizon, N
    for i = 1:N
        % Set the control forces
        u = initialguess(i,1:3);
        
        % Calculate the state after this control force
        [~,y] = ode45(@(t,y) theStateEquations(t,y,u,m,Iz),tspan,yi);
        
        % Set the [x,y,theta,dx,dy,omega] values in the solution as the 
        % final propagated state, corresponding to the solution at Td. This
        % is for a single step in the horizon.
        initialguess(i,4:9) = y(end,:);
        
        % Set the next "initial condition" for the next step in the horizon
        yi = y(end,:);
    end
    
    % Reshape the solution
    initialguess = reshape(initialguess.',[9*N,1]);
    
end