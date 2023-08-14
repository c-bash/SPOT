function px = theReferenceTrajectory(x0,N)
    % Extract from the initial guess the x,y positions of the trajectory in
    % order to make a 2xN matrix of position values denoting the reference
    % trajectory. This is used for the Linearized obstacle constraints.
    %
    % Inputs:
    %   x0: the initial guess, 9N x 1
    %   N: the horizon length
    % Outputs:
    %   px: the reference trajectory of the chaser; 2 x N
    
    xreshape = reshape(x0,9,N);     % 9xN
    px = xreshape(4:5,:);           % 2xN
end