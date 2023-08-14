function [pobj1hz, ptar] = theMovingObstacleHorizonParameters(Td,N,ct,pobj1,vtar,vobj1)
    % This function calculates the center position of the obstacles over
    % the horizon. It uses the known 'current' position of the obstacle and
    % velocity (assumed to be constant).
    %
    % Inputs:
    %   Td: the discretization time, s
    %   N: the horizon length
    %   ct:    the current position of the target; (2 x 1)
    %   pobjx: the current position of obstacle x; (2 x 1)
    %   vtar:  the current velocity of the target; (2 x 1)
    %   vobjx: the current velocity of obstacle x; (2 x 1)
    % Outputs:
    %   ptar: the position of the target at the current step and the
    %            next N steps; (2 x (N+1))
    %   pobjxhz: the position of obstacle x at the current step and the
    %            next N steps; (2 x (N+1))
    
   
    % Create the multiplying array in terms of Td
    % Tdarray = 0:Td:N*Td; % 1 x (N+1)
    Tdarray =  linspace(0,N,N+1)*Td;
    
    % Create arrays of the positions matching the output size, (2 x (N+1))
    ctarray = repmat(ct,1,N+1);
    pobj1array = repmat(pobj1,1,N+1);
    
    % Compute the new target/obstacle positions over the horizon
    ptar = ctarray + vtar * Tdarray;  % (2 x (N+1))
    pobj1hz = pobj1array + vobj1 * Tdarray;  % (2 x (N+1))   
    
end






