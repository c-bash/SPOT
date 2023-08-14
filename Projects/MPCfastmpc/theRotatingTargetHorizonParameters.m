function [thetat,m1,m2,nc1,nc2,pdock,xt,cdock] = theRotatingTargetHorizonParameters(Td,N,ptar,vtar,pdockr,rhold,rdock,cdockoa,pdockoa,ecoa,thetat_dot,thetah,thetat0,cdockdeltax)
    % This function calculates angle of the target spacecraft over the
    % horizon, then uses this to calculate the new slopes for the entry
    % cone constraints, the new docking state, and the new target
    % state for the chaser to go to, xt.
    %
    % Inputs:
    %   Td: the discretization time, s
    %   N: the horizon length
    %   ptar: the position of the target over the horizon; (2 x (N+1))
    %   vtar: the velocity of the target over the horizon; (2 x 1)
    %   rt: the width of the target spacecraft, m
    %   rhold: the holding radius at the current iteration, m
    %   rc: the width of the chaser spacecraft, m
    %   wdock: the width of the docking apparatus, m
    %   thetat_dot: the rotation rate of the target spacecraft, rad/s
    %   thetah: the half-angle for the entry cone, rads
    %   thetat0: the current orientation of the target, rads
    % Outputs:
    %   thetat: the target spacecraft orientation (current and over N horizon); 1 x (N+1)
    %   mX: the slopes of the entry cone hyperplanes (FOR PLOTTING); 1 x (N+1)
    %   ncX: the normal to the entry cone hyperplanes; (N+1) x 2
    %   pdock: the apex of the entry cones; 2 x (N+1)
    %   xt: the chaser holding state over the horizon; 6 x (N+1)
    %   cdock: the chaser docking state over the horizon; 6 x (N+1)

    
    % Create the multiplying array in terms of Td
    % Tdarray = 0:Td:N*Td; % 1 x (N+1)
    Tdarray =  linspace(0,N,N+1)*Td;
    
    % Compute the new target orientation over the horizon
    thetat = thetat0 .* ones(1,length(Tdarray)) + thetat_dot .* Tdarray; % (1 x (N+1))
    
    % Calculate the entry cone slopes over the horizon
    m1 = tan(thetat + ecoa + thetah);
    m2 = tan(thetat + ecoa - thetah);
    
    % The normal vector to the entry cone hyperplanes over the horizon
    nc1 = [sin(thetat + ecoa + thetah); -cos(thetat + ecoa + thetah)].';
    nc2 = [sin(thetat + ecoa - thetah); -cos(thetat + ecoa - thetah)].';
    
    % Calculate the state of the docking port on the target (apex), used for
    % the entry cone constraints
    pdockx = ptar(1,:) + pdockr*cos(thetat + pdockoa); 
    pdocky = ptar(2,:) + pdockr*sin(thetat + pdockoa); 
    
    % Calculate the state of the center of the chaser spacecraft during
    % docking; i.e. docking condition
    vtararray = repmat(vtar,1,N+1); % , (2 x (N+1))
    
    cdockx = ptar(1,:) + rdock*cos(thetat + cdockoa);
    cdocky = ptar(2,:) + rdock*sin(thetat + cdockoa);
    cdockxdot = vtararray(1,:) - rdock*thetat_dot*sin(thetat + cdockoa);
    cdockydot = vtararray(2,:) + rdock*thetat_dot*cos(thetat + cdockoa);
    
    % Calculate the desired state for the chaser over the horizon, to be
    % used in the objective function; i.e. holding state
    xtdockoa = atan(rhold/cdockdeltax);
    xtx = ptar(1,:) + rhold*cos(thetat + xtdockoa);
    xty = ptar(2,:) + rhold*sin(thetat + xtdockoa);
    xtxdot = vtararray(1,:) - rhold*thetat_dot*sin(thetat + xtdockoa);
    xtydot = vtararray(2,:) + rhold*thetat_dot*cos(thetat + xtdockoa);
    
    % Combine into arrays for output
    pdock = [pdockx; pdocky]; % 2xN+1 vector
    xt = [xtx; xty; thetat + pi; xtxdot; xtydot; thetat_dot*ones(1,N+1)]; % 6xN+1 vector
    cdock = [cdockx; cdocky; thetat + pi; cdockxdot; cdockydot; thetat_dot*ones(1,N+1)]; % 6xN+1 vector
end



