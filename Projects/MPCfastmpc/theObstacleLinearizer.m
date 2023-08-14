function ro = theObstacleLinearizer(px,pobj,KOZ)
    % Returns the point on the obstacle KOZ, ro, that is between the
    % obstacle position and the chaser trajectory over the course of the
    % horizon at the current iteration step. 
    %
    % Inputs:
    %   px: the chaser reference trajectpry; 2 x N
    %   pobj: the position of the object over the trajectory; 2 x N
    %   KOZ: the KOZ radius of the obstacle, m
    % Outputs:
    %   ro: the position on the KOZ of the obstacle between the obstacle
    %       center and the reference trajctory of the chaser; 2 x N
    
    
    % Calculate the vector between the chaser and obstacle positions over
    % the horizon
    dr = px - pobj;   % 2xN matrix
    deltax = dr(1,:); % 1xN matrix
    deltay = dr(2,:); % 1xN matrix
    
    % Calculate the slope of the line between the chaser and obstacle
    % positions over the horizon
    m =  deltay ./ deltax ;
    
    % Calculate the angle of the vector wrt the obstacle center. Use the
    % fact that arctan only ranges from (-pi) to pi, so any angles in
    % quadrants 2 and 3 must have a factor of pi added.
    inQ2Q3 = deltax < 0; % logical array, whether the angle is in q2/q3 or not
    theta = pi*inQ2Q3 + atan(m);
    
    % Calculate the point ro, which is the point on the KOZ of the obstacle
    % inbetween the chaser and obstacle positions
    ro = pobj + KOZ*[cos(theta);sin(theta)];

end
