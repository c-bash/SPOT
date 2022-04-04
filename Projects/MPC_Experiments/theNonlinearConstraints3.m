function nl = theNonlinearConstraints3(x,Wmat)
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
    % Inputs:
    %   x: the optimization state and control vector; (9N x 1)
    %   Wmat: the large constraint matrix; (w/o EC: 12N x 9N; w/ EC: 14N x 9N)
    % Outputs:
    %   nl: the constraint equations
        
    % Constraints for the dynamics of the problem, dotx-Ax-Bu  
    nl = Wmat * x;

end