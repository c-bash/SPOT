function J = theObjectiveFunction123(x,xt,H)
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
  % This function defines the cost function for the optimization problem.
  %
  % Inputs:
  %   x: the optimization state and control vector; (9N x 1)
  %   xt: the chaser holding state over the horizon (ignoring the current 
  %       state at index 1) reshaped from 6 x (N+1) to be in the format
  %       [03 xt(k+1) 03 xt(k+2) 03 ... 03 xt(k+N)]; (9N x 1)
  %   H: the block diagonal matrix of [R,Q,R,Q,...,R,Q,R,P]; (9N x 9N)
  % Outputs:
  %   J: the cost at the current iteration (scalar)
  
  ztilde = x - xt;
  
  J = ztilde.' * H * ztilde;

end