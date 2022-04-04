function J = theObjectiveFunction0(x,N,xt,P,Q,R)
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
  % xt is the target position for the chaser at rhold from the target
  % spacecraft. It is an array of size 6xN+1, where the first column is
  % target for the current state and can therefore be ignored in the
  % objective function.
  
  [n,~] = size(xt); % the length of the initial conditions
  Nu = 3; % the number of control forces
  Nxu = n + Nu; % the combined size of the x and u parameters for a single horizon step
  
  %J = (x((Nu+1)+Nxu*(N-2):Nxu*(N-1)) - xt).' * P * (x((Nu+1)+Nxu*(N-2):Nxu*(N-1)) - xt) + (xi-xt).' * Q * (xi-xt) + x(1:Nu).' * R * x(1:Nu);
  %for i = 0:N-2 % for the rest of the horizon
  %    iN = Nxu*i; % this is the iteration number
  %    J = J + (x((Nu+1)+iN:Nxu+iN) - xt).' * Q * (x((Nu+1)+iN:Nxu+iN) - xt) + x((Nxu+1)+iN:(Nxu+Nu)+iN).' * R * x((Nxu+1)+iN:(Nxu+Nu)+iN);
  %end
  
  J = (x(4+9*(N-1):9*N) - xt(:,N+1)).' * P * (x(4+9*(N-1):9*N) - xt(:,N+1)) + x(1:3).' * R * x(1:3);
  for i = 0:N-2 % for the rest of the horizon
      %iN = Nxu*i; % this is the iteration number
      J = J + (x(i*9+4:(i+1)*9) - xt(:,i+2)).' * Q * (x(i*9+4:(i+1)*9) - xt(:,i+2)) + x(i*9+10:i*9+12).' * R * x(i*9+10:i*9+12);
  end

end