function dydt = theStateEquations(t,y,prop_u,prop_m,prop_Iz)
    % This function is used in thePropagator and theInitialGuesser. 
    % It simply defines the derivative of the state vector, i.e. 
    %       Y = [x,y,theta,dx,dy,omega]
    % so    dY/dt = [dx,dy,omega,dx/dt,dy/dt,domega/dt]

    dydt = [y(4);y(5);y(6);prop_u(1)/prop_m;prop_u(2)/prop_m;prop_u(3)/prop_Iz];
    
end