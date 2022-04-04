function [com_totalcontrol,com_totalcost] = theComparator(rc,ux,uy,tau,solJ)
    % Calculates the total control force and total cost for the solution
    % that spans k iterations.
    % Inputs:
    %       rc: width of the chaser spacecraft, m
    %       ux: 1 x (k+1) vector of x control forces, N
    %       uy: 1 x (k+1) vector of y control forces, N
    %       tau: 1 x (k+1) vector of torques, N m
    %       solJ: 1 x (k+1) vector of costs

    % ---------- force summation ------------ %
    com_ux = ux;
    com_uy = uy;
    com_tau = tau;
    com_ut = com_tau / rc;
    
    com_totalcontrol = sum(abs(com_ux)) + sum(abs(com_uy)) + sum(abs(com_ut));
    
    % ---------- total cost ---------- %
    com_totalcost = sum(solJ);
    
end