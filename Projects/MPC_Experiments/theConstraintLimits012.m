function [cl,cu] = theConstraintLimits012(con0,con,umax,taumax)
    % Returns the limits of the constraint equations.
    %
    % Inputs:
    %   con0: an array of the number of constraints in each of the
    %         following categories (in order): state equations, control 
    %         forces, control torques, obstacle constraints
    %   umax: scalar of the maximum control force
    %   taumax: scalar of the maximum control torque
    % Outputs:
    %   cl,cu: arrays of the lower and upper limits, respectively, of the
    %          constraints

    % Defining the lower limits of the constraints
    cl = zeros(sum(con), 1); 
    cl(1:con0(1)) = 0;
    cl(con0(1)+1 : sum(con0(1:2))) = -umax;
    cl(sum(con0(1:2))+1 : sum(con0(1:3))) = -taumax;
    cl(sum(con0(1:3))+1 : end) = -inf;

    % Defining the upper limits of the constraints
    cu = zeros(sum(con), 1);
    cu(1:con0(1)) = 0;
    cu(con0(1)+1 : sum(con0(1:2))) = umax;
    cu(sum(con0(1:2))+1 : sum(con0(1:3))) = taumax;
    cu(sum(con0(1:3))+1 : end) = 0;
end