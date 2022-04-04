function [cl,cu] = theConstraintLimits3(con0,con,bvec,Fmat,D1,D2,Dt,...
    ro_1_9,ro_2_9,ro_t_9,robj1_9,robj2_9,rtar_9,G1,G2,pdock_9)
    % Returns the limits of the constraint equations.
    %
    % Inputs:
    %   con0: an array of the number of constraints in each of the
    %         following categories (in order): state equations, control 
    %         forces, control torques, obstacle constraints
    %   con: an array of the number of constraint. Is the same as con0 when
    %        there are no entry cone constraints. Otherwise, will be 2N
    %        longer.
    %   Fmat: 3Nx1 array of umax (first 2N) and taumax (last N)
    %   DX: Nx9N matrix of linearized obstacle constraints
    %   ro_X_9: 9Nx1 array of expansion states on the KOZ over the horizon
    %   robjX_9: 9Nx1 array of the obstacle states over the horizon
    %   GX: Nx9N entry cone constraint matrix
    %   pdock9: 9Nx1 array of the docking point states over the horizon
    % Outputs:
    %   cl,cu: arrays of the lower and upper limits, respectively, of the
    %          constraints

    % Defining the lower limits of the constraints
    cl = zeros(sum(con), 1); 
    cl(1:con0(1)) = bvec;
    cl(con0(1)+1 : sum(con0(1:3))) = -Fmat;
    cl(sum(con0(1:3))+1 : end) = -inf;

    % Defining the upper limits of the constraints
    cu = zeros(sum(con), 1);
    cu(1:con0(1)) = bvec;
    cu(con0(1)+1 : sum(con0(1:3))) = Fmat;
    cu(sum(con0(1:3))+1 : sum(con0(1:4))) = -ones(con0(4),1) - D1*(ro_1_9 + robj1_9);
    cu(sum(con0(1:4))+1 : sum(con0(1:5))) = -ones(con0(5),1) - D2*(ro_2_9 + robj2_9);
    cu(sum(con0(1:5))+1 : sum(con0(1:6))) = -ones(con0(6),1) - Dt*(ro_t_9 + rtar_9);
    
    if sum(con) > sum(con0) % if entry cone constraints turned on
        cu(sum(con0(1:6))+1 : sum(con(1:7))) = -G1 * pdock_9;
        cu(sum(con(1:7))+1 : sum(con(1:8))) = G2 * pdock_9;
    end
    
end