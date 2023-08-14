function [cl,cu,bineqvec] = theConstraintLimits4(Umat,Fmat,Wmat,con0,con,D1,Dt,...
    ro_1_9,ro_t_9,robj1_9,rtar_9,G1,G2,pdock_9)
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

    % State-Control Limits
    %   cu = [umax umax taumax Inf Inf Inf Inf Inf Inf umax umax ... ]
    Fmat_9 = Umat.' * Fmat; % 9N x 1
    temp = inf*(Fmat_9 == 0);
    temp(isnan(temp)) = 0;
    cu = temp + Fmat_9;
    cl = -cu;
    
    % Defining the upper limits of the constraints, obstacles and entry cone
    bineqvec = zeros(size(Wmat,1), 1);
    % bineqvec(1 : con0(4)) = -ones(con0(4),1) - D1*(ro_1_9 + robj1_9);
    % bineqvec(con0(4)+1 : sum(con0(4:5))) = -ones(con0(5),1) - Dt*(ro_t_9 + rtar_9);
    N = 15;
    bineqvec(1 : N) = -ones(N,1) - D1*(ro_1_9 + robj1_9);
    bineqvec(N+1 : 2*N) = -ones(N,1) - Dt*(ro_t_9 + rtar_9);
    
    if sum(con) > sum(con0) % if entry cone constraints turned on
        bineqvec(sum(con0(4:5))+1 : sum(con(4:6))) = -G1 * pdock_9;
        bineqvec(sum(con(4:6))+1 : sum(con(4:7))) = G2 * pdock_9;
    end
    
end