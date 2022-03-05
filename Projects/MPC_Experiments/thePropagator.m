% ------------ The Propagator ------------ %

% This function takes the control forces of the output solution and
% propagates the dynamics of the soacecraft, xdot = Ax+Bu using ODE45, 
% the Runge-Kutta 4(5) algorithm. The propagated solution is then checked
% with the output trajectory of the solver to see if the dynamics
% constraints were preserved.
%
% Note: Run theAnimator before running this script.
%
% This file relies upon the following variables being present in the
% workspace: (note that N is the number of horizons, k is the number of iterations)
%   sol: 1 x (9k+6) vector of state and control variables
%   Td: the discretized time, s
%   m: the mass of the chaser spacecraft, kg
%   Iz: the moment of inertia of the chaser spacecraft, kg m^2


% ================================
% ========== VARIABLES ===========
% ================================

% Define a stepsize for the RK4(5) algorithm (less than Td). Use the
% stepsize to get an array of time values up until Td.
prop_stepsize = 0.1;    % s
prop_tspan = 0:prop_stepsize:Td;

% Define the chaser mass/intertia properties
prop_m = m;
prop_Iz = Iz;

% Assign the initial starting position of the chaser to a variable
prop_y0 = sol(1:6).';
prop_y00 = prop_y0.'; % 1 x 6; to be prepended to the final solution

% Assign the control forces to a set of variables
prop_ux = sol(7:9:end-3);  % 1 x k
prop_uy = sol(8:9:end-3);  % 1 x k
prop_tau = sol(9:9:end-3); % 1 x k

% Make a vector of the control forces
prop_controls = [prop_ux;prop_uy;prop_tau]; % 3 x k

% The propagated solution will be the size of the steps within Td, 
% (prop_sizetspan), multiplied by k
prop_sizetspan = length(prop_tspan)-1;
prop_sizeux = length(prop_ux);

% ================================
% ========= PROPAGATION ==========
% ================================

% Initialize array for solution
prop_sol = zeros(prop_sizetspan*prop_sizeux,6); % (k*Td/deltat) x 6

% Loop through the k iterations
for i = 1:prop_sizeux
    % Set the control forces at the ith iteration
    prop_u = prop_controls(:,i); % 3 x 1
    
    % Propagate the control forces for the Td/deltat timesteps within the
    % iterval
    [prop_t,prop_y] = ode45(@(t,y) theStateEquations(t,y,prop_u,prop_m,prop_Iz),prop_tspan,prop_y0);

    % Add the propagated solution (from index 2 and onwards) into the
    % solution array. The reason why we don't include the first index is
    % that we don't want overlap between the final point of the ith
    % iteration and the first point of the (i+1)th iteration.
    prop_sol(1+(i-1)*prop_sizetspan:i*prop_sizetspan,:) = prop_y(2:end,:);
    
    % Update the current state to be the final state of this solution
    prop_y0 = prop_y(end,:);
end

% Prepend the initial starting state
prop_sol = [prop_y00;prop_sol]; % (k*Td/deltat + 1) x 6


% ================================
% =========== PLOTTING ===========
% ================================

% Assign each of the propagated solutions parameters a variable
prop_posx = prop_sol(:,1);          % (k*Td/deltat + 1) x 1
prop_posy = prop_sol(:,2);          % (k*Td/deltat + 1) x 1
prop_theta = prop_sol(:,3);         % (k*Td/deltat + 1) x 1
prop_velx = prop_sol(:,4);          % (k*Td/deltat + 1) x 1
prop_vely = prop_sol(:,5);      	% (k*Td/deltat + 1) x 1
prop_omega = prop_sol(:,6)*180/pi;  % (k*Td/deltat + 1) x 1; deg/s

% To compare the propagated solution and the solver solution, we need to
% compare them over the same x-scale, so we write the equivalent elapsed
% time for the solution
ani_t = 0:Td:Td*length(ani_posx)-1; % 0:Td:Td*k
prop_t = 0:prop_stepsize:(length(prop_posx)-1)*prop_stepsize; % 0:deltat:Td*k

% Start a figure
figure('Position',[300,150,700,400])

% -------- X-position ---------- %
subplot(3,2,1)
plot(ani_t,ani_posx,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_posx,'r--','LineWidth',1)
grid on
set(gca, 'XTickLabel', [])
ylabel('X-position [m]')
legend('MPC','ODE45')

% -------- Y-position ---------- %
subplot(3,2,3)
plot(ani_t,ani_posy,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_posy,'r--','LineWidth',1)
grid on
set(gca, 'XTickLabel', [])
ylabel('Y-position [m]')

% -------- THETA-angle ---------- %
subplot(3,2,5)
plot(ani_t,ani_theta,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_theta,'r--','LineWidth',1)
grid on
ylabel('Attitude [rads]')
xlabel('Time [s]')

% -------- X-velocity ---------- %
subplot(3,2,2)
plot(ani_t,ani_velx,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_velx,'r--','LineWidth',1)
grid on
set(gca, 'XTickLabel', [])
ylabel('X-velocity [m/s]')

% -------- Y-velocity ---------- %
subplot(3,2,4)
plot(ani_t,ani_vely,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_vely,'r--','LineWidth',1)
grid on
set(gca, 'XTickLabel', [])
ylabel('Y-velocity [m/s]')

% -------- THETA-velocity ---------- %
subplot(3,2,6)
plot(ani_t,ani_omega,'k','LineWidth',0.5)
hold on
plot(prop_t,prop_omega,'r--','LineWidth',1)
grid on
ylabel('Angular velocity [deg/s]')
xlabel('Time [s]')

savefig(strcat(PRENAME, 'propagator.fig'))


    