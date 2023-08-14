% The following script is the initializer for SPOT 4.0; in this script,
% users define all initials parameters and/or constants required for
% simulation and experiment.

clear;
clc;
close all force;

warning('off','all')

%% Start the graphical user interface:

run('GUI_v4_0_Main');

%% Place any custom variables or overwriting variables in this section

% As an example, here are the control parameters for the manipulator.

% Set torque limits on joints

Tz_lim_sharm                   = .1; % Shoulder Joint [Nm]

Tz_lim_elarm                   = .1; % Elbow Joint [Nm]

Tz_lim_wrarm                   = .1; % Wrist Joint [Nm]

% Transpose Jacobian controller gains:

Kp = [0.08 0 0
      0    0.08 0
      0    0    0.002];
Kv = [0.05 0 0
      0    0.05 0
      0    0    0.005];

% Initialize the PID gains for the ARM:

Kp_sharm                       = 1.5;
Kd_sharm                       = 1.0;

Kp_elarm                       = 1.2;
Kd_elarm                       = 0.8;

Kp_wrarm                       = 2;
Kd_wrarm                       = 0.6;


% Define the model properties for the joint friction:
% Based on https://ieeexplore.ieee.org/document/1511048

%Shoulder
Gamma1_sh = 0.005; 
Gamma2_sh = 5;
Gamma3_sh = 40;
Gamma4_sh = 0.015; 
Gamma5_sh = 800; 
Gamma6_sh = 0.005;

%Elbow
Gamma1_el = 0.12; 
Gamma2_el = 5;
Gamma3_el = 10;
Gamma4_el = 0.039; 
Gamma5_el = 800;
Gamma6_el = 0.000001;

%Wrist
Gamma1_wr = 0.025;
Gamma2_wr = 5;
Gamma3_wr = 40;
Gamma4_wr = 0.029;
Gamma5_wr = 800; 
Gamma6_wr = 0.02;

%% Initialize

ALGORITHM = 7;
isCONTROL = 0;

% Converting from degrees to radians and vis versa:

d2r                            = pi/180;
r2d                            = 180/pi;

% Key variables of the optimization problem
Td = 2;  % s
N = 15;  % number of horizons

% Run the file that defines some generic problem variables
theInitializer;

% Position and velocity of obstacles
testcase = 'B'; % 'A' or 'B' or 'C' ...

if testcase == 'A'
    cc = [3.0; 2.0]; % chaser center (initial)
    thetac0 = 0; % rad, initial chaser angle

    ct = [0.5; 0.5]; % target center
    vtar = [0.01;0.01]; % target COM velocity, m/s
    thetat0 = -pi/4; % rad, initial target angle
    thetat_dot = 1.0 * d2r; % rad/s, rotation rate of target spacecraft
elseif testcase == 'B'
    cc = [3.3; 1.8]; % chaser center (initial)
    thetac0 = 0; % rad, initial chaser angle

    ct = [2.0; 2.2]; % target center
    vtar = [-0.01;-0.005]; % target COM velocity, m/s
    thetat0 = -pi/2; % rad, initial target angle
    thetat_dot = -1.25 * d2r; % rad/s, rotation rate of target spacecraft

elseif testcase == 'C'
    cc = [1.75; 0.8]; % chaser center (initial)
    thetac0 = 3.5-2*pi; % rad, initial chaser angle

    ct = [2.9; 1.25]; % target center
    vtar = [-0.01;0.005]; % target COM velocity, m/s
    thetat0 = -1.25; % rad, initial target angle
    thetat_dot = -2.0 * d2r; % rad/s, rotation rate of target spacecraft
    
elseif testcase == 'D'
    cc = [1.5; 2.0]; % chaser center (initial)
    thetac0 = pi; % rad, initial chaser angle

    ct = [2.2; 0.2]; % target center
    vtar = [-0.005;0.005]; % target COM velocity, m/s
    thetat0 = -pi/2; % rad, initial target angle
    thetat_dot = 2.5 * d2r; % rad/s, rotation rate of target spacecraft
end

if ~isCONTROL
    if testcase == 'A'
        pobj1 = [1.5; 1.25];          % position of obstacle 1, Test Case A
        vobj1 = [0.01;0.01];          % velocity of the moving obstacle 1
    elseif testcase == 'B'
        pobj1 = [2.5; 1.5];          % position of obstacle 1, Test Case B
        vobj1 = [-0.01;-0.01];        % velocity of the moving obstacle 1
    elseif testcase == 'C'
        pobj1 = [2.6; 0.6];          % position of obstacle 1, Test Case C
        vobj1 = [-0.006;0.0030];        % velocity of the moving obstacle 1
    elseif testcase == 'D'
        pobj1 = [2.5; 1.2];          % position of obstacle 1, Test Case D
        vobj1 = [-0.01;0.00];        % velocity of the moving obstacle 1
    end
else
    pobj1 = [10.; 10.];          
    vobj1 = [0.00;0.00];
end

% Initial conditions of the chaser
xi = [cc(1),cc(2),thetac0,0,0,0].'; % chaser current state, this will be overwritten in loop

% Current docking condition, the state of the COM of the chaser spacecraft  
cdock = [ct(1)+rdock*cos(thetat0+cdockoa),...
         ct(2)+rdock*sin(thetat0+cdockoa),...
         thetat0+pi,...
         vtar(1)-rdock*thetat_dot*sin(thetat0+cdockoa),...
         vtar(2)+rdock*thetat_dot*cos(thetat0+cdockoa),...
         thetat_dot].';   

init_states_RED           = xi(1:3).'; % [m; m; rad]
init_states_BLACK         = [ct;thetat0].';      % [m; m; rad]
init_states_BLUE          = [pobj1;0].';      % [m; m; rad]

% Setting the upper and lower limits of the constraint equations
% There are:
%     State equation constraints (xdot = Ax+Bu): n*N
%     Thrust constraints (umax): 2*N
%     Thrust constraints (taumax): N
%     Obstacle avoidance w/ 1 moving obstacles: N
%     Obstacle avoidance w/ rotating target: N
[n,~] = size(xi);
con0 = [n*N, 2*N, N, N, N]; % The initial number of constraints
n_constraints0 = sum(con0);

% Initial Starting Guess (Required for Nonlinear Problems)
guess = [-0.14, 0, 0].';
x0 = theInitialGuesser(xi,guess,N,Td,m,Iz);
maxiterations = 10;

% Values needed for the loop. The tolerance sets whether or not the chaser
% has 'docked'. The counting variable is used for 'saving data' purposes.
tol = 0.00015; % m^2