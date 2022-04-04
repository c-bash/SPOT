% The following script is the initializer for SPOT 3.0; in this script,
% users define all initials parameters and/or constants required for
% simulation and experiment.

% Version: 3.07 (Beta Release)

% Authors: Alexander Crain
% Legacy: David Rogers & Kirk Hovell

clear;
clc;
close all force;
addpath(genpath('../../Custom_Library'))

warning('off','all')

fprintf('|----------------------------------------------------------------|\n')
fprintf('|----------------------------------------------------------------|\n')
fprintf('|                       -------------------                      |\n')
fprintf('|                     | Welcome to SPOT 3.0 |                    |\n')
fprintf('|                       -------------------                      |\n')
fprintf('|                                                                |\n')
fprintf('|Authors (v3.0): Alex Crain                                      |\n')
fprintf('|Authors (v2.0): Alex Crain and Kirk Hovell                      |\n')
fprintf('|Authors (Legacy): Dave Rogers and Kirk Hovell                   |\n')
fprintf('|                                                                |\n')
fprintf('|Current Version: 3.07 (Beta Release)                            |\n')
fprintf('|                                                                |\n')
fprintf('|Last Edit: 2021-03-07                                           |\n')
fprintf('|                                                                |\n')
fprintf('|----------------------------------------------------------------|\n')
fprintf('|----------------------------------------------------------------|\n')

%% User-defined constants:

% Converting from degrees to radians and vis versa:

d2r                            = pi/180;
r2d                            = 180/pi;

% Initialize the table size for use in the GUI (don't delete):

xLength                        = 3.51155;   % [m]
yLength                        = 2.41935;   % [m]

% Initialize the PID gains for the RED platform:

Kp_xr                          = 2;
Kd_xr                          = 5;

Kp_yr                          = 2;
Kd_yr                          = 5;

Kp_tr                          = 0.1;
Kd_tr                          = 0.4;

% Initialize the PID gains for the BLACK platform:

Kp_xb                          = 2;
Kd_xb                          = 5;

Kp_yb                          = 2;
Kd_yb                          = 5;

Kp_tb                          = 0.1;
Kd_tb                          = 0.4;

% Initialize the PID gains for the BLUE platform:

Kp_xblue                       = 2;
Kd_xblue                       = 5;

Kp_yblue                       = 2;
Kd_yblue                       = 5;

Kp_tblue                       = 0.1;
Kd_tblue                       = 0.4;

% Set the noise variance level for the RED and BLACK platforms:

noise_variance_RED             = 0;
noise_variance_BLACK           = 0;
noise_variance_BLUE            = 0;

%% Set the base sampling rate: 

% This variable will change the frequency at which the template runs. If
% the frequency of the template changes, the frequency of the server must
% also be changed, i.e. open the StreamData.sln under the PhaseSpace Server
% folder, and change line 204 from owl.frequency(10) to 
% owl.frequency(serverRate):

baseRate                       = 0.05;      % 20 Hz

%% Set the frequency that the data is being sent up from the PhaseSpace:

% This variable must be less then the baseRate; in simulation, setting this
% equal to the baseRate causes the simulation fail, while in experiment
% setting this equal to or higher then the baseRate causes the data to
% buffer in the UDP send.

serverRate                     = 0.1;       % 10 Hz

%% Set the duration of each major phase in the experiment, in seconds:

Phase0_Duration                = 5;        % [s]
Phase1_Duration                = 5;         % [s]
Phase2_Duration                = 26;        % [s]
Phase3_Duration                = 130;        % [s]
Phase4_Duration                = 20;        % [s]
Phase5_Duration                = 5;         % [s]

% Set the duration of the sub-phases. Sub-phases occur during the
% experiment phase (Phase3_Duration) and must be manually inserted into the
% diagram. The total duration of the sub-phases must equal the length of
% the Phase3_Duration.

Phase3_SubPhase1_Duration      = 4;        % [s]
Phase3_SubPhase2_Duration      = 2;        % [s]
Phase3_SubPhase3_Duration      = 2;        % [s]
Phase3_SubPhase4_Duration      = 122;        % [s]

% Determine the total experiment time from the durations:

tsim                           = Phase0_Duration + Phase1_Duration + ...
                                 Phase2_Duration + Phase3_Duration + ...
                                 Phase4_Duration + Phase5_Duration;        

% Determine the start time of each phase based on the duration:

Phase0_End                     = Phase0_Duration;
Phase1_End                     = Phase0_Duration + Phase1_Duration;           
Phase2_End                     = Phase0_Duration + Phase1_Duration + ...
                                 Phase2_Duration;         
Phase3_End                     = Phase0_Duration + Phase1_Duration + ...
                                 Phase2_Duration + Phase3_Duration;      
Phase4_End                     = Phase0_Duration + Phase1_Duration + ...
                                 Phase2_Duration + Phase3_Duration + ...
                                 Phase4_Duration; 
Phase5_End                     = Phase0_Duration + Phase1_Duration + ...
                                 Phase2_Duration + Phase3_Duration + ...
                                 Phase4_Duration + Phase5_Duration;                              
                             
% Determine the start time of each sub-phase based on the duration:  

Phase3_SubPhase1_End           = Phase2_End + Phase3_SubPhase1_Duration;
Phase3_SubPhase2_End           = Phase2_End + Phase3_SubPhase1_Duration + ...
                                 Phase3_SubPhase2_Duration;
Phase3_SubPhase3_End           = Phase2_End + Phase3_SubPhase1_Duration + ...
                                 Phase3_SubPhase2_Duration +...
                                 Phase3_SubPhase3_Duration;
Phase3_SubPhase4_End           = Phase2_End + Phase3_SubPhase1_Duration + ...
                                 Phase3_SubPhase2_Duration +...
                                 Phase3_SubPhase3_Duration +...
                                 Phase3_SubPhase4_Duration;                             
                          
%% Load in any required data:

% Define the mass properties for the RED, BLACK, and BLUE platforms:

model_param(1)                 = 16.9478; % RED Mass
model_param(2)                 = 0.2709;  % RED Inertia;
model_param(3)                 = 12.3341; % BLACK Mass
model_param(4)                 = 0.1880;  % BLACK Inertia
model_param(5)                 = 12.7621; % BLUE Mass
model_param(6)                 = 0.1930;  % BLUE Inertia

% Initialize the thruster positions for the RED, BLACK, and BLUE platforms,
% as well as the expected maximum forces. The expected forces will only 
% affect the simulations.

F_thrusters_RED               = 0.25.*ones(8,1);
F_thrusters_BLACK             = 0.25.*ones(8,1);
F_thrusters_BLUE              = 0.25.*ones(8,1);
thruster_dist2CG_RED          = [49.92;-78.08;70.46;-63.54;81.08;-50.42;57.44;-75.96];
thruster_dist2CG_BLACK        = [83.42;-52.58;55.94;-60.05;54.08;-53.92;77.06;-55.94];
thruster_dist2CG_BLUE         = [83.42;-52.58;55.94;-60.05;54.08;-53.92;77.06;-55.94];
                                              
%% FROM MAIN

% ================================
% ========== ALGORITHM ===========
% ================================

% Choose which algorithm configuration to use to solve the problem
%    NMPC (IPOPT): 0
%    Compact NMPC (IPOPT): 1
%    Linearized MPC (IPOPT): 2
%    Linearized Compact MPC (IPOPT): 3
%    Linearized Compact MPC (Fmincon): 4
%    Linearized Compact MPC (Quadprog): 5
ALGORITHM = 5;

if ALGORITHM == 0
    % NMPC
    % IPOPT
    % Objective function = for loop 
    % Nonlinear constraints = for loop
    algorithmname = 'NMPC';
elseif ALGORITHM == 1
    % NMPC
    % IPOPT
    % Objective function = BIG matrix multiply 
    % Nonlinear constraints = small matrix multiply and assignment
    algorithmname = 'CompactNMPC';
elseif ALGORITHM == 2
    % MPC
    % IPOPT
    % Objective function = BIG matrix multiply 
    % Linear constraints = small matrix multiply and assignment
    algorithmname = 'LinearizedMPC';
elseif ALGORITHM == 3
    % MPC
    % IPOPT
    % Objective function = BIG matrix multiply 
    % Linear constraints = BIG matrix multiply
    algorithmname = 'LinearizedCompactMPC';
elseif ALGORITHM == 4
    % MPC
    % Fmincon
    % Objective function = BIG matrix multiply 
    % Linear constraints = big matrix multiply, split into = and <= constraints
    algorithmname = 'Fmincon';
elseif ALGORITHM == 5
    % MPC
    % Quadprog
    % Objective function = Split into sections, 2(0.5*x.'*H*x - xt.'*H*x) + xt.'*H*xt
    % Linear constraints = big matrix multiply, split into = and <= constraints
    algorithmname = 'Quadprog';
end

% ================================
% ========== VARIABLES ===========
% ================================

% Key variables of the optimization problem
Td = 2;  % s
N = 15;  % number of horizons

% Run the file that defines some generic problem variables
theInitializer;

% Initial conditions of the chaser
cc = [3; 2]; % chaser center (initial)
thetac0 = pi; % rad, initial chaser angle
xi = [cc(1),cc(2),thetac0,0,0,0].'; % chaser current state, this will be overwritten in loop

% Initial conditions of the target
ct = [1.; 1.]; % target center
vtar = [0.01;0.00]; % target COM velocity, m/s
thetat0 = pi/4; % rad, initial target angle
thetat_dot = 2.0 * d2r; % rad/s, rotation rate of target spacecraft

% Current docking condition, the state of the COM of the chaser spacecraft  
cdock = [ct(1)+rdock*cos(thetat0),...
         ct(2)+rdock*sin(thetat0),...
         thetat0+pi,...
         -rdock*thetat_dot*sin(thetat0),...
         rdock*thetat_dot*cos(thetat0),...
         thetat_dot].';     

% Position and velocity of obstacles
testcase = 'A'; % 'A' or 'B' or 'C' ...

vobj1 = [-0.01;0.01];             % velocity of the moving obstacle 1 (0.01,0.01)
vobj2 = [0.0;0.0];            % velocity of the moving obstacle 2 (-0.02,0.02)
if testcase == 'A'
    pobj1 = [2.25; 1.25];          % position of obstacle 1, Test Case A
    pobj2 = [4.; 4.];          % position of obstacle 2, Test Case A
elseif testcase == 'B'
    pobj1 = [2.0; 2.5];          % position of obstacle 1, Test Case B
    pobj2 = [3.0; 2.5];          % position of obstacle 2, Test Case B
elseif testcase == 'C'
    pobj1 = [2.5; 2.5];          % position of obstacle 1, Test Case C
    pobj2 = [1.5; 2.0];          % position of obstacle 2, Test Case C
elseif testcase == 'D'
    pobj1 = [1.0; 2.5];          % position of obstacle 1, Test Case D
    pobj2 = [3.0; 2.0];          % position of obstacle 2, Test Case D
elseif testcase == 'E'
    pobj1 = [0.0; 0.0];          % position of obstacle 1, Test Case E
    pobj2 = [0.0; 0.0];          % position of obstacle 2, Test Case E
elseif testcase == 'F'
    pobj1 = [4.0; 4.0];          % position of obstacle 1, Test Case F
    pobj2 = [4.0; 4.0];          % position of obstacle 2, Test Case F
elseif testcase == 'G'
    pobj1 = [3.0; 1.0];          % position of obstacle 1, Test Case G
    pobj2 = [2.0; 0.5];          % position of obstacle 2, Test Case G
end

% Filename for output files, animations, figures, etc.
PRENAME = strcat('1-TestCase',testcase,'-',algorithmname,'-');

% Setting the upper and lower limits of the constraint equations
% There are:
%     State equation constraints (xdot = Ax+Bu): n*N
%     Thrust constraints (umax): 2*N
%     Thrust constraints (taumax): N
%     Obstacle avoidance w/ 2 moving obstacles: 2*N
%     Obstacle avoidance w/ rotating target: N
[n,~] = size(xi);
con0 = [n*N, 2*N, N, N, N, N]; % The initial number of constraints
n_constraints0 = sum(con0);

% Initial Starting Guess (Required for Nonlinear Problems)
guess = [-0.15, -0.15, 0].';
x0 = theInitialGuesser(xi,guess,N,Td,m,Iz);

% ================================
% ===== SOLVING THE PROBLEM ======
% ================================

% MOVED OPTIONS ...

% Values needed for the loop. The tolerance sets whether or not the chaser
% has 'docked'. The counting variable is used for 'saving data' purposes.
tol = 0.0001; % m^2
count = 1;

% This loop runs until the current position of the chaser matches the desired
% position of the COM as set by cdock (to within the set tolerance)
% COMMENT: should this loop also include the angular difference between
% chaser and target?
isdocked = 0;

%%  Set the drop, initial, and home positions for each platform:

init_states_RED           = xi(1:3); % [m; m; rad]
init_states_BLACK         = [ct;thetat0];      % [m; m; rad]
init_states_BLUE          = [pobj1;0];      % [m; m; rad]

drop_states_RED           = rand([3,1])/10+init_states_RED; % [m; m; rad]
drop_states_BLACK         = rand([3,1])/10+init_states_BLACK;  % [m; m; rad]
drop_states_BLUE          = rand([3,1])/10+init_states_BLUE;         % [m; m; rad]

home_states_RED           = [ xLength/2+0.7; yLength/2; pi]; % [m; m; rad]
home_states_BLACK         = [ xLength/2; yLength/2; 0];  % [m; m; rad]
home_states_BLUE          = [ xLength/2-0.9; yLength/2+0.5; 0];  % [m; m; rad]

%% Load in BLACK control data from previous experiment 
% Go to SimulationData_2022_4_2_15_5
load('Saved Data/SimulationData_2022_4_2_15_5/dataPacket_SIM.mat')
controlsblack = dataPacket(:,[1,65,66,67]);

%% Start the graphical user interface:

run('GUI_v3_07')


