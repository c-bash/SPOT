% ------------ The Initializer ------------ %

% This file will create the generic variables used for the optimized
% control problem. It is called by main.m.

% The testbed constants are defined by the SRCL GitHub wiki page, 
% https://github.com/Carleton-SRCL/SPOT/wiki/Overview-of-the-Spacecraft-
% Proximity-Operations-Testbed

% ================================
% ====== TESTBED VARIABLES =======
% ================================

% Conversion from degrees to radians
d2r = pi/180;

% Testbed constants
m = 16.95; % kg; mass of testbed chaser (RED)
Iz = 0.271; % kg m^2; moment of inertia of chaser (RED)
rt = 0.15; % m, target (BLACK) radius to docking port
rc = 0.15; %m, chaser (RED) radius to docking port
wdock = 0.05; %m, the length of the docking port
rdock = rt+rc+wdock; %m, the distance of the docking condition from target COM
umax = 0.15; % N, or perhaps closer to 0.25 N from the SPOT wiki page
taumax = umax*rc; % Nm

% ================================
% ======== MODEL DYNAMICS ========
% ================================

% Continuous state space model (double-integrator model)
Ac = [0,0,0,1,0,0; 0,0,0,0,1,0; 0,0,0,0,0,1; 0,0,0,0,0,0; 0,0,0,0,0,0; 0,0,0,0,0,0];
Bc = [0,0,0; 0,0,0; 0,0,0; 1/m,0,0; 0,1/m,0; 0,0,1/Iz];
Cc = [1,1,1,1,1,1];
Dc = [1,1,1];
css = ss(Ac,Bc,Cc,Dc);

% Discrete state space model, we want the Ad and Bd matrices
dss = c2d(css,Td);
Ad = dss.A;
Bd = dss.B;

% Create the large C matrix of the form
% C = [-B  I  0  0  0 ...  0  0  0
%       0 -A -B  I  0 ...  0  0  0
%       0  0  0 -A -B ...  0  0  0
%       0  0  0  0  0 ...  I  0  0
%       0  0  0  0  0 ... -A -B  I]
row1 = zeros(6,9*N);
row1(:,1:3) = -Bd;
row1(:,4:9) = eye(6);

rowi = zeros(6,9*N);
rowi(:,4:9) = -Ad;
rowi(:,10:12) = -Bd;
rowi(:,13:18) = eye(6);

Cmat = zeros(6*N,9*N);
Cmat(1:6,:) = row1;
for i=2:N
    Cmat((i-1)*6+1:i*6,:) = rowi;
    rowi(:,9*(N-1)+1:9*N) = [];
    rowi = [zeros(6,9),rowi];
end

clear row1 rowi

% ================================
% ===== COST FUNCTION WEIGHTS ====
% ================================

% Define the P,Q,R matrices (references to Park 2016,2017)
Q = [10,0,0,0,0,0; 0,10,0,0,0,0; 0,0,10,0,0,0; 0,0,0,300,0,0; 0,0,0,0,300,0; 0,0,0,0,0,500];
R = [10,0,0; 0,10,0; 0,0,10];
[P,~,~] = idare(Ad,Bd,Q,R,[],[]);

% Create the large H matrix of the form 
% H = diag[R,Q,R,Q,R,Q,...,R,Q,R,P]
RQ = blkdiag(R,Q);  % 9x9
RP = blkdiag(R,P);  % 9x9
pattern = repmat({RQ},N,1);
pattern{N} = RP;
H = blkdiag(pattern{:});

clear RQ RP pattern

% ================================
% ======= ROTATING TARGET ========
% ================================

% Rotating target parameters
rhold = 1.15; % m, initial dynamic holding radius
gamma = 0.95; % holding radius scaling factor

% ================================
% ========= ENTRY CONE ===========
% ================================

% Entry cone parameters
rcone = 0.75; % m, entry cone corridor length
thetah = 20.0 * d2r; % rad, entry cone half-angle

% ================================
% ========== OBSTACLES ===========
% ================================

% Obstacle parameters
KOZ = 0.3; % m, the keep out zone of the obstacle
KOZ_actual = KOZ + rc*sqrt(2); % taking into account the width of the chaser
Sj = [KOZ_actual^(-2),0; 0,KOZ_actual^(-2)]; % obstacle shape matrix; circle

% ================================
% ====== THRUST CONSTRAINTS ======
% ================================
Umat = zeros(3*N,9*N);
for row = 1:N
    Umat(row,1+9*(row-1)) = 1;
    Umat(row+N,2+9*(row-1)) = 1;
    Umat(row+2*N,3+9*(row-1)) = 1;
end

Fmat = reshape(repmat([umax,umax,taumax],N,1),3*N,1); % 3N x 1
