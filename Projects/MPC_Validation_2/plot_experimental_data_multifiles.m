% ------------ The Animator (for Experimental Data) ------------ %

% This function takes the data output from the SPOT experiment and produces
% an animation of the trajectory, velocities, and control forces for the
% particular solution. The animation is then saved in .avi format.
%
% This file relies upon the following variables being present in the
% workspace: (note that N is the number of horizons, k is the number of iterations)
%   sol: 1 x (9k+6) vector of state and control variables
%   soliter: (9N+6) x k vector of state and control variable solutions over
%              k iterations
%   solthetat: (k+1) x 1 vector of target spacecraft theta angles
%   solpobj: (k+1) x 4 vector of the x,y positions of the 2 obstacles
%   solptar: (k+1) x 2 vector of the x,y positions of the target
%   solmx: (k+1) x 1 vector of the entry cone hyperplane slopes
%   solJ: 1 x k vector of the cost for each iteration
%   solCPUtime: 1 x k vector of the solver time each iteration
%   KOZ: the KOZ of the obstacles, m
%   rc: the width of the chaser spacecraft, m
%   rt: the width of the target spacecraft, m
%   rcone: the length of the entry cone, m
%   thetat_dot: the angular rotation rate of the target, rad/s
%   wdock: the width of the docking mechanism, m

clear all

DIRECTORY = 'Saved Data/';
EXPDIRECTORY = strcat(DIRECTORY,'ExperimentData_');

% Start a figure
figure

% ================================
% ========== LOAD DATA ===========
% ================================
% filename = '6_16_20_51';
% PRENAME = strcat('ExperimentData_RED_2022_6_20_11_59',filename);
% load(strcat(strcat('Saved Data/',PRENAME),'/dataPacket_SIM.mat'))
%PRENAME = strcat(PRENAME,'_2');

SCENARIO = 'TESTCASEC-QUADPROG-CTRL';

switch SCENARIO
    
    case 'TESTCASEA-QUADPROG-OBST'
        DATA_TIMESTAMP = ["2022_6_21_13_45",...
                          "2022_6_21_13_54",...
                          "2022_6_21_14_7",...
                          "2022_6_21_14_35",...
                          "2022_6_21_14_42"]; 
        DATA_MATNAME = ["MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_2__stitched.mat",...
                        "MPC_Validation_2_2__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = true;
                    
    case 'TESTCASEA-QUADPROG-CTRL'
        DATA_TIMESTAMP = ["2022_6_21_15_6",...
                          "2022_6_21_15_11",...
                          "2022_6_21_15_15",...
                          "2022_6_21_15_30",...
                          "2022_6_21_15_39"]; 
        DATA_MATNAME = ["MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = false;
                    
    case 'TESTCASEB-QUADPROG-OBST'
        
        DATA_TIMESTAMP = ["2022_6_21_16_48",...
                          "2022_6_21_16_55",...
                          "2022_6_21_17_14",...
                          "2022_6_21_17_19",...
                          "2022_6_21_17_24"]; 
        DATA_MATNAME = ["MPC_Validation_2_2__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = true;
        
    case 'TESTCASEB-QUADPROG-CTRL'
        
        DATA_TIMESTAMP = ["2022_6_22_11_9",...
                          "2022_6_22_11_17",...
                          "2022_6_22_11_37",...
                          "2022_6_22_11_41",...
                          "2022_6_22_11_47"]; 
        DATA_MATNAME = ["MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_3__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = true;
        
    case 'TESTCASEC-QUADPROG-OBST'
        
        DATA_TIMESTAMP = ["2022_6_22_12_24",...
                          "2022_6_22_12_33",...
                          "2022_6_22_12_47",...
                          "2022_6_22_12_53",...
                          "2022_6_22_13_5"]; 
        DATA_MATNAME = ["MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = true;
        
    case 'TESTCASEC-QUADPROG-CTRL'
        
        DATA_TIMESTAMP = ["2022_6_22_13_16",...
                          "2022_6_22_13_32",...
                          "2022_6_22_13_44",...
                          "2022_6_22_13_55"]; 
        DATA_MATNAME = ["MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_2__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat",...
                        "MPC_Validation_2_1__stitched.mat"];
        PLOT_OBST = true;

end





startind = 881; % corresponds to 44 s in st_tout, which is the start of Phase3_SubPhase4
step = 20;


for ind = 1:length(DATA_TIMESTAMP)
    
    clearvars -except DIRECTORY EXPDIRECTORY SCENARIO DATA_TIMESTAMP ...
        DATA_MATNAME PLOT_OBST ind startind step
    
    
    TIMESTAMP = DATA_TIMESTAMP(ind);
    MATNAME = DATA_MATNAME(ind);
    DIRECTORY_RED = strcat(strcat(strcat(EXPDIRECTORY,'RED_'),TIMESTAMP),'/');
    FILEPATH = strcat(DIRECTORY_RED,MATNAME);
    load(FILEPATH)
    
    lastind = size(rt_dataPacket,1);
    
    % ================================
    % ========== VARIABLES ===========
    % ================================

    % Define an arbitrary set of angles to be used for plotting the obstacle
    % KOZs
    angles = linspace(0, 2*pi, 500);

    % Assign the state and control vectors, x and u, into separate variables 
    % for each component.
    ani_posx = rt_dataPacket(startind:step:lastind,5);   % 1 x (k+1)
    ani_posy = rt_dataPacket(startind:step:lastind,6);    % 1 x (k+1)
    ani_theta = rt_dataPacket(startind:step:lastind,7);   % 1 x (k+1)
    ani_velx = rt_dataPacket(startind:step:lastind,8);    % 1 x (k+1)
    ani_vely = rt_dataPacket(startind:step:lastind,9);    % 1 x (k+1)
    ani_omega = rt_dataPacket(startind:step:lastind,10)*180/pi; % converted to deg/s
    ani_ux = rt_dataPacket(startind:step:lastind,75);      % 1 x (k+1)
    ani_uy = rt_dataPacket(startind:step:lastind,76);      % 1 x (k+1)
    ani_tau = rt_dataPacket(startind:step:lastind,77);     % 1 x (k+1)

    % Assign the docking mechanism lengths into separate variables
    ani_rc = 0.15;        % m
    rt = 0.15; % m

    % Assign the parameters for the target spacecraft position, attitude, and
    % velocities into separate variables
    ani_ptar(:,1) = rt_dataPacket(startind:step:lastind,25); % (k+1) x 2
    ani_ptar(:,2) = rt_dataPacket(startind:step:lastind,26);
    ani_targett = rt_dataPacket(startind:step:lastind,27);                % (k+1) x 1
    ani_tar_velx = rt_dataPacket(startind:step:lastind,28);   % 1 x (k+1); m/s
    ani_tar_vely = rt_dataPacket(startind:step:lastind,29);   % 1 x (k+1); m/s
    ani_tar_omega = rt_dataPacket(startind:step:lastind,30)*180/pi; % 1 x (k+1); deg/s

    ani_pobj(:,1) = rt_dataPacket(startind:step:lastind,45); % (k+1) x 2
    ani_pobj(:,2) = rt_dataPacket(startind:step:lastind,46);
    ani_thetaobj(:,2) = rt_dataPacket(startind:step:lastind,47);

    % ================================
    % ========== PLOTTING ============
    % ================================



    
    %plot(ani_targetx(:,i),ani_targety(:,i),'k--', 'LineWidth', 0.5) % Plot predicted solution
    %hold on
    %plot(KOZ*cos(angles)+ani_pobj1(i,1), KOZ*sin(angles)+ani_pobj1(i,2), 'b-'); % Obstacle 1
    %hold on
    %plot(KOZ*cos(angles)+ani_pobj2(i,1), KOZ*sin(angles)+ani_pobj2(i,2), 'b-'); % Obstacle 2
    %hold on
    
    % Calculate endpoints of each of the entry cone hyerplanes
    %entryconex1 = ani_pdockx(i) + rcone*sin(ani_targett(i)+ecoa+thetah)/ani_m1(i); % point for entry cone 1
    %entryconey1 = ani_pdocky(i) + rcone*cos(ani_targett(i)+ecoa+thetah)*ani_m1(i);
    %entryconex2 = ani_pdockx(i) + rcone*sin(ani_targett(i)+ecoa-thetah)/ani_m2(i); % point for entry cone 2
    %entryconey2 = ani_pdocky(i) + rcone*cos(ani_targett(i)+ecoa-thetah)*ani_m2(i);   
    
    %plot([ani_pdockx(i),entryconex1],[ani_pdocky(i),entryconey1],'Color', '#999999') % Entry cone 1
    %hold on
    %plot([ani_pdockx(i),entryconex2],[ani_pdocky(i),entryconey2],'Color', '#999999') % Entry cone 2
    %hold on
    
    chaser = theChaserSpacecraft(ani_posx(end), ani_posy(end), ani_theta(end), ani_rc);
    patch(chaser(:,1), chaser(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'r') % Plot the chaser spacecraft
    hold on
    plot(ani_posx, ani_posy, 'r.-'); % Plot the solution up until i
    hold on
    target = theTargetSpacecraft(ani_ptar(end,1), ani_ptar(end,2), ani_targett(end), rt);
    patch(target(:,1), target(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'k') % Plot the target spacecraft
    hold on    
    plot(ani_ptar(:,1), ani_ptar(:,2), 'k.-'); % Plot the solution up until i
    hold on
    if PLOT_OBST
        obstac = theObstacleSpacecraft(ani_pobj(end,1), ani_pobj(end,2), ani_thetaobj(end), rt);
        patch(obstac(:,1), obstac(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'b') % Plot the target spacecraft
        plot(ani_pobj(:,1), ani_pobj(:,2), 'b.-'); % Plot the solution up until i
        hold on
    end
    
    grid on
    ax = gca;
    ax.FontSize = 6;
    %ylim([min([ani_posy.',ani_ptar(:,2).'])-ani_rc*2 max([ani_posy.',ani_ptar(:,2).'])+ani_rc*2])
    %xlim([min([ani_posx.',ani_ptar(:,1).'])-ani_rc*2 max([ani_posx.',ani_ptar(:,1).'])+ani_rc*2])
    xlim([0,3.5])
    ylim([0,2.4])
    
    %daspect([1 1 1])
    %axis square;
    

    
end

xlabel('x-position [m]')
ylabel('y-position [m]')

saveas(gcf,strcat(strcat(strcat(DIRECTORY,'FIG-'),SCENARIO),'.png'))


  