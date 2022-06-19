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

% ================================
% ========== LOAD DATA ===========
% ================================
filename = '6_16_20_51';
PRENAME = strcat('SimulationData_2022_',filename);
load(strcat(strcat('Saved Data/',PRENAME),'/dataPacket_SIM.mat'))
%PRENAME = strcat(PRENAME,'_2');

startind = 881; % corresponds to 44 s in st_tout, which is the start of Phase3_SubPhase4
step = 20;
lastind = 3323;

% ================================
% ========== VARIABLES ===========
% ================================

% Define an arbitrary set of angles to be used for plotting the obstacle
% KOZs
angles = linspace(0, 2*pi, 500);

% Ensure that the sol variable is 1 x (9N), by adding 3 NaN values in the
% place of the last 3 control variables (i.e., when the chaser has docked)
%if rem(length(sol),9) > 0
%    sol(end+1:end+3) = [NaN,NaN,NaN]; % now 1 x (9(k+1))
%end

% Assign the state and control vectors, x and u, into separate variables 
% for each component.
ani_posx = dataClass.RED_Px_m(startind:step:lastind);   % 1 x (k+1)
ani_posy = dataClass.RED_Py_m(startind:step:lastind);    % 1 x (k+1)
ani_theta = dataClass.RED_Rz_m(startind:step:lastind);   % 1 x (k+1)
ani_velx = dataClass.RED_Vx_mpers(startind:step:lastind);    % 1 x (k+1)
ani_vely = dataClass.RED_Vy_mpers(startind:step:lastind);    % 1 x (k+1)
ani_omega = dataClass.RED_RzD_radpers(startind:step:lastind)*180/pi; % converted to deg/s
ani_ux = dataClass.RED_Fx(startind:step:lastind);      % 1 x (k+1)
ani_uy = dataClass.RED_Fy(startind:step:lastind);      % 1 x (k+1)
ani_tau = dataClass.RED_Tz(startind:step:lastind);     % 1 x (k+1)

% Reshape the solution into an array
%reshapesol = reshape(sol,[9,length(ani_posx)]); % 9 x (k+1)

% Assign the trajectory and attitude solutions (over the horizon) for each 
% iteration into separate variables. This is xt.
%ani_targetx = soliter(1:9:end,:); %  (N+1) x k
%ani_targety = soliter(2:9:end,:); %  (N+1) x k
%ani_targetx(:,end+1) = NaN; % The final iteration is a docked chaser; (N+1) x (k+1)
%ani_targety(:,end+1) = NaN; % (N+1) x (k+1)

% Assign the paths of the obstacles into new variables
%ani_pobj1 = solpobj(:,1:2); % (k+1) x 2
%ani_pobj2 = solpobj(:,3:4); % (k+1) x 2

% Assign the docking mechanism lengths into separate variables
ani_rc = 0.15;        % m
rt = 0.15; % m

% Assign the parameters for the target spacecraft position, attitude, and
% velocities into separate variables
ani_ptar(:,1) = dataClass.BLACK_Px_m(startind:step:lastind); % (k+1) x 2
ani_ptar(:,2) = dataClass.BLACK_Py_m(startind:step:lastind);
ani_targett = dataClass.BLACK_Rz_m(startind:step:lastind);                % (k+1) x 1
ani_tar_velx = dataClass.BLACK_Vx_mpers(startind:step:lastind);   % 1 x (k+1); m/s
ani_tar_vely = dataClass.BLACK_Vy_mpers(startind:step:lastind);   % 1 x (k+1); m/s
ani_tar_omega = dataClass.BLACK_RzD_radpers(startind:step:lastind)*180/pi; % 1 x (k+1); deg/s

ani_pobj(:,1) = dataClass.BLUE_Px_m(startind:step:lastind); % (k+1) x 2
ani_pobj(:,2) = dataClass.BLUE_Py_m(startind:step:lastind);
ani_thetaobj(:,2) = dataClass.BLUE_Rz_m(startind:step:lastind);

% Assign the slopes of the entry cones to separate variables
%ani_m1 = solm1(:,1);    % (k+1) x 1
%ani_m2 = solm2(:,1);    % (k+1) x 1

% Assign the x,y positions of the docking port position to separate
% variables
%ani_pdockx = solpdock(:,1); % (k+1) x 1
%ani_pdocky = solpdock(:,2); % (k+1) x 1

% Calculate the total control force and total cost for the solution of k
% iterations. The total control force parameter (in N) can be compared
% across solutions as a performance metric.
%[com_totalcontrol,com_totalcost] = theComparator(ani_rc,ani_ux(1:end-1),ani_uy(1:end-1),ani_tau(1:end-1),solJ);

% Write a string for plotting purposes.
%ani_analysisstring = strcat('Total control force (N): ', string(com_totalcontrol),...
%    '; Total Cost: ', string(com_totalcost));

% Write a string for plotting purposes.
%ani_CPUstring = strcat('Solver time (s): ', string(min(solCPUtime)),...
%    ' (min), ',string(max(solCPUtime)),' (max), ',string(mean(solCPUtime)),' (mean)') ; 

% Write a string for plotting purposes.
%if ALGORITHM == 0
%    ani_algorithm = 'NMPC' ; 
%elseif ALGORITHM == 1
%    ani_algorithm = 'Compact NMPC' ; 
%elseif ALGORITHM == 2
%    ani_algorithm = 'Linear MPC' ; 
%elseif ALGORITHM == 3
%    ani_algorithm = 'Linear Compact MPC' ; 
%elseif ALGORITHM == 4
%    ani_algorithm = 'Linear Compact Fmincon' ; 
%elseif ALGORITHM == 5
%    ani_algorithm = 'Linear Compact Quadprog' ; 
%end

% ================================
% ========== PLOTTING ============
% ================================

% This gives the sizing of each of the components in the figure,
%  - trajectory (x vs y)
%  - velocity vs time
%  - angular velocity vs time
%  - control forces vs time
sp_nrows = 3; % Number of rows in total
sp_ncols = 5; % Number of columns in total
sp_wlin = 2;  % Number of columns for [] vs time plots

% Assigning the correct indexing to each component in the figure
sp_traj = [1:sp_nrows, 
           1+sp_ncols:sp_ncols+sp_nrows, 
           1+2*sp_ncols:2*sp_ncols+sp_nrows];
sp_vel = [sp_nrows+1:sp_nrows+sp_wlin];
sp_avel = [sp_ncols+sp_nrows+1:sp_ncols+sp_nrows+sp_wlin];
sp_ctrl = [sp_ncols*2+sp_nrows+1:sp_ncols*2+sp_nrows+sp_wlin];

% Start a figure
figure

% Looping through each of the iterations (with a nonzero control force,
% i.e. right up until docking)
for i = 1 : length(ani_posx)-1
    
    % ------- Trajectory Plotting ------- %
    subplot(sp_nrows,sp_ncols,sp_traj)
    plot(ani_posx(1:i), ani_posy(1:i), 'k.-'); % Plot the solution up until i
    hold on
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
    
    chaser = theChaserSpacecraft(ani_posx(i), ani_posy(i), ani_theta(i), ani_rc);
    patch(chaser(:,1), chaser(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'r') % Plot the chaser spacecraft
    hold on
    target = theTargetSpacecraft(ani_ptar(i,1), ani_ptar(i,2), ani_targett(i), rt);
    patch(target(:,1), target(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'k') % Plot the target spacecraft
    hold on
    obstac = theObstacleSpacecraft(ani_pobj(i,1), ani_pobj(i,2), ani_thetaobj(i), rt);
    patch(obstac(:,1), obstac(:,2), 'w', 'facealpha', 0.5, 'edgecolor', 'b') % Plot the target spacecraft
    
    grid on
    ax = gca;
    ax.FontSize = 6;
    %ylim([min([ani_posy.',ani_ptar(:,2).'])-ani_rc*2 max([ani_posy.',ani_ptar(:,2).'])+ani_rc*2])
    %xlim([min([ani_posx.',ani_ptar(:,1).'])-ani_rc*2 max([ani_posx.',ani_ptar(:,1).'])+ani_rc*2])
    xlim([0,3.5])
    ylim([0,2.4])
    
    %daspect([1 1 1])
    %axis square;
    
    %text(min([ani_posx,ani_ptar(:,1).'])-ani_rc*2, 1.06*(max([ani_posy,ani_ptar(:,2).'])+ani_rc*2), ani_CPUstring,'FontSize',8)
    %hold on
    %text(min([ani_posx,ani_ptar(:,1).'])-ani_rc*2, 1.02*(max([ani_posy,ani_ptar(:,2).'])+ani_rc*2), ani_analysisstring,'FontSize',8)
    %hold on
    %text(min([ani_posx,ani_ptar(:,1).'])-ani_rc*2 + 1.8*(max([ani_posx,ani_ptar(:,1).']) - min([ani_posx,ani_ptar(:,1).'])), 1.06*(max([ani_posy,ani_ptar(:,2).'])+ani_rc*2), ani_algorithm,'FontSize',8)
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    
    
    % ------- Velocity Plotting ------- %
    subplot(sp_nrows,sp_ncols,sp_vel)
    plot(1:i, ani_velx(1:i), '.-r','LineWidth', 0.5);
    hold on
    plot(1:i, ani_vely(1:i), '.-b','LineWidth', 0.5);
    hold on
    plot(1:i, ani_tar_velx(1:i), ':r','LineWidth', 1)
    hold on
    plot(1:i, ani_tar_vely(1:i), ':b','LineWidth', 1)
    %hold on
    %text(length(ani_posx)-35, 2.5*(max([ani_velx,ani_vely,ani_tar_velx,ani_tar_vely])), ani_algorithm,'FontSize',7)
    
    grid on
    ax = gca;
    ax.FontSize = 6;
    leg1 = legend('v_{x,c}','v_{y,c}','v_{x,t}','v_{y,t}','Location','northoutside', 'NumColumns',4);
    leg1.ItemTokenSize = [15,9];
    ylim([min([ani_velx;ani_vely;ani_tar_velx;ani_tar_vely]) max([ani_velx;ani_vely;ani_tar_velx;ani_tar_vely])])
    xlim([1 length(ani_posx)])
    xlabel('Iteration Number', 'FontSize', 5)
    ylabel('Velocity [m/s]')
    
    
    % ------- Angular Velocity Plotting ------- %
    subplot(sp_nrows,sp_ncols,sp_avel)
    plot(1:i, ani_omega(1:i), '.-k','LineWidth', 0.5);
    hold on
    plot(1:i, ani_tar_omega(1:i), ':k','LineWidth', 0.5);
    
    grid on
    ax = gca;
    ax.FontSize = 6;
    leg2 = legend('\omega_{c}','\omega_{t}','Location','northoutside', 'NumColumns',2);
    leg2.ItemTokenSize = [15,9];
    ylim([min([ani_omega]) max([ani_omega])])
    xlim([1 length(ani_posx)])
    xlabel('Iteration Number')
    ylabel('Angular Velocity [deg/s]')
    
    
    % ------- Control Force Plotting ------- %
    subplot(sp_nrows,sp_ncols,sp_ctrl)
    plot(1:i, ani_ux(1:i), '.-r','LineWidth', 0.5);
    hold on
    plot(1:i, ani_uy(1:i), '.-b','LineWidth', 0.5);
    hold on
    plot(1:i, ani_tau(1:i), '.-k','LineWidth', 0.5);
    
    grid on
    ax = gca;
    ax.FontSize = 6;
    leg3 = legend('u_{x}','u_{y}','\tau','Location', 'northoutside', 'NumColumns',3);
    leg3.ItemTokenSize = [15,9];
    ylim([min([ani_ux;ani_uy; ani_tau]) max([ani_ux;ani_uy; ani_tau])])
    xlim([1 length(ani_posx)])
    xlabel('Iteration Number')
    ylabel('Control Force [N],[Nm]')
    
    % Save the current figure as a "frame" for the animation
    movieVector(i) = getframe(gcf) ;
    
    if i == length(ani_posx)-1
        % last still
        saveas(gcf,strcat(strcat('Saved Data/',PRENAME),'/finalstill.png'))
%         saveas(gcf,strcat(strcat('Saved Data/',PRENAME),'/finalstill'))
    end
    
    clf
end

%%
% ================================
% ========= ANIMATING ============
% ================================

% Compile the figures from previous into an animation and save (note, only
% the .avi format was saving properly)
myWriter = VideoWriter(strcat(PRENAME,'mov'));
myWriter.FrameRate = (2/0.1/step) * 10;
open(myWriter);
for i=1:length(movieVector)
    frame = movieVector(i) ;    
    writeVideo(myWriter, frame);
end
close(myWriter);



  