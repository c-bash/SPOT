% ---------------- The Bellman's Principle MPC Style -------------------

bel_targetx = soliter(1:9:end,:);
bel_targety = soliter(2:9:end,:);
bel_targettheta = soliter(3:9:end,:);
bel_targetvx = soliter(4:9:end,:);
bel_targetvy = soliter(5:9:end,:);
bel_targetomega = soliter(6:9:end,:);

bel_t = 0:2:2*(length(bel_targetx)+N+1)-1;

% -------------- plot -----------------
ms = 5;
figure('Position',[300,150,700,400])
for i = 1:length(bel_targetx)
    subplot(3,2,1)
    plot(bel_t(i:i+N),bel_targetx(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targetx(1,i),'r.','MarkerSize',ms)
    grid on
    set(gca, 'XTickLabel', [])
    ylabel('X-position [m]')

    subplot(3,2,3)
    plot(bel_t(i:i+N),bel_targety(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targety(1,i),'r.','MarkerSize',ms)
    grid on
    set(gca, 'XTickLabel', [])
    ylabel('Y-position [m]')

    subplot(3,2,5)
    plot(bel_t(i:i+N),bel_targettheta(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targettheta(1,i),'r.','MarkerSize',ms)
    grid on
    ylabel('Attitude [rads]')
    xlabel('Time [s]')

    subplot(3,2,2)
    plot(bel_t(i:i+N),bel_targetvx(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targetvx(1,i),'r.','MarkerSize',ms)
    grid on
    set(gca, 'XTickLabel', [])
    ylabel('X-velocity [m]')

    subplot(3,2,4)
    plot(bel_t(i:i+N),bel_targetvy(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targetvy(1,i),'r.','MarkerSize',ms)
    grid on
    set(gca, 'XTickLabel', [])
    ylabel('Y-velocity [m/s]')

    subplot(3,2,6)
    plot(bel_t(i:i+N),bel_targetomega(:,i),'k','LineWidth',0.5)
    hold on
    plot(bel_t(i),bel_targetomega(1,i),'r.','MarkerSize',ms)
    grid on
    ylabel('Angular velocity [deg/s]')
    xlabel('Time [s]')
end

savefig(strcat(PRENAME, 'bellmans.fig'))