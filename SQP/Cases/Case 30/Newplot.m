figure('Renderer', 'painters', 'Position', [10 10 1000 5000])
fig = gcf; 
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
set(0, 'DefaultAxesFontSize',15)
subplot(3,2,1); x = iteration(:,1); y = iteration(:,2);
    %+iteration(:,5)+iteration(:,6); 
    plot(x,y,'-.x');
    title('progress of trust region for case 9', 'interpreter', 'latex'); 
    axis([0 15 0 10]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('trust region', 'interpreter', 'latex');

    subplot(3,2,2); x = iteration1(:,1); y = iteration1(:,2);
    %+iteration(:,5)+iteration(:,6); 
plot(x,y,'-.x');
    title('progress of trust region for case 30', 'interpreter', 'latex'); 
    %axis([1 3 0 10]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('trust region', 'interpreter', 'latex');

    subplot(3,2,3); x = iteration(:,1); y = iteration(:,3);
    %+iteration(:,5)+iteration(:,6); 
  plot(x,y,'-.x');
    title('progress of objective value for case 9', 'interpreter', 'latex'); 
    axis([0 15 1000 5500]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('objective value', 'interpreter', 'latex');
    
    % case 30
    subplot(3,2,4); x = iteration1(:,1); y = iteration1(:,5);
    %+iteration(:,5)+iteration(:,6); 
plot(x,y,'-.x');
    title('progress of objective value for case 30', 'interpreter', 'latex'); 
    axis([1 3 560 580]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('objective value', 'interpreter', 'latex');

    subplot(3,2,5); x = iteration(:,1); y = iteration(:,6);
    %+iteration(:,5)+iteration(:,6); 
plot(x,y,'-.x');
    title('progress of constraint violation for case 9', 'interpreter', 'latex'); 
    axis([0 15 0 6]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('constraint violation', 'interpreter', 'latex');

    subplot(3,2,6); x = iteration1(:,1); y = iteration1(:,6);
    %+iteration(:,5)+iteration(:,6); 
  plot(x,y,'-.x');
    title('progress of constraint violation for case 30', 'interpreter', 'latex'); 
    %axis([1 3 0 2.5]);
    xlabel('iterations', 'interpreter', 'latex');
    ylabel('constraint violation', 'interpreter', 'latex');