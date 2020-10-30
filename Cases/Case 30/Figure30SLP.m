load('iteration.mat')


subplot(3,1,1);
plot(iteration(:,1),iteration(:,2),'b')
xlabel('iterations'),ylabel('trust region'), title('progress of trust region')

% the constraint ratio against iteration plot does not mean alot.
% plot(iteration(:,1),iteration(:,8))
% xlabel('iterations'),ylabel('constraint ratio'), title('case 9')
% figure;
subplot(3,1,2);
plot(iteration(:,1),iteration(:,5),'r')
xlabel('iterations'),ylabel('objective value'), title('progress of objective value')
subplot(3,1,3);
plot(iteration(:,1),iteration(:,6),'r')
xlabel('iterations'),ylabel('constraint violation'), title('progress of constraint violation'), ylim([0 5])