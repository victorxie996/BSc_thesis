clear all;

%load data from AMPL Model
[x,bl,bu,v,cl,cu] = amplfunc('case118.nl');

% define inital point for SLP:
x_i = zeros(length(bl),1);
x_i(1:118)= ones(118,1);

% define inital size of trustregion
trust = 5;

%define lower and upper bounds for variables d
lbound = -min(trust, abs(bl-x_i));
ubound = min(trust, abs(bu-x_i));

% define an inital large current constraint violation
d=1;
lam = zeros(size(cl));


% set iteration counter to 1
i = 1;
fid = fopen('iteration.txt','w');
fprintf(fid,'Iteration & Trust region &Current Objective & predicted new objective & actual new objective &Current Constraint violation & new constraint violation & Constraint ratio & objective ratio \\\\ \n');

tic
while norm(d,inf) > 1e-5
    
    % find the next trial point x_i_temp
    [x_i_temp, g_i, f_i,predicted_obj,d,lam] = SQP(x_i,cu,cl,lbound,ubound,lam);
    
    trust_old = trust;
    % adjust the trust region to fit the potential changes made during the
    % testing conditions.
    lbound = -min(trust, abs(bl-x_i));
    ubound = min(trust, abs(bu-x_i));
    % test progress of temporary x_i value
    progress_test_SQP;
    
    % write a .txt file with all the iteration information.
    fprintf(fid,'%3.2d',i);
    fprintf(fid,' & ');
    fprintf(fid,'%3.8f',trust_old);
    fprintf(fid,' & ');
    fprintf(fid,'%5.2f',f_i);
    fprintf(fid,' & ');
    fprintf(fid,'%5.2f',predicted_obj);
    fprintf(fid,' & ');
    fprintf(fid,'%5.2f',f_i_new);
    fprintf(fid,' & ');
    fprintf(fid,'%3.8e',cv_old);
    fprintf(fid,' & ');
    fprintf(fid,'%3.8e',cv_new);
    fprintf(fid,' & ');
    fprintf(fid,'%3.6f',constraint_ratio);
    fprintf(fid,' & ');
    fprintf(fid,'%3.6f',objective_ratio);
    fprintf(fid,' \\\\ \n');
    % increase the iteration number by one
    i = i + 1;
    
end
t=toc;
fprintf(fid,'Time = %3.12f',t);
fclose(fid);
% evaluate the objective and the constraints at the solution found by SLP
[f_i, g_i] = amplfunc(x_i,0);
[nabla_f_i, nabla_g_i] = amplfunc(x_i,1);

% evaluate and define the Coefficient matrix of the constraints and the
% upper bound vector b
A = [ nabla_g_i;-nabla_g_i];
b = [cu-g_i;-cl+g_i];
W = amplfunc(-v);
%remove infinity values
infx = (b<inf);
A2 = A(infx,:);
b2 = b(infx);
% evaluate the Lagrangian multipliers
[d,fval,exitflag,output,lambda] = quadprog(W, nabla_f_i.',A2,b2,[],[],lbound,...
    ubound, zeros(size(x_i)));

% evaluate the testing condition for Lagrangian duality
dual_test = (nabla_f_i) + (lambda.ineqlin' *A2)' - lambda.lower + lambda.upper;

%print a table with the testing condtion, the lagrangian multipliers and
%the bounds of the problem.
table(dual_test, lambda.lower,lambda.upper, x_i-bl, bu-x_i, lbound, ubound,...
    'VariableNames',{'dual_test','lambdalower','lambdaupper',...
    'x_iminusbl', 'buminusx_i','lbound','ubound'})
