clear all;

%load data from AMPL Model
[x,bl,bu,v,cl,cu] = amplfunc('case9.nl');

% define inital point for SLP:
x_i = zeros(length(bl),1);
x_i(1:9)= ones(9,1);

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

while norm(d,inf) > 1e-6
    
    % find the next trial point x_i_temp
    [x_i_temp, g_i, f_i,predicted_obj,d,lam] = SQP(x_i,cu,cl,lbound,ubound,lam);
    
    trust_old = trust;
    % test progress of temporary x_i value
    progress_test_SQP;
    
    % display the important values at the current iteration
    disp(sprintf('%4d %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g\n', ...
        i, trust_old, f_i, predicted_obj, f_i_new, cv_old, cv_new, ...
        constraint_ratio, objective_ratio));
    
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
    fprintf(fid,'%3.8f',cv_old);
    fprintf(fid,' & ');
    fprintf(fid,'%3.8f',cv_new);
    fprintf(fid,' & ');
    fprintf(fid,'%3.6f',constraint_ratio);
    fprintf(fid,' & ');
    fprintf(fid,'%3.6f',objective_ratio);
    fprintf(fid,' \\\\ \n');
    % increase the iteration number by one
    i = i + 1;
    
end
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

fid = fopen('dualtest.txt','w');
fprintf(fid,'Entry &Dual Test &Entry &Dual Test&Entry &Dual Test \\\\ \n');
for i = 1:length(dual_test)/3
    fprintf(fid,'%2.0f',i);
    fprintf(fid,' & ');
    fprintf(fid,'%3.3e',dual_test(i));
    fprintf(fid,' & ');
    fprintf(fid,'%2.0f',i+length(dual_test)/3);
    fprintf(fid,' & ');
    fprintf(fid,'%3.3e',dual_test(i+length(dual_test)/3));
    fprintf(fid,' & ');
    fprintf(fid,'%2.0f',i+2*length(dual_test)/3);
    fprintf(fid,' & ');
    fprintf(fid,'%3.3e',dual_test(i+2*length(dual_test)/3));
    fprintf(fid,' \\\\ \n');
end
fclose(fid);

fid = fopen('lowerbound.txt','w');
fprintf(fid,'Entry & k & $x_i - bl$& lbound&Entry & k & $x_i - bl$& lbound \\\\ \n');
for i = 1:length(dual_test)/2
    fprintf(fid,'%2.0f',i);
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lambda.lower(i));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',x_i(i)-bl(i));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lbound(i));
    fprintf(fid,' & ');
    fprintf(fid,'%2.0f',i+length(dual_test)/2);
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lambda.lower(i+length(dual_test)/2));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',x_i(i+length(dual_test)/2)-bl(i+length(dual_test)/2));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lbound(i+length(dual_test)/2));
    fprintf(fid,' \\\\ \n');
end
fclose(fid);
fid = fopen('upperbound.txt','w');
fprintf(fid,'Entry & l &  $bu - x_i$  & ubound &Entry & l &  $bu - x_i$  & ubound \\\\ \n');
for i = 1:length(dual_test)/2
    fprintf(fid,'%2.0f',i);
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lambda.upper(i));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',bu(i)-x_i(i));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',ubound(i));
    fprintf(fid,' & ');
    fprintf(fid,'%2.0f',i+length(dual_test)/2);
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',lambda.upper(i+length(dual_test)/2));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',bu(i+length(dual_test)/2)-x_i(i+length(dual_test)/2));
    fprintf(fid,' & ');
    fprintf(fid,'%3.2e',ubound(i+length(dual_test)/2));
    fprintf(fid,' \\\\ \n');
end
fclose(fid);
