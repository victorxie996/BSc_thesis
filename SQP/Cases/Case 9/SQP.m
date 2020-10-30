function [x_i_temp, g_i, f_i, predicted_obj,d,lam] = SQP(x_i,cu,cl,lbound, ubound, v)
% This function performs one iteration of SLP.

% evaluate the objective and the constraints at the current trial point
[f_i, g_i] = amplfunc(x_i,0);
[nabla_f_i, nabla_g_i] = amplfunc(x_i,1);
W = amplfunc(-v);
% evaluate and define the Coefficient matrix of the constraints and the
% upper bound vector b
A = [ nabla_g_i;-nabla_g_i];
b = [cu-g_i;-cl+g_i];

%remove infinity values
infx = (b<inf);
A2 = A(infx,:); 
b2 = b(infx);

%change settings of quadprog so no outputs are diplayed
options = optimset('quadprog');
options.Display = 'off';

% solve the SQP
[d, fval, exitflag, output, lambda] = quadprog(W, nabla_f_i.',A2,b2,[],[],lbound,...
    ubound, zeros(size(x_i)), options);

% define the new dual variables at d, that we will need to evaluate the
% Hessian W in the next iteration.
lam = zeros(size(b));
m = length(b);
lamtmp(infx) = lambda.ineqlin';
lam = lamtmp(1:m/2)+lamtmp(m/2+1:m);
lam = lam';

% define the predicted objective at x_i_temp
predicted_obj = f_i+nabla_f_i'*d + 0.5*d'*W*d;

% define temporary x_i value to test progress before taking step
x_i_temp = d + x_i;
end

