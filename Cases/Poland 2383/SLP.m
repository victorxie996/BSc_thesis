function [x_i_temp, g_i, f_i, predicted_obj,d,lam] = SLP(x_i,cu,cl,lbound, ubound, v)
% This function performs one iteration of SLP.

% perform iteration of SLP
[f_i, g_i] = amplfunc(x_i,0);
[nabla_f_i, nabla_g_i] = amplfunc(x_i,1);
W = amplfunc(-v);
A = [ nabla_g_i;-nabla_g_i];
b = [cu-g_i;-cl+g_i];
%remove_inf;

infx = (b<inf);
A2 = A(infx,:); b2 = b(infx);

options = optimset('linprog');
options.Display = 'off';
[d, fval, flag, out, lam2] = linprog(nabla_f_i.',A2,b2,[],[],lbound, ubound, options);

flag;
lam = zeros(size(b));
m = length(b);
lamtmp(infx) = lam2.ineqlin';
lam = lamtmp(1:m/2)+lamtmp(m/2+1:m);
lam = lam';

predicted_obj = f_i+nabla_f_i'*d;

% define temporary x_i value to test progress before taking step
x_i_temp = d + x_i;
end

