function [x_i_temp, g_i, f_i, predicted_obj,d] = SLP(x_i,cu,cl,lbound, ubound)
% This function performs one iteration of SLP.

% perform iteration of SLP
[f_i, g_i] = amplfunc(x_i,0);
[nabla_f_i, nabla_g_i] = amplfunc(x_i,1);
A = [ nabla_g_i;-nabla_g_i];
b = [cu-g_i;-cl+g_i];
%remove_inf;

infx = (b<inf);
A2 = A(infx,:); b2 = b(infx);

options = optimset('linprog');
options.Display = 'off';
[d, fval, flag, out, lam2] = linprog(nabla_f_i.',A2,b2,[],[],lbound, ubound, options);

predicted_obj = f_i+nabla_f_i'*d;

% define temporary x_i value to test progress before taking step
x_i_temp = d + x_i;
end

