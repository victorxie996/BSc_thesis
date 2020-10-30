[f_i_new, g_i_new] = amplfunc(x_i_temp,0);
cv_old = 0;
cv_new = 0;


% evaluate constraint violations
cv_old = sum(abs(max(g_i-cu,0))) + sum(abs(min(g_i-cl,0)));
cv_new = sum(abs(max(g_i_new-cu,0))) + sum(abs(min(g_i_new-cl,0)));

% test if the new constraint violation is smaller.
constraint_ratio = (cv_old - cv_new)/cv_old;
% test if objective value is better (i.e. if new f is lower than old f)
objective_ratio = (f_i - f_i_new)./(f_i-predicted_obj);

%if (cv_old > 1e-6)
if (cv_old > 1e-12)
if min(constraint_ratio, objective_ratio) > 0.75
    trust = trust*2;
    x_i = x_i_temp;
elseif min(constraint_ratio, objective_ratio) < 0.05
            trust = trust/2;
else
    x_i = x_i_temp;
end
else
  if objective_ratio>0.75 
    trust = trust*2;
    x_i = x_i_temp;
elseif objective_ratio < 0.05
            trust = trust/2;
  else
    x_i = x_i_temp;
  end
end
trustregion = trust*ones(length(bl),1);
    lbound = -min(trustregion, abs(bl-x_i));
    ubound = min(trustregion, abs(bu-x_i));
