% evaluate the objective function and the constraints at the new trial
% point x_i_temp
[f_i_new, g_i_new] = amplfunc(x_i_temp,0);


% evaluate constraint violations
cv_old = sum(abs(max(g_i-cu,0))) + sum(abs(min(g_i-cl,0)));
cv_new = sum(abs(max(g_i_new-cu,0))) + sum(abs(min(g_i_new-cl,0)));

% define the constraint ratio
constraint_ratio = (cv_old - cv_new)/cv_old;

% define the objective ratio
objective_ratio = (f_i - f_i_new)./(f_i-predicted_obj);

% start testing conditions

if cv_old > 1e-6
    % if both ratios are greater than 0.75 we take the step and increase the
    % trust region
    if min(constraint_ratio, objective_ratio) > 0.75
        trust = trust*2;
        x_i = x_i_temp;
        
        % if either of the ratios is lower than 0.05 we reduce the trust region and
        % reject the new trial point
    elseif min(constraint_ratio, objective_ratio) < 0.05
        trust = trust/2;
        
        % if both constraints are within the range of 0.05 to 0.75 we take the step
        % without cahnging the trust region.
    else
        x_i = x_i_temp;
    end
else
    % if both ratios are greater than 0.75 we take the step and increase the
    % trust region
    if objective_ratio > 0.75
        trust = trust*2;
        x_i = x_i_temp;
        
        % if either of the ratios is lower than 0.05 we reduce the trust region and
        % reject the new trial point
    elseif objective_ratio < 0.05
        trust = trust/2;
        
        % if both constraints are within the range of 0.05 to 0.75 we take the step
        % without cahnging the trust region.
    else
        x_i = x_i_temp;
    end
end