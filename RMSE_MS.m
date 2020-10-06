function [uni_d, RMSE_M, RMSE_S] = RMSE_MS(d, RMSE)
%RMSE_MS This function calculates the mean and standard deviation of RMSE
% at each measurement length (d)
% Input parameters: 
        % d: measurement length
        % RMSE: calculated RMSE (as a function of d)
% Output parameters:    
        % uni_d: unique measurement length d
        % RMSE_M: mean value of all RMSE values at length d
        % RMSE_S: standard deviation of all RMSE values at length d
        
uni_d = unique(d);
L = length(uni_d);
RMSE_M = zeros(L,1);
RMSE_S = zeros(L,1);

for i = 1:1:L
    d_i=uni_d(i);
    pos_i = d==d_i;
    RMSE_i = RMSE(pos_i);
    RMSE_M(i) = mean(RMSE_i);
    RMSE_S(i) = std(RMSE_i,1);
end


end

