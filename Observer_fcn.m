function [x_hat_dot] = Observer_fcn( i , A, gamma, aij ,H_i ,M_i ,L_i ,x_hat_i ,x_hat_ij , y_i )

% i = which node we working on
% x_hat_ij = exepct a matrix [x1 x2 x3 x4]
% ------
% x_hat_i = estimate value of the node we working on

% Initialize commu

% commu = zeros(size(x_hat_i));
commu = 0;


% Calculate commu for each adjacency system
for j = 1:size(aij, 2)
    if j ~= i
        commu = commu + aij(i, j) * (x_hat_ij(j) - x_hat_i);
    end
end


x_hat_dot = A*x_hat_i + L_i*(y_i - H_i*x_hat_i) + gamma*M_i*commu;


end

