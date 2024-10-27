function [G, Gmv, Ge] = costFcnJacobian(X, U, ~, data, Q, R, Qt)
%% Jacobian of Cost Function for computational efficiency

% Extract prediction horizon
p = data.PredictionHorizon;

% Control inputs
u1 = U(1:p-1, data.MVIndex(1));
u2 = U(1:p-1, data.MVIndex(2));

% Output variables
x = X(2:p+1, 1);
y = X(2:p+1, 2);

% Error vector
err_x = data.References(1, 1) - x;
err_y = data.References(1, 2) - y;

% Initialise gradients
G = zeros(p,6);
Gmv = zeros(p,length(data.MVIndex));
Ge = 0;

% Compute partial derivatives for x, y
for k = 1:p-1
    G(k,1) = -2 * Q * err_x(k);
    G(k,2) = -2 * Q * err_y(k);
end
G(p,1) = -2 * Qt * err_x(p);
G(p,2) = -2 * Qt* err_y(p);

% Compute partial derivatives for u1, u2
for k = 1:p-1
    Gmv(k,1) = 2 * R * u1(k);
    Gmv(k,2) = 2 * R * u2(k);
end

end
