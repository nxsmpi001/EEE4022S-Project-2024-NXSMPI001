function J = costFcn(X,U,~,data,Q, R, Qt)
%% Cost Function

% Extract prediction horizon
p = data.PredictionHorizon;

% Control inputs
u1 = U(1:p-1,data.MVIndex(1));
u2 = U(1:p-1,data.MVIndex(2));

% State variables
x = X(2:p+1,1);
y = X(2:p+1,2);

% Output variables
err = [data.References(1,1)-x, data.References(1,2)-y];

% Compute cost function
J = Q*sum(err(1:p-1,1).^2 + err(1:p-1,2).^2) + R*sum(u1.^2+u2.^2) ...
    + Qt*(err(p,1)^2 + err(p,2)^2); 
end