function C = outputFcnJacobian(~, ~, ~, ~, ~)
%% Jacobian of the Output function for computational efficiency 

% Initialise gradient
C = zeros(4,6); 

% Compute partial derivatives for output - C(i,j) = dy_i/dx_j
C(1,1) = 1;
C(2,2) = 1;
C(3,3) = 1;
C(4,6) = 1;
end