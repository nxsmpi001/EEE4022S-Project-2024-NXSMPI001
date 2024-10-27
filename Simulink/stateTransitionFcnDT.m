function xk1 = stateTransitionFcnDT(xk, u)
%% Discrete State Transition Function
uk = u(1:2);  % Control input
dt = u(3);    % Time step

% Compute RK4 intermediate steps
k1 = stateTransitionFcn(xk, [uk; 0; 0]);
k2 = stateTransitionFcn(xk + 0.5 * dt * k1, [uk; 0; 0]);
k3 = stateTransitionFcn(xk + 0.5 * dt * k2, [uk; 0; 0]);
k4 = stateTransitionFcn(xk + dt * k3, [uk; 0; 0]);

% Compute next state using RK4 formula
xk1 = xk + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
end
