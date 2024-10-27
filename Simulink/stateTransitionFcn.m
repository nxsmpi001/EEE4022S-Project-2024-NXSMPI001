function dxdt = stateTransitionFcn(x, u, ~, ~, ~)
%% auto-generated state function of nonlinear grey box
%# codegen
% Load the parameters
p = {0.699999982041001, 43.6722524126937, 63.0144902620752, 29.999999915303, 35.3670238623635, 74.8832173268989, 125.558050313511, 3.46823335263162};
dxdt = SystemDynamics(0, x, u, p{:});
end
