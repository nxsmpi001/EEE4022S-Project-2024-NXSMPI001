function [y, bounds] = magMeasurementFcn(x)
%% Magnetometer Measurement Function
y = x(3);           % Heading is measured
bounds = [0 2*pi];  % Heading is bounded between 0 and 360 degrees
end