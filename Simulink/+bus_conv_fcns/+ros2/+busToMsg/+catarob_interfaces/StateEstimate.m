function rosmsgOut = StateEstimate(slBusIn, rosmsgOut)
%#codegen
%   Copyright 2021 The MathWorks, Inc.
    rosmsgOut.x = double(slBusIn.x);
    rosmsgOut.y = double(slBusIn.y);
    rosmsgOut.psi = double(slBusIn.psi);
    rosmsgOut.u = double(slBusIn.u);
    rosmsgOut.v = double(slBusIn.v);
    rosmsgOut.r = double(slBusIn.r);
    rosmsgOut.covariance = double(slBusIn.covariance);
end
