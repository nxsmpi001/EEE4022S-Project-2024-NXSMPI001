function [dx, y] = SystemDynamics(t, x, u, x_g,m11,m22,m33,d11,d22,d33,k, varargin)
%%   SystemDynamics
%    F = SystemDynamics(D11,D22,D33,K,M11,M22,M33,NU,PSI,R,U1,U2,UPSILON,W1,W2,X_G)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    15-Oct-2024 00:22:47

t2 = cos(x(3));
t3 = sin(x(3));
t4 = x_g.^2;
t5 = 1.0./m11;
t6 = 1.0./m22;
t7 = 1.0./m33;

% Compute derivative of states
dx = [-x(5).*t3+t2.*x(4);x(5).*t2+t3.*x(4);x(6);t2.*u(3)+t3.*u(4)-d11.*t5.*x(4)+k.*t5.*u(1)+k.*t5.*u(2)+m22.*x(5).*x(6).*t5+m22.*x(6).^2.*t5.*x_g;t2.*u(4).*1.0-t3.*u(3).*1.0-d22.*x(5).*t6-d22.*x(5).*t4.*t7+d33.*x(6).*t7.*x_g-m11.*x(6).*t6.*x(4)-k.*t7.*u(1).*x_g.*4.15e-1+k.*t7.*u(2).*x_g.*4.15e-1-m11.*x(6).*t4.*t7.*x(4)+m22.*x(6).*t4.*t7.*x(4)-m11.*x(5).*t7.*x(4).*x_g+m22.*x(5).*t7.*x(4).*x_g+m22.*t2.*t4.*t7.*u(4).*1.0-m22.*t3.*t4.*t7.*u(3).*1.0;d33.*x(6).*t7.*-1.0+k.*t7.*u(1).*4.15e-1-k.*t7.*u(2).*4.15e-1+d22.*x(5).*t7.*x_g.*1.0+m11.*x(5).*t7.*x(4).*1.0-m22.*x(5).*t7.*x(4).*1.0+m11.*x(6).*t7.*x(4).*x_g.*1.0-m22.*x(6).*t7.*x(4).*x_g.*1.0-m22.*t2.*t7.*u(4).*x_g.*1.0+m22.*t3.*t7.*u(3).*x_g.*1.0];

% Output Function
y = [x(1);x(2);x(3);x(6)];
end
