function deriv = gravrk(s,t,M1,M2,rM1,rM2)
%  Returns right-hand side of Kepler ODE; used by Runge-Kutta routines
%  Inputs
%    s              State vector [r(1) r(2) v(1) v(2)]
%    t              Time (not used)
%    M1,M2,rM1,rM2  Parameters related to the primary masses
%  Output
%    deriv  Derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]

%* Compute acceleration
r = [s(1) s(2)];  % Unravel the vector s into position and velocity
v = [s(3) s(4)];
r1=r-rM1;
r2=r-rM2;
accel = -M1*r1/norm(r1)^3-M2*r2/norm(r2)^3-2*[-v(2),v(1)]+r;

%* Return derivatives [dr(1)/dt dr(2)/dt dv(1)/dt dv(2)/dt]
deriv = [v(1) v(2) accel(1) accel(2)];
return;