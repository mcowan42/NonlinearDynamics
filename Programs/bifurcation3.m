clear all;
tic;
function xout = rk4(x,t,tau,derivsRK,M1,M2,rM1,rM2)
%  Runge-Kutta integrator (4th order)
% Input arguments -
%   x = current value of dependent variable
%   t = independent variable (usually time)
%   tau = step size (usually timestep)
%   derivsRK = right hand side of the ODE; derivsRK is the
%             name of the function which returns dx/dt
%             Calling format derivsRK(x,t,param).
%   param = extra parameters passed to derivsRK
% Output arguments -
%   xout = new value of x after a step of size tau
half_tau = 0.5*tau;
F1 = feval(derivsRK,x,t,M1,M2,rM1,rM2);  
t_half = t + half_tau;
xtemp = x + half_tau*F1;
F2 = feval(derivsRK,xtemp,t_half,M1,M2,rM1,rM2);  
xtemp = x + half_tau*F2;
F3 = feval(derivsRK,xtemp,t_half,M1,M2,rM1,rM2);
t_full = t + tau;
xtemp = x + tau*F3;
F4 = feval(derivsRK,xtemp,t_full,M1,M2,rM1,rM2);
xout = x + tau/6.*(F1 + F4 + 2.*(F2+F3));
endfunction;

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
endfunction;

%MAIN PROGRAM
c=20; %initial parameter value
ccount=1; %count the number of different values of c
plotcount=1;

while c<40
  mu = c/(1+c);
  M1 = mu;
  M2 = 1-mu;
  %Positions of the two primary masses
  rM1 = [-(1-mu),0];
  rM2 = [mu,0];
  %Set initial position and velocity of the object
  r0 = [.5*(M1-M2),sqrt(3)/2]; %initial position is the L4 Lagrange point 
  r = r0;
  v = [.01,.01];
  state = [ r(1) r(2) v(1) v(2) ];   % Used by R-K routines

  %Initialize time
  time = 0;
  %Loop over desired number of steps
  tau = 0.01; 
  nStep = 1e6;
  for istep=1:nStep
    %* Calculate new position and velocity using RK4.
    state = rk4(state,time,tau,'gravrk',M1,M2,rM1,rM2);
    r = [state(1) state(2)];
    v = [state(3) state(4)];
    time = time + tau;
    
    %Save the last point for plotting
    if (nStep==istep)
      xplot(plotcount)=r(1)-r0(1);
      yplot(plotcount)=r(2)-r0(2);
      rplot(plotcount)=norm(r-r0);
      cplot(plotcount)=c;
      plotcount=plotcount+1;
    end
  end
  
  %change the parameter value
  if c>=24.5 && c<25.5
    c=c+.01;  %small increment near the bifurcation point
  else
    c=c+.5;   %larger increments further from the bifurcation point
  end
  ccount=ccount+1;
end

figure(1); clf;
plot(cplot,xplot,'.');
xlabel('M1/M2'); ylabel('x position');

figure(2); clf;
plot(cplot,yplot,'.');
xlabel('M1/M2'); ylabel('y position');

figure(3); clf;
plot(cplot,rplot,'.');
xlabel('M1/M2'); ylabel('radial distance from L4');

disp(toc);