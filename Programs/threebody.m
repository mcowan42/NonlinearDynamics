clear all;

%Masses of the two primary masses
c = input('Enter the ratio of the primary masses M1/M2: ');
mu = c/(1+c);
M1 = mu;
M2 = 1-mu;
%Positions of the two primary masses
rM1 = [-(1-mu),0];
rM2 = [mu,0];
%Set initial position and velocity of the object
r = [.5*(c-1)/(c+1),sqrt(3)/2];  
v = [.01,.01];
state = [ r(1) r(2) v(1) v(2) ];   % Used by R-K routines
%Initialize time
time = 0;
j=1;
k=1;
xtemp(1)=r(1);
ytemp(1)=r(1);

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
%Loop over desired number of steps
tau = .01; 
nStep = input('Enter the number of steps to calculate: ');
for iStep=1:nStep  
  %* Record position for plotting
  xplot(iStep) = r(1);
  yplot(iStep) = r(2);
  tplot(iStep) = time;
  %Poincare section - plot position once per period of primaries
  if (mod(time,2*pi)<tau)
    Pxplot(j)=r(1);
    Pyplot(j)=r(2);
    j=j+1;
  end
  
  if (mod(iStep,nStep/200)<1)
    k=k+1;
    xtemp(k)=r(1);
    ytemp(k)=r(2);
    figure(1);
    %plot(xtemp,ytemp,'.',rM1(1),rM1(2),'+',rM2(1),rM2(2),'+');
    %drawnow;
  endif

  %* Calculate new position and velocity using RK4.
  state = rk4(state,time,tau,'gravrk',M1,M2,rM1,rM2);
  r = [state(1) state(2)];
  v = [state(3) state(4)];
  time = time + tau;
  
end

figure(2);
plot(xplot,yplot,'.',rM1(1),rM1(2),'+',rM2(1),rM2(2),'+',xplot(1),yplot(1),'o');
xlabel('x position'); ylabel('y position');
title('Particle Trajectory');

figure(3); clf;
plot(Pxplot,Pyplot,'.');
title('Poincare section');

figure(4);
plot(xplot,yplot,'.',rM1(1),rM1(2),'+',rM2(1),rM2(2),'+',xplot(1),yplot(1),'o');
axis([-2,2,-2,2]);
xlabel('x position'); ylabel('y position');
title('Particle Trajectory');