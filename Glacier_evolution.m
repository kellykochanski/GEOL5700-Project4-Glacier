%% 1D Valley Glacier
% Kelly Kochanski 23-02-2016, for GEOL 5700-004
%   (worked with Marko Visnjic)

%%%%% DO GLACIERS RECORD A CHANGING CLIMATE? %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This simulation investigates a glacier's response to an oscillating ELA.
% Specifically, it shows that glaciers do not respond to all signals
% equally.
%
% The most obvious part of a glacier's history is it's terminus. This is
% marked by moraines, ponds, suffocated lichen... 
% This code simulates glacier movement with an oscillating ELA, and tracks
% the position of the terminus as a function of time.
%
% The period of the oscillation is controlled by parameter 'omega'.
% I recommend you run the program with at least two values of omega:
%  -  omega = 0.0006 (ten-thousand-year period, ie Laurentide Ice Sheet)
%  -  omega = 0.006 (one-thousand-year period)
%  -  omega = 0.01 (six-hundred-year period, ie little ice age,
%  post-industrial warming)
%
% You'll see that the glacier preserves these signals very differently.
%
% If you cut the program off by hitting cntrl-C, you will not see all the
% plots. The ten-thousand year period cannot run fast and be stable, sorry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters
clf
figure (1)
rho_i = 917; %Density of Ice
gamma = 0.01;   %variation of ablation with altitude, y^-1
s = 0.4;        %slope, m/m
g = 9.81;%*(pi*10^7)^2;       %gravity, m/yr^2
A = 2.1*10^-16; %ratio between glacier stress and strain

%% Arrays
dx = 1;
xmax = 1000;
x = 0:dx:xmax; %Distance Array m

dz = 0.1;
zmax = 100;
zb =zmax:-dz:0; %Elevation Array m

%% Climate variation
ELA0 = (2/3)*zmax;
omega = 0.01; %1/years
ELA = @(t) 15.*sin(2*pi*omega*t) + ELA0;

b0 = gamma*(zb-ELA0);  % ablation of glacier surface, m/year
b = b0;
G0 = 30+30.*(-(x-(xmax/3)).^2)/(xmax/3).^2; %Initial Glacial thickness m
G0 = max(0, G0);
G = G0.*ones(size(x)); %Initial Glacial Array

%Elevation of glacier surface
h = (zmax-s*x)+G; %m

%Time
dt = 1/52; %yrs
tmax = 1.5*pi/omega;
t = 0:dt:tmax;
imax = length(t);


%% Solver and time advancement

%% 1. Run with constant climate for 300 years (approx time to steady state)
disp('Running simulation to steady state, please wait...')
for i =1:300/dt;
    t = i*dt;
    b = gamma*(zb-ELA0);
   % array of fluxes of ice
   q = A*(rho_i*g*s).^3.*G.^5/5;
   q = [0 q];
   dqdx = diff(q)/dx;
   dGdt = b - dqdx;
   % update glacier thickness and elevation
   G = G+dGdt*dt;
   G = max(G, 0);
   h = G + zb;
end

disp('Steady state reached. Running simulation...')

%% 2. Start varying climate and watch the glacier respond
tvec = zeros(1,round(imax/400)-1);          
terminus_vec = tvec; % store terminus position as fn of time
for i =1:imax
    t = i*dt;
    b = gamma*(zb-ELA(t));
   % array of fluxes of ice
   q = A*(rho_i*g*s).^3.*G.^5/5;
   q = [0 q];
   dqdx = diff(q)/dx;
   dGdt = b - dqdx;
   % update glacier thickness and elevation
   G = G+dGdt*dt;
   G = max(G, 0);
   h = G + zb;
   % plot all the things every 10 steps
   % Keep track of glacier extend as a function of time
   if round(mod(i, 400)) == 0;
       % End of glacier:
       no_glacier     = find(G==0); %index of glacier terminus
       terminus_vec(i/400)= zb(no_glacier(1));
       tvec(i/400)    = t;
        plot(x,h,'b')
        hold on
        plot(x,zb,'k')
        plot([0,xmax], [ELA(t), ELA(t)], ':r')
        hold off
        axis([0 x(end) 0 150])
        title(sprintf('Time: %d years', round(t)))
        xlabel('Range, m'); ylabel('Elevation, m')
        drawnow
   end
end

figure;
title('How does glacier length respond to changes in ELA?')
plot(tvec, terminus_vec, 'b', 'linewidth', 2)
hold on
plot(tvec, ELA(tvec), 'r', 'linewidth', 2)
xlabel('Time, years'); ylabel('Elevation, m')
legend('Elevation of glacier terminus', 'ELA')