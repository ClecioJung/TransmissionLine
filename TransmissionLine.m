%---------------------------------------------
% Transmission line simulation (coaxial cable)
%---------------------------------------------

clc;
close all;

%---------------------------------------------
% Parameters
%---------------------------------------------

% Transmission line parameters
R = 45.3e-3;        % ohm/m
G = 0;              % S/m (information not available in the catalog)
L = 0.205e-6;       % H/m
C = 82e-12;         % F/m
Z = sqrt(L/C);      % ohm
vel = 1/sqrt(L*C);  % m/s

% Load resistor
Rl = 10;            % ohm

% Arbitrary function generator
Rg = 50;            % ohm
Ton = 10e-9;        % s
Toff = 110e-9;      % s
Vamp = 10;          % V
fg = 10e6;          % Hz

% Type of signal applied by the function generator
% Options: 'p'ulse, 's'ine wave, 't'riangle wave or s'q'uare wave
signal = 'p';

% Finite difference method used
% Options: 'a'dvance, 'r'etrograde or 'c'entered
finiteDiff = 'c';

% Simulation method: 'e'xplicit or 'i'mplicit
method = 'i';

% Simulation intervals
dt = 1e-10;         % s
t0 = 0;             % s
tf = 1000e-9;       % s
dx = 1e-1;          % m
x0 = 0;             % m
xf = 100;           % m

% Animation settings
filename = 'wave.gif';  % File to save the animation
Nmax = 200;             % Maximum amount of points
axisInc = Vamp/10;      % Quantization of axes

%---------------------------------------------
% Initial calculations
%---------------------------------------------

% Vetctors
x = x0:dx:xf;
t = t0:dt:tf;
J = length(x);
N = length(t);
v = zeros(N,J);     % Null initial conditions

% Function generator
switch (signal)
    case 's'
        Vg = @(t) Vamp*heaviside(t-Ton)*sin(2*pi()*fg*t);
    case 't'
        Vg = @(t) Vamp*heaviside(t-Ton)*sawtooth(2*pi()*fg*t,1/2);
    case 'q'
        Vg = @(t) Vamp*heaviside(t-Ton)*square(2*pi()*fg*t);
    otherwise %case 'p'
        Vg = @(t) Vamp*(heaviside(t-Ton)-heaviside(t-Toff));
end

% Parameters used in the simulation
nu = (dt/dx)^2*(1/(L*C));
tau = (G/C + R/L)*dt;
switch (finiteDiff)
    case 'a'
        alfa = nu/(1+tau);
        beta = -1/(1+tau);
        a = -(2 + 1/nu + tau/nu);
        b = (-2/nu - tau/nu + R*G*dx^2);
        c = 1/nu;
    case 'r'
        alfa = nu;
        beta = -1+tau;
        a = -(2 + 1/nu);
        b = (-2/nu + tau/nu + R*G*dx^2);
        c = (1/nu - tau/nu);
    otherwise % case 'c'
        alfa = nu/(1+tau/2);
        beta = -1+tau/(1+tau/2);
        a = -(2 + 1/nu + tau/(2*nu));
        b = (-2/nu + R*G*dx^2);
        c = (1/nu - tau/(2*nu));
end
gama = 1-beta-(2+R*G*dx^2)*alfa;
% Expressions for the frontier: forward method, because backwards results
% in implicit expressions and centered makes the simulation unstable
delta = Rg*dt/(L*dx);
ro = 1 - R*dt/L - Rg*dt/(L*dx);
theta = -1 + R*dt/L;
phi = Rl*dt/(L*dx);
epslon = 1 - R*dt/L - Rl*dt/(L*dx);
if (Rg == 0)
    d = 0;
    e = 0;
    f = 0;
    g = 1;
else
    d = -1 + L*dx/(Rg*dt);
    e = -R*dx/Rg + L*dx/(Rg*dt);
    f = R*dx/Rg - L*dx/(Rg*dt);
    g = L*dx/(Rg*dt);
end
if (Rl == 0)
    h = 0;
    i = 0;
else
    h = -1 + L*dx/(Rl*dt);
    i = -R*dx/Rl + L*dx/(Rl*dt);
end

%---------------------------------------------
% Simulation
%---------------------------------------------

if (method == 'e') % Explicit method
    for n = 2:(N-1)
        v(n+1,1) = delta*v(n,2) + ro*v(n,1) + Vg((n+1)*dt) + theta*Vg(n*dt);
        for j = 2:(J-1)
            v(n+1,j) = alfa*(v(n,j+1) + v(n,j-1)) + beta*v(n-1,j) + gama*v(n,j);
        end
        v(n+1,J) = phi*v(n,J-1) + epslon*v(n,J);
    end
else % Implicit method
    % Initialize matrices and vectors
    A = zeros(J,J);
    B = zeros(J,J);
    F = zeros(J,1);
    G = zeros(J,1);
    for j = 2:(J-1)
        A(j,j-1) = 1;
        A(j,j) = a;
        A(j,j+1) = 1;
        B(j,j) = b;
    end
    A(1,1) = (a-d);
    A(1,2) = 1;
    B(1,1) = (b-e);
    A(J,J) = (a-h);
    A(J,J-1) = 1;
    B(J,J) = (b-i);
    F(1) = -f;
    G(1) = -g;
    mA = A\B;
    mB = c*inv(A);
    mC = A\F;
    mD = A\G;
    
    % Apply numerical method
    for n = 2:(N-1)
        v(n+1,:) = mA*v(n,:)' + mB*v(n-1,:)' + mC*Vg(n*dt) + mD*Vg((n+1)*dt);
    end
end

%---------------------------------------------
% Graphs
%---------------------------------------------

% Calculations to simplify graph plotting
supV = ceil(max(v(:))/axisInc)*axisInc;
infV = floor(min(v(:))/axisInc)*axisInc;
if (N >= Nmax)
    inc = floor(N/Nmax);
    N2 = floor(N/inc)*inc;
else
    N2 = N;
    inc = 1;
end

% Graph 3D
figure('WindowState','maximized'); % Only available from version R2018a
gca.FontSize = 14;
surf(x,t(1:inc:N2),v(1:inc:N2,:),'edgecolor','none');
xlabel('Distance (m)','FontSize',16);
ylabel('Time (s)','FontSize',16);
zlabel('Voltage (V)','FontSize',16);

% Waveforms at the generator and load
figure('WindowState','maximized'); % Only available from version R2018a
gca.FontSize = 14;
plot(t,v(:,1));
hold on;
grid;
plot(t,v(:,end));
xlabel('Time (s)','FontSize',16);
ylabel('Voltage (V)','FontSize',16);
legend({'Generator', 'Load'},'FontSize',14);

% Waveform animation along the cable
figure('WindowState','maximized'); % Only available from version R2018a
gca.FontSize = 14;
handle_line = plot(x,v(1,:));
axis([x0,xf,infV,supV]);
grid;
xlabel('Distance (m)','FontSize',16);
ylabel('Voltage (V)','FontSize',16);
title('Propagation along the cable','FontSize',16);
frame = getframe(gcf);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
for n = (1+inc):inc:N2
    handle_line.YData = v(n,:);
    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
end

%---------------------------------------------
% End
%---------------------------------------------