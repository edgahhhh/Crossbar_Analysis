%{
    modeling car as source + sink flow to see what the air is doing at
    location of crossbars. Taking that information and analyzing the Re at
    that location for the different levels of freestream velocity. Using
    this to find the CD and ultimately power required to move crossbars.

    You can run the whole script altogether but I would recomend running
    each section on their own. 
%}

clear
%% modeling car as source+sink potential flow
% These first lines of code were from the given code to plot the rankine
% oval using source and sinks. The numbers were just tweaked slightly
mind = -2;
maxd = 6;
n = 1000;

x = linspace(mind, maxd, n);
y = linspace(mind, maxd, n);
[X, Y] = meshgrid(x, y);

% I didn't tweak any of the mask stuff because I've never used it before
R = sqrt(X.^2 + Y.^2);
dec = 1:10:n;
mask = false(n, n);
mask(1:10:n, 1:10:n) = true;
mask(R < 0.1) = false;

% define freestream velocity
Vinfinity =5;     % m/s

% choose points where the sources and sinks will be places
x1s1 = 0.25;            % source 1 @ x = 0.25 m
x2s1 = x1s1 + 1.4;      % sink 1 will be 1.4 m after source

% locations of second source and sink
x1s2 = 3.75;            % source 2 @ x = 3.75 m
x2s2 = x1s2 + 0.125;    % sink 2 located 0.125 m after source

% lambdas I found to be most accurate model through trial and error
lambda1 = Vinfinity * 0.5;    % pair 1
lambda2 = Vinfinity * 4.0;    % pair 2

% defining different 'x' distances from sources and sinks
x1 = X - x1s1; x2 = X - x1s2; x3 = X - x2s1; x4 = X - x2s2;
% define the different potentials
phi3 = Vinfinity * X;
phi1 = (lambda1/(2*pi)).*log(((x1.^2+Y.^2)./(x2.^2+Y.^2)).^0.5);
phi2 = (lambda2/(2*pi)).*log(((x3.^2+Y.^2)./(x4.^2+Y.^2)).^0.5);
phi = phi3 + phi1 + phi2;

% u from derivative of psy, d/dy(psy)
u = (Vinfinity * ones(size(X))) + ...
    (lambda1/(2*pi)).* ...
    ((x1./(x1.^2+Y.^2))-(x2./(x2.^2+Y.^2))) +...
    (lambda2/(2*pi)).* ...
    ((x3./(x3.^2+Y.^2))-(x4./(x4.^2+Y.^2)));

% v from derivate of psy, d/dx(psy)
v = -1 .*( (lambda1/(2*pi)) .* ...
    ((Y./(x2.^2+Y.^2))-(Y./(x1.^2+Y.^2))) + ...
    (lambda2/(2*pi)) .* ...
    ((Y./(x4.^2+Y.^2))-(Y./(x3.^2+Y.^2))) );

% plotting streamlines, velocity, and potential. This is directly from the
% rankine oval code as well, I just tweaked it to plot more streamlines
start_x = ones(50, 1) * -1;
start_y = linspace(-1, 1, 50);
normalized_phi = phi - min(phi(:));
normalized_phi = normalized_phi / max(normalized_phi(:));
figure
hold on
contour(X, Y, phi, 50, '--', LineWidth=2)
quiver(X(mask), Y(mask), u(mask), v(mask), LineWidth=2)
streamline(X, Y, u, v, start_x, start_y)
legend('\phi', 'velocity', '\psi', fontsize=16)
axis equal
xlabel('x position (m)')
ylabel('y position (m)')
xlim([-1, 4])
ylim([-1.5,1.5])
title('freestream of 5 m/s')

% my attempt at obtaining the velocity field at point of interests
%{
myX = 2.15 ; myY = 1.14;
for i = 1:n
    tempX = X(i,i);
    testX = myX - tempX;
    if abs(testX) < 0.0001
        indexX = i;
        break
    end
end
for b = 1:n
    tempY = Y(b,b);
    testY = myY - tempY;
    if abs(testY) < 0.0001
        indexY = b;
        break
    end
end
%}

%% Finding Re and plotting the found properties
% freestream velocities
Vfstm = linspace(5, 30 , 6);
% 'u' found at point of first  crossbar
U = [7.37776 14.7555 22.1333 29.511 36.8888 44.2665];
% 'v' found at point of first 
V = [1.50424 3.00847 4.51271 6.01694 7.52118 9.02541];


% finding Re for each 'u'
% define constants
rho = 1.293; d = 0.03 ; mu = 1.729e-5;
Re = (rho * d)/mu ; Re = U .* Re;

% subplots
%{

figure(2);
% plotting u vs Vinf in top row
subplot(2,1,1)
plot(U ,Vfstm , 'mo', LineWidth = 2)
ylim([0 , 35])
xlabel('u (m/s)'); ylabel('Vinf (m/s)')
title('horizontal velocity components at crossbars')

% plotting v vs Vinf in bottom row
hold on;
subplot(2,1,2)
plot(V , Vfstm , 'mo' , LineWidth = 2)
ylim([0 , 35])
xlabel('v (m/s)'); ylabel('Vinf (m/s)')
title('vertical velocity components at crossbars')
%}
% plotting u and v vs Vinf on same plot with legend
figure(3);
subplot(2 , 1 ,1 )
plot(U ,Vfstm , 'mo', LineWidth = 2)
hold on
plot(V , Vfstm , 'ro' , LineWidth = 2)
hold off
ylim([0 , 35]); title('u & v vs. freestream')
legend('u','v'); ylabel('Vinf (m/s)'); xlabel('(m/s)')

% plotting Re under previous plot
hold on
subplot(2, 1 , 2)
plot(Re , U , 'ro' , LineWidth = 2)
ylim([0,50]); xlim([0,10.5^(5)])
ylabel('u (m/s)'); xlabel('Re'); title('Re vs. u')
hold off

% CD for each Re from graphs, lower and upper
upperCD = [1.0 1.0 1.0 1.0 1.0 1.0];
lowerCD = [0.7 0.7 0.7 0.7 0.7 0.7];
% calculating drag
qbar = 1/2 .* rho .* U.^2; s = 1.7 * 0.03;
upperD = upperCD .* qbar .* s; lowerD = lowerCD .* qbar .* s;
upperP = upperD .* U; lowerP = lowerD .* U;
% finding mean values for the cD, drag, and power
meanCD = (upperCD + lowerCD) ./ 2; meanD = (upperD + lowerD) ./ 2;
meanP = (upperP + lowerP) ./ 2 ;

% Plotting CD, D, and P vs freestream
figure(4);
% subplot top left corner for CD
subplot(2, 2, 1); plot(Vfstm, lowerCD, 'ro', LineWidth = 2)
hold on; plot(Vfstm, upperCD, 'mo', LineWidth = 2); hold off
hold on ; plot(Vfstm, meanCD, LineWidth = 2); hold off
xlim([0,35]); ylim([0,1.2]); title('freestream vs. cD')
xlabel('Vinf (m/s)'); ylabel('cD');legend('lower cD','upper cD','avg cD')

% subplot top right corner for drag
subplot(2, 2, 2); plot(Vfstm, lowerD, 'ro', LineWidth = 2)
hold on; plot(Vfstm, upperD, 'mo', LineWidth = 2); hold off
hold on; plot(Vfstm, meanD, LineWidth = 2); hold off
xlim([0,35]); ylim([0,70]); title('freestream vs. drag')
xlabel('Vinf (m/s)');ylabel('D (N)');legend('lower D','upper D','avg D')

% subplot bottom left corner for power
subplot(2,2,3); 
plot (Vfstm, lowerP, 'ro',Vfstm, upperP, 'mo',Vfstm,meanP,lineWidth = 2)
xlim([0 , 35]); ylim([0 , 3*10^3]); title('freestream vs. power')
xlabel('Vinf (m/s)');ylabel('P (W)');legend('lower P','upper P','avg P')




