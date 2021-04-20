%%Marche aléatoire sur réseau 1D
%incrément dt cst et pas equiprobables
clear all, close all, clc
Npas=1000;
for i=1:10000
ax=randi(2,1,Npas); % pour la position des particules, à chaque pas
soit la particule avance ou recule dans un tirage entre 1 et 2
bx=2*ax-3; %coef directeur pour tirer entre -1 et 1
posx=[0 cumsum(bx)]; %position des particules
ay=randi(2,1,Npas);
by=2*ay-3;
posy=[0 cumsum(by)];
t=0:Npas; %vecteur temps
pos{i}=[posx; posy]; %structure qui regroupe l'ensemble des
trajectoires des particules selon x et y
end
for i=1:50 %la position x de 50 particules
plot(t,pos{i}(1,:),'-')
hold on
end
Q=zeros(1,1000) %matrice à 3 dimensions Q= matrice 2 lignes avec un
deplacement de 1:100 et 1000 = nombre de particules
for i=1:1000
Q(1,i)=pos{i}(1,10); %le temps =10, la matrice est remplie par les
positions des particules en x et y pour un temps donné qu'on met dans la
première page de la matrice
end
plot(t,Q(1,:,1:20),'o') %on trace le temps pour les particules de 1 à 20
hist(Q(1,400,:)) %histogramme à t=400
N = 1000;
displacement = randn(1,N);
plot(displacement)
hist(displacement, 25);
x = cumsum(displacement);
plot(x);
ylabel('position');
xlabel('time step');
title('Position of 1D Particle en fonction du Temps');
particle = struct();
particle.x = cumsum( randn(N, 1) );
particle.y = cumsum( randn(N, 1) );
plot(particle.x, particle.y);
ylabel('Y Position');
xlabel('X Position');
title('position in 2D');
dsquared = particle.x .^ 2 + particle.y .^ 2;
plot(dsquared);
d = 1.0e-6; % diameter in meters
eta = 1.0e-3; % viscosity of water in SI units (Pascal-seconds)
kB = 1.38e-23; % Boltzmann constant
T = 293; % Temperature in degrees Kelvin
D = kB * T / (3 * pi * eta * d)
tau = .1; % time interval in seconds
time = tau * 1:N;
% create a time vector for plotting
k = sqrt(2*D*tau);
dx = k * randn(N,1);
dy = k * randn(N,1);
x = cumsum(dx);
y = cumsum(dy);
dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
squaredDisplacement = ( x .^ 2) + ( y .^ 2);
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3); %plot theoretical
line
plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time for 1 Particle in 2 Dimensions');
plot(x,y); title('Particle Track of a Single Simulated Particle');
clf; hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3);
% plot theoretical line
plot(time, squaredDisplacement);
hold off;
xlabel('Time');
ylabel('Displacement Squared');
title('Displacement Squared versus Time with Bulk Flow');
% compute the ensemble average
rsquaredSum = zeros(1,N);
for i = 1:particleCount
rsquaredSum = rsquaredSum + particle{i}.rsquared;
end
ensembleAverage = rsquaredSum / particleCount;
% create the plot
clf; hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'b', 'LineWidth', 3); % plot theoretical
line
plot(time, ensembleAverage , 'k', 'LineWidth', 3); % plot ensemble average
legend('Theoretical','Average','location','NorthWest');
for i = 1:particleCount
plot(time, particle{i}.rsquared, 'color', rand(1,3)); % plot each particle
track
end
xlabel('Time (seconds)');
ylabel('Displacement Squared (m^2)');
title('Displacement Squared vs Time'); hold off;
N=1000;
d = 1.0e-6; % diameter in meters
eta = 1.0e-3; % viscosity of water in SI units (Pascal-seconds)
kB = 1.38e-23; % Boltzmann constant
T = 293; % Temperature in degrees Kelvin
D = kB * T / (3 * pi * eta * d)
tau = .1; % time interval in seconds
time = tau * 1:N; % create a time vector for plotting
k = sqrt(2*D*tau);
dx = k * randn(N,1);
dy = k * randn(N,1);
x = cumsum(dx);
y = cumsum(dy);
dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
squaredDisplacement = ( x .^ 2) + ( y .^ 2);
plot(x,y); title('Particle Track of a Single Simulated Particle');
clf;
hold on;
plot(time, (0:1:(N-1)) * 2*k^2 , 'k', 'LineWidth', 3); % plot theoretical
line
plot(time, squaredDisplacement);
hold off;
xlabel('Tau');
ylabel('Displacement Squared');
title('MSD en fonction de tau');