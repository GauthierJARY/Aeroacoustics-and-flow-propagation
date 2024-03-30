%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin COTTE - 02/04/2023
% MF208 Aeroacoustic and acoustic propagation in moving media - 2023
% Practical work 3 - Ray-tracing code in a stratified moving atmosphere
% Calculation for a point source in the x-z plane
% Numerical solution based on a 4th order Runge-Kutta scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% program used to run the
code = 1; % run with Matlab
% code = 2; % run with Octave

% Input parameters
Nrays = 11; % number of rays (even number)
tetamax = 30; % aperture angle of the ray beam (degrees)
zs = 1300;      % source height (m)
L = 150000; % calculation distance (curvilinear distance along rays in m)

% profile 1: logarithmic increase for V - constant c
% profile 2: logarithmic increase for c - V = 0
% profile 3: logarithmic decrease for V - constant c
% profile 4: logarithmic decrease for c - V = 0
% profile 5: linear increase for V - constant c
% profile 6: linear increase for c - V = 0
iprofile = 1;

% save parameters
savedir = ''; % directory where results are saved
fname = ['rays_zS',num2str(round(zs)),'m_tetamax',num2str(round(tetamax)),'deg_profile',num2str(iprofile)]

tetamax=tetamax*pi/180;    % aperture angle of the ray beam (radians)
dteta=2*tetamax/(Nrays-1); % angular step between 2 rays

% plot sound speed and wind speed profiles
zprofile = 0:5000;
A = SSP_rays_moving(zprofile,iprofile);
cprofil = A(1,:); % sound speed profile
gradprofil = A(2,:); % wind speed profile
clear A

h = figure(2);
set(h,'Position',[200 200 400 400])
plot(zprofile,cprofil,'k','LineWidth',2)
set(gca,'FontSize',15)
ylabel('c(z) (m/s)')
xlabel('z (m)')
title('Profil de la célérité du son en fonction de la profondeur')
grid on

h = figure(3);
set(h,'Position',[800 200 400 400])
plot(zprofile,gradprofil,'k','LineWidth',2)
set(gca,'FontSize',15)
ylabel('dc/dz (s-1)')
xlabel('z (m)')
title('Profil du gradient de la célérité en fonction de la profondeur')
grid on

% sound speed c and horizontal wind speed Vx at source height
A = SSP_rays_moving(zs,iprofile);
cs = A(1);
Vxs = A(3); 
clear A

t0=0; % initial time
tmax=L/cs; % maximum travel time (s)
dt = 1/16; % time step (s)
% converge a partir de 10e-5
niter = ceil(tmax/dt); % number of time iterations

% number of reflections for each ray
zmax_rays = zeros(1,Nrays);
nb_refl = zeros(1,Nrays);
ray_length = zeros(1,Nrays);
travel_time = zeros(1,Nrays);

% initialize storage variables
xrays = zeros(Nrays,niter);
zrays = zeros(Nrays,niter);


tic
for in=1:Nrays % loop over rays
    % initialize variables
    x = zeros(1,niter); % horizontal distance
    z = zeros(1,niter); % vertical distance
    kz = zeros(1,niter); % wavenumber k projected over z
    
    % initial parameters for a source close to the ground
    teta0=-tetamax+(in-1)*dteta; % initial ray direction
    k0x = cos(teta0)/(cs+Vxs*cos(teta0)); % wavenumber k projected over x (omega arbitrarily set to 1)
    U = [0 zs sin(teta0)/(cs+Vxs*cos(teta0))]'; % vector U=[x,z,kz] at t=0
    x(1) = U(1);
    z(1) = U(2);
    kz(1)= U(3);

    for it=1:niter-1 % loop over time
        % time integration with 4th order Runge Kutta scheme
        k1 = equations_rays_moving(U        ,k0x,iprofile);
        k2 = equations_rays_moving(U+dt/2*k1,k0x,iprofile);
        k3 = equations_rays_moving(U+dt/2*k2,k0x,iprofile);
        k4 = equations_rays_moving(U+dt*k3  ,k0x,iprofile);
        U = U + dt*(k1 + 2*k2 + 2*k3 + k4)/6.;
        
        % solution at iteration it+1
        it = it+1;
        x(it)  = U(1);
        z(it)  = U(2);
        kz(it) = U(3);
        % test if there is a reflection between time steps it and it+1
        test_reflection = z(it)*z(it-1);
        if test_reflection < 0 % change of sign = reflection
            nb_refl(in)=nb_refl(in)+1; % count number of reflections of ray in
            % position of reflection obtained by interpolation
            slope_inter = -z(it-1)/(z(it)-z(it-1));
            x_inter = x(it-1) + slope_inter*(x(it)-x(it-1)); 
            kz_inter = kz(it-1) + slope_inter*(kz(it)-kz(it-1)); 
            x(it) = x_inter;
            z(it) = 0.;
            kz(it) = -kz_inter; % direction of specular reflection
            % new vector U at iteration it+1
            U(1) = x(it);
            U(2) = z(it);
            U(3) = kz(it);
        end
        % group velocity at z(it)
        v=SSP_rays_moving(z(it),iprofile);
        c=v(1);
        Vx=v(3);
        vg = sqrt( (Vx+c*k0x/sqrt(k0x^2+kz(it)^2))^2 + c^2*kz(it)^2/(k0x^2+kz(it)^2) );
        
        % update ray length and travel time
        dL = sqrt((x(it)-x(it-1))^2 + (z(it)-z(it-1))^2); % ray length
        ray_length(in) = ray_length(in) + dL;
        travel_time(in) = travel_time(in) + dL/vg;
    end
    % convergence test based on the maximum height
    zmax_rays(in) = max(z);  

    % store ray coordinates    
    xrays(in,:) = x;
    zrays(in,:) = z;
end

cputime = toc


color_rays = jet(Nrays);
figure(11);
hold on
for k=1:Nrays % loop over rays
  plot(xrays(k,:),zrays(k,:),'Color',color_rays(k,:),'LineWidth',2)
end

% % traits noirs sur les courbes de rayon : front d'onde
% for it=1:200:niter % loop over rays
%   plot(xrays(:,it),zrays(:,it),'k-','LineWidth',2)
% end


set(gca,'FontSize',15)
xlabel('x (m)')
ylabel('z (m)')
grid on
xlim([0 100000])
ylim([0 15000])
axis ij
title('Tracé de 11 rayons pour une source à la profondeur 200m, angle 90°');
disp(['maximum height of the 1st ray: ',num2str(zmax_rays(1),'%.2f'),'m - travel time = ',num2str(travel_time(1)*1000,'%.2f'),'ms'])
disp(['maximum height of the last ray: ',num2str(zmax_rays(Nrays),'%.2f'),'m - travel time = ',num2str(travel_time(Nrays)*1000,'%.2f'),'ms'])

if code == 1 % Matlab
    save([savedir fname],'Nrays','xrays','zrays','niter','zs','tetamax','L','dt')
elseif code == 2 % Octave
    save('-mat-binary',[savedir fname '.mat'],'Nrays','xrays','zrays','niter','zs','tetamax','L','dt')
end


