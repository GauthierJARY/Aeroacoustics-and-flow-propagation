%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin COTTE - 02/04/2023
% MF208 Aeroacoustic and acoustic propagation in moving media - 2023
% Practical work 3 - Ray-tracing code in a stratified moving atmosphere
% Calculation for a point source in the x-z plane
% Numerical solution based on a 4th order Runge-Kutta scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


% Input parameters
Nrays = 30; % number of rays (even number)
tetamax = 0.1; % aperture angle of the ray beam (degrees)
nb_souhaite_rebond=2; % number of reflection wanted before receptor
teta_init=2*pi/180; % angle autour duquel on cherche le rayon propre (désymétrise le cône de dispersion
zs = 40;  % source height (m)
L = 3000; % calculation distance (curvilinear distance along rays in m)
zR=0.2;
xR=2000;

% profile 1: exercice 1
% profile 2: exercice 2 v1
% profile 3: exercice 2 v2

iprofile = 3;

% save parameters
savedir = ''; % directory where results are saved
fname = ['rays_zS',num2str(round(zs)),'m_tetamax',num2str(round(tetamax)),'deg_profile',num2str(iprofile)]

tetamax=tetamax*pi/180;    % aperture angle of the ray beam (radians)
dteta=2*tetamax/(Nrays-1); % angular step between 2 rays

% sound speed c and horizontal wind speed Vx at source height
A = SSP_rays_moving(zs,iprofile);
cs = A(1);
Vxs = A(3); 
clear A
t1=1.9849;
t2=8.08;
t3=8.13;
t=[t1,t2,t3];
t=t*pi/180;
t0=0; % initial time
tmax=L/cs; % maximum travel time (s)
dt = 0.002 % time step (s)
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
for in=1:3 % loop over rays
    % initialize variables
    x = zeros(1,niter); % horizontal distance
    z = zeros(1,niter); % vertical distance
    kz = zeros(1,niter); % wavenumber k projected over z
    
    % initial parameters for a source close to the ground
    teta0=t(in); % initial ray direction
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
        if test_reflection < 0  % change of sign = reflection
            if x(it)<xR % on ne veut compter que les réflexions avant le récepteur pour pouvoir discriminer notre recherche
                % disp(['coord rebonds', num2str(z(it))])
                nb_refl(in)=nb_refl(in)+1; % count number of reflections of ray in
            end
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
    xR=2000;
    Deltax=1;
    Deltaz=0.1;
    [temp,indx] = min(abs(real(x-xR)));
%     fprintf('impact abs: %i, hauteur: %i\n',x(indx),z(indx));
    if ((abs(x(indx)-xR) < Deltax) && (abs(real(z(indx))-zR) < Deltaz) && nb_souhaite_rebond==nb_refl(in) )% rayon propre qui touche le récepteur avec les bonnes conditions souhaitées
        disp('Hourra!!!!!!!!!!!!!!!!!!!!!!')
        figure();
        plot(x,z,'b-','LineWidth',2)
        title('Rayon propre avec rebond pour m=0.55 avec zs=40m');
        xlim([0 2100])
        hold on
        plot(xR,zR,'ro-','LineWidth',2)
        set(gca,'FontSize',15)
        xlabel('x (m)')
        ylabel('z (m)')
        grid on
        disp(['eigenray found: teta0 = ',num2str(teta0*180/pi),'deg'])
        disp(['ray length = ',num2str(ray_length(in)),'m and travel time = ',num2str(travel_time(in)*1000,'%.1f'),'ms']);
    end
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
hold on
plot(xR,zR,'ro-','LineWidth',2)
set(gca,'FontSize',15)
xlabel('x (m)')
ylabel('z (m)')
grid on

set(gca,'FontSize',15)
xlabel('x (m)')
ylabel('z (m)')
grid on
xlim([0 2002])
%ylim([0 15000])
title('Tracé de rayons pour une source à hauteur éolienne');


