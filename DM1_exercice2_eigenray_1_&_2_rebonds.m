%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin COTTE - 29/03/2022
% MF208 Aeroacoustic and acoustic propagation in moving media - 2022
% Practical work 3 - Ray-tracing code in a stratified moving atmosphere
% Calculation for a point source in the x-z plane
% Numerical solution based on a 4th order Runge-Kutta scheme
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% program used to run the
code = 1; % run with Matlab

% Input parameters for eigenray search
zS = 40;  % source height (m)
zR = 1.5;  % receiver height (m)
xS = 0;     % source range (m)
xR = 2000; % receiver range (m)
nb_refl_eigenray = 1; % number of reflections on the surface z=0

% profile 1: logarithmic increase for V - constant c
% profile 2: logarithmic increase for c - V = 0
iprofile = 3;

tetamin = -30*pi/180; % initial ray angle (radians)
tetamax = 30*pi/180; % initial ray angle (radians)
teta0 = (tetamin+tetamax)/2;

% plot sound speed and wind speed profiles
zprofile = 0:0.1:150;
A = SSP_rays_moving(zprofile,iprofile);
Vprofil = A(3,:); % wind speed profile
dVprofil_dz = A(4,:); % dV/dz
clear A

h = figure(1);
set(h,'Position',[200 200 400 400])
subplot(121)
plot(zprofile,Vprofil,'k','LineWidth',2)
hold on
%plot([1500 1560],[1.3 1.3],'k--','LineWidth',2)
set(gca,'FontSize',15)
ylabel('V(z) (m/s)')
xlabel('z (m)')
xlim([0 150])
grid on

% h = figure(2);
% set(h,'Position',[200 200 400 400])
subplot(122)
plot(zprofile,dVprofil_dz,'k','LineWidth',2)
hold on
%plot([1500 1560],[1.3 1.3],'k--','LineWidth',2)
set(gca,'FontSize',15)
ylabel('dV/dz (/s)')
xlabel('z (m)')
xlim([0 150])
grid on
title('Profils dela vitesse du milieu et de son gradient vertical')

% sound speed c and horizontal wind speed Vx at source height
A = SSP_rays_moving(zS,iprofile);
cs = A(1); 
Vxs = A(3); 
clear A

t0=0; % initial time
tmax=(xR-xS)/cs*1.1; % maximum travel time
dt = 5e-4; % time step
Nt = ceil(tmax/dt); % number of time iterations

% precision desired
Deltax = 1; % <1m
Deltaz = 0.1; % <0.1m
eigenrayfound = 0;
niter = 1;

while ((eigenrayfound == 0) & (niter < 30))
    % initialize variables
    x = zeros(1,Nt); % horizontal distance
    z = zeros(1,Nt); % vertical distance
    kz = zeros(1,Nt); % wavenumber k projected over z
    nb_refl = 0;
    xrefl = [];
    ray_length = zeros(1,Nt);
    travel_time = zeros(1,Nt);

    % initial parameters for a source close to the ground
%     teta0=-tetamax+(in-1)*dteta; % initial ray direction
    k0x = cos(teta0)/(cs+Vxs*cos(teta0)); % wavenumber k projected over x (omega arbitrarily set to 1)
    U=[0 zS sin(teta0)/(cs+Vxs*cos(teta0))]'; % vector U=[x,z,kz] at t=0
    x(1) = U(1);
    z(1) = U(2);
    kz(1)= U(3);

    for it=1:Nt-1 % loop over time
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
            nb_refl=nb_refl+1; % count number of reflections of ray in
            % position of reflection obtained by interpolation
            slope_inter = -z(it-1)/(z(it)-z(it-1));
            x_inter = x(it-1) + slope_inter*(x(it)-x(it-1)); 
            kz_inter = kz(it-1) + slope_inter*(kz(it)-kz(it-1)); 
            x(it) = x_inter;
            z(it) = 0.;
            kz(it) = -kz_inter; % direction of specular reflection
            xrefl = [xrefl x_inter];
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
        ray_length(it) = ray_length(it-1) + dL;
        travel_time(it) = travel_time(it-1) + dL/vg;
    end;
    figure(11);
    plot(xR,zR,'ro-','LineWidth',1)
    hold on
    plot(x,z,'k-','LineWidth',1)
    set(gca,'FontSize',13)
    xlabel('x (m)')
    ylabel('z (m)')
    grid on
    
    % eigenray found
    [temp,indx] = min(abs(x-xR));
    if (nb_refl > nb_refl_eigenray)%there are more reflections than wanted to reach receiver
        if (xrefl(end) < xR) % reflection before the receiver, as we want it after, we have to increase teta0
            tetamin = teta0; % initial ray angle (radians)
            teta0 = (tetamin+tetamax)/2;
            disp(['depth too low: teta0 increased to ',num2str(teta0*180/pi),'deg'])
        else  % reflection after the receiver
            if ((abs(x(indx)-xR) < Deltax) && (abs(z(indx)-zR) < Deltaz))%receiver in the ray's path before 2nd reflection
                plot(x,z,'b-','LineWidth',2)
                title('Recherche du rayon propre pour m=0.55 et zs=40m');
                xlim([0 2100])
                disp(['eigenray found: teta0 = ',num2str(teta0*180/pi),'deg'])
                disp(['ray length = ',num2str(ray_length(indx)),'m and travel time = ',num2str(travel_time(indx)*1000,'%.1f'),'ms']);
                return
            elseif z(indx) < zR % increase teta0, 
                tetamin = teta0; % initial ray angle (radians)
                teta0 = (tetamin+tetamax)/2;
                disp(['depth too low: teta0 increased to ',num2str(teta0*180/pi),'deg'])
            elseif z(indx) > zR % decrease teta0 
                tetamax = teta0; % initial ray angle (radians)
                teta0 = (tetamin+tetamax)/2;
                disp(['depth too high: teta0 decreased to ',num2str(teta0*180/pi),'deg'])
            end
        end
        
    elseif (nb_refl_eigenray>0 && nb_refl==nb_refl_eigenray) %right number of reflections
        if (xrefl(end) < xR) %last reflection before receiver
            if ((abs(x(indx)-xR) < Deltax) && (abs(z(indx)-zR) < Deltaz))%receiver in the ray's path before 2nd reflection
                plot(x,z,'b-','LineWidth',2)
                title('Recherche du rayon propre pour m=0.55 et zs=40m');
                xlim([0 2100])
                disp(['eigenray found: teta0 = ',num2str(teta0*180/pi),'deg'])
                disp(['ray length = ',num2str(ray_length(indx)),'m and travel time = ',num2str(travel_time(indx)*1000,'%.1f'),'ms']);
                return
            elseif z(indx) > zR % increase teta0, 
                tetamin = teta0; % initial ray angle (radians)
                teta0 = (tetamin+tetamax)/2;
                disp(['depth too low zx > zR: teta0 increased to ',num2str(teta0*180/pi),'deg'])
            elseif z(indx) < zR % decrease teta0 
                tetamax = teta0; % initial ray angle (radians)
                teta0 = (tetamin+tetamax)/2;
                disp(['depth too high zx < zR: teta0 decreased to ',num2str(teta0*180/pi),'deg'])
            end
        else % (xrefl(end) >= xR), last reflection after receiver, we want it before, decrease teta0
            tetamax = teta0; % initial ray angle (radians)
            teta0 = (tetamin+tetamax)/2;
            disp(['depth too high: teta0 decreased to ',num2str(teta0*180/pi),'deg'])                 
        end
        
    elseif(nb_refl_eigenray ==0)%necessarily nb_refl==0 because >=0 and not >nb_refl_eigenray
        if ((abs(x(indx)-xR) < Deltax) && (abs(z(indx)-zR) < Deltaz)) 
            plot(x,z,'b-','LineWidth',2)
            title('Recherche du rayon propre pour m=0.55 et zs=40m');
            xlim([0 2100])
            disp(['eigenray found: teta0 = ',num2str(teta0*180/pi),'deg'])
            disp(['ray length = ',num2str(ray_length(indx)),'m and travel time = ',num2str(travel_time(indx)*1000,'%.1f'),'ms']);
            return
        elseif z(indx) < zR % increase teta0
            tetamin = teta0; % initial ray angle (radians)
            teta0 = (tetamin+tetamax)/2;
            disp(['depth too low: teta0 increased to ',num2str(teta0*180/pi),'deg'])
        elseif z(indx) > zR % decrease teta0
            tetamax = teta0; % initial ray angle (radians) 
            teta0 = (tetamin+tetamax)/2;
            disp(['depth too high: teta0 decreased to ',num2str(teta0*180/pi),'deg'])
        end
    else %nb_refl_eigenray>0 and nb_refl < nb_refl_eigenray : decrease teta0
        tetamax = teta0; % initial ray angle (radians) 
            teta0 = (tetamin+tetamax)/2;
            disp(['depth too high: teta0 decreased to ',num2str(teta0*180/pi),'deg'])
    end
    
    niter = niter + 1;
end
disp(['niter = ',num2str(niter)])

figure(12);
hold on
plot(x(1:indx),z(1:indx),'b-','LineWidth',2) % eigenray
plot([0 xR],[zS zR],'k-','LineWidth',2) % eigenray in homogeneous medium
plot(xR,zR,'ro-','LineWidth',2)
set(gca,'FontSize',15)
xlabel('x (m)')
ylabel('z (m)')
grid on
xlim([0 xR])
legend('downward refraction','homogeneous conditions')

% values in a homogeneous atmosphere with c0=340m/s
c0 = 343;
teta0_hom = atan((zR-zS)/xR)*180/pi
ray_length_hom = sqrt(xR^2 + (zR-zS)^2)
travel_time_hom = ray_length_hom/c0