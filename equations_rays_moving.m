function Up=equations_rays_moving(U,k0x,iprofile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Benjamin COTTE - 29/03/2022
% MF208 Aeroacoustic and acoustic propagation in moving media - 2022
% Practical work 3 - Ray-tracing code in a stratified moving medium
% System of equations to solve for ray-tracing for vector U=[x,z,kz]
% output parameters :
% Up: vector F(U)
% input parameters :
% U: vector U=[x,z,kz]
% k0x: wavenumber k projected over x
% iprofile: index defining the sound speed and velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zloc=U(2); % height z

% calculate c(z), dcdz(z), V(z) and dVdz(z)
v=SSP_rays_moving(zloc,iprofile); 
c=v(1);
dcdz=v(2);
Vx=v(3);
dVxdz=v(4);

% ray-tracing equations 
k = sqrt(k0x^2 + U(3)^2); % U(3) = kz
Up(1) =  k0x*c/k + Vx;
Up(2) = U(3)*c/k; % U(3) = kz
Up(3) = -k*dcdz - k0x*dVxdz;
Up=Up';
end
