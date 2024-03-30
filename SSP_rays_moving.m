function v = SSP_rays_moving(z,iprofile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MF208 Aeroacoustic and acoustic propagation in moving media - 2021
% Practical work 3 - Ray-tracing code in a moving layered atmosphere
% Calculation of the sound speed and wind speed profiles
% output parameters :
% v(1,:): sound speed profile c(z)
% v(2,:): vertical sound speed gradient dc/dz
% v(3,:): velocity profile Vx(z)
% v(4,:): vertical velocity gradient dVx/dz
% input parameters :
% z: vector of heights
% iprofile: index defining the sound speed and velocity profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iprofile == 1 % profile 1: logarithmic increase for Vx - constant c
    z1= 1300;
    e=0.0074;
    cm=1500;
    v = zeros(4,length(z));
    v(1,:)= cm*(1+e*((z-z1)/(z1/2)+exp(-(z-z1)/(z1/2))-1));
    v(2,:)= cm*(0+e*(2/z1-2/z1*exp(-(z-z1)/(z1/2))));
    v(3,:)= 0;
    v(4,:)= 0;
elseif iprofile == 2 % profile 2: logarithmic increase for c - Vx=0
    c0 = 343;
    hv=0.2;
    zref=80;
    Vref=8;
    m=0.55;
    bv=Vref*((hv/zref)^m)*(1/log(2));
    v = zeros(4,length(z));
    v(1,:)= c0;
    v(2,:)= 0.;
    v(3,:)= (Vref.*(z/zref).^m).*(z>=hv) + (bv.*log(1+z/hv)).*(z<hv) ;
    v(4,:)=(m*Vref/zref.*(z/zref).^(m-1)).*(z>=hv) +  (bv/hv.*1./(1+z/hv)).*(z<hv); 
elseif iprofile == 3 % profile 2: logarithmic increase for c - Vx=0
    vref = 8;
    c0 = 343;
    zref = 80;
    hv = 0.2;
    m = 0.55;
    bv = (1/log(2)) * vref * ((hv/zref)^m);
    [p,q] = size(z);
    c = zeros(p,q);
    grad = zeros(p,q);
    for i = 1:q
        if z(i)>=hv
            c(i) = vref * ((z(i)/zref)^m);
            grad(i) = m*vref*(z(i)^(m-1)) / (zref^(m));
        else
            c(i) = bv * log(1 + (z(i)/hv));
            grad(i) = bv / (hv + z(i));
        end
    end
    v = zeros(4,length(z));
    v(1,:)= c0;
    v(2,:)= 0.;
    v(3,:)= c;
    v(4,:)= grad; 
else
    disp('profile not recognized')
    return
end

