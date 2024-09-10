function ad_x = ADIM(x,m0)
% v1 (master)

% dimensional inputs
% adimensionalized outputs
% x [km, km/s]    
    
    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    MU=m0;                                      % m0 [kg]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1 

    ad_x=zeros(size(x));
    ad_x(1:3,:)=x(1:3,:)./LU;       % rr
    ad_x(4:6,:)=x(4:6,:)./(LU/TU);  % vv
    ad_x(7,:)=x(7,:)./MU;           % m

end