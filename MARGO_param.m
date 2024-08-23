function [sc_param, dr_sc_param, P] = MARGO_param(r)
% v1 (master)

% adimensional input [AU]
% adimensional output
% sc_param [kg*km/s^2, km/s]

% spice g0 computation
%     mu_E=cspice_bodvrd('Earth','GM',1);
%     R_E=cspice_bodvrd('Earth','RADII',3);
%     R=mean(R_E);
%     g0=mu_E/R^2;

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    MU=22.3;                                    % m0 [kg]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    g0=9.80665e-3;      % [km/s^2]

    a=[-0.7253 0.02481 0 0 0];
    b=[2652 -18.123 0.3887 -0.00174 0];
    c=[840.11 -1754.3 1625.01 -739.87 134.45];

    da=(1:4).*a(2:end);
    db=(1:4).*b(2:end);
    dc=(1:4).*c(2:end);
    
%     Tmax=@(P) dot(P.^(0:4),a);
%     Isp=@(P) dot(P.^(0:4),b); 
%     Pin=@(r) dot(r.^(0:4),c); 

    Tmax=@(P) sum(P.^(0:4).*a,2); % [mN]
    Isp=@(P) sum(P.^(0:4).*b,2);  % [s]
    Pin=@(r) sum(r.^(0:4).*c,2);  % [W]

%     dTmaxdPin=@(P) dot(P.^(0:3),da);
%     dIspdPin=@(P) dot(P.^(0:3),db); 
%     dPindr=@(r) dot(r.^(0:3),dc);   

    dTmaxdPin=@(P) sum(P.^(0:3).*da,2); % [mN/W]
    dIspdPin=@(P) sum(P.^(0:3).*db,2);  % [s/W]
    dPindr=@(r) sum(r.^(0:3).*dc,2);    % [W/AU]

    sc_param=zeros(2,length(r));
    dr_sc_param=zeros(2,length(r));
    P=zeros(1,length(r));

    sc_param(1,:)=TU^2/(MU*LU)*Tmax(Pin(r))*1e-6;  % T [-]
    sc_param(2,:)=TU/LU*Isp(Pin(r))*g0;            % c [-]

    dr_sc_param(1,:)=TU^2/(MU*LU)*dTmaxdPin(Pin(r)).*dPindr(r)*1e-6; % dT/dr
    dr_sc_param(2,:)=TU/LU*dIspdPin(Pin(r)).*dPindr(r)*g0;           % dIsp/dr

    P=Pin(r); % [W]

end