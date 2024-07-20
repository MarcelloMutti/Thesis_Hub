function sc_param=MARGO_param(r)
% v1 (master)

% adimensional input [AU]
% dimensional output
% sc_param [kg*km/s^2, km/s]

% spice g0 computation
%     mu_E=cspice_bodvrd('Earth','GM',1);
%     R_E=cspice_bodvrd('Earth','RADII',3);
%     R=mean(R_E);
%     g0=mu_E/R^2;

    g0=9.80665e-3;      % [km/s^2]

    a=[-0.7253 0.02481 0 0 0];
    b=[2652 -18.123 0.3887 -0.00174 0];
    c=[840.11 -1754.3 1625.01 -739.87 134.45];
    
    Tmax=@(P) dot(P.^(0:4),a); % [mN]
    Isp=@(P) dot(P.^(0:4),b);  % [s]
    Pin=@(r) dot(r.^(0:4),c);  % [W]

    sc_param=zeros(2,1);
    sc_param(1)=Tmax(Pin(r))*1e-6;  % T [kg*km/s^2]
    sc_param(2)=Isp(Pin(r))*g0;     % c [km/s]

end