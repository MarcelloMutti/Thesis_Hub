function [sc_param, dr_sc_param, P, dPdr] = MARGO_param(r)
% v1 (master)

% adimensional input [AU]
% adimensional output
% sc_param [kg*km/s^2, km/s]

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    MU=22.6;                                    % m0 [kg]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    g0=9.80665e-3;      % [km/s^2]

    a=[-0.7253 0.02481 0 0 0];
    b=[2652 -18.123 0.3887 -0.00174 0];
    c=[840.11 -1754.3 1625.01 -739.87 134.45];

    da=(1:4).*a(2:end);
    db=(1:4).*b(2:end);
    dc=(1:4).*c(2:end);

    Tmax=@(P) sum(P.^(0:4).*a,2); % [mN]
    Isp=@(P) sum(P.^(0:4).*b,2);  % [s]
    Sp=@(r) sum(r.^(0:4).*c,2);   % [W]
    % Pin=@(r) min(120,Sp(r));      % [W] to add prob.Plim(2)
    Pin=@(r) 120;

    dTmaxdPin=@(P) sum(P.^(0:3).*da,2); % [mN/W]
    dIspdPin=@(P) sum(P.^(0:3).*db,2);  % [s/W]
    dPindr=@(r) sum(r.^(0:3).*dc,2);    % [W/AU]

    sc_param=zeros(2,length(r));
    dr_sc_param=zeros(2,length(r));

    sc_param(1,:)=TU^2/(MU*LU)*Tmax(Pin(r))*1e-6;  % T [-]
    sc_param(2,:)=TU/LU*Isp(Pin(r))*g0;            % c [-]

    % P=Sp(r); % [W]
    P=120*ones(size(r));
    % dPdr=dPindr(r); % [W/AU]
    dPdr=0;

    % if Sp(r)<120
    % 
    %     dr_sc_param(1,:)=TU^2/(MU*LU)*dTmaxdPin(Pin(r)).*dPindr(r)*1e-6; % dT/dr
    %     dr_sc_param(2,:)=TU/LU*dIspdPin(Pin(r)).*dPindr(r)*g0;           % dIsp/dr
    % 
    % elseif Sp(r)>=120

        dr_sc_param(1,:)=0;
        dr_sc_param(2,:)=0;

    % end

    
    

end