function [dydPhi] = TwBP_EL(~, z)

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1
    MU=22.3;                                    % m0 [kg]
    cf=1e-6;    % conversion factor from [mN] to [kg*km/s^2]

    % adimensionalization param
    U=[LU TU MU cf];

    g0=9.80665e-3;                              % [km/s^2]
    
    % MARGO prop parameters
    MP=[[-0.7253 0.02481 0       0        0];       % a
        [2652    -18.123 0.3887  -0.00174 0];       % b
        [840.11  -1754.3 1625.01 -739.87  134.45]]; % c

    y=z(1:14);
    Phi=reshape(z(15:210),[14,14]);

    dy=F(y,MP,U,g0);
    A=DyF(y,MP,U,g0);

    dydPhi=zeros(210,1);
    dydPhi(1:14)=dy;
    dydPhi(15:210)=reshape(A*Phi,[14*14,1]);
end