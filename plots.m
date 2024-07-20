function plots(llt,t0,sc_param)

% testing branch
% dimensional input (but costates)
% adimensional output

    odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    l0_g=llt(1:7);
    tf=llt(8);

    x0=[SEL2_ND(t0); sc_param(3)];

    [x0_ad, sc_param_ad]=ADIM(x0,sc_param);

    y0=[x0_ad; l0_g];

    [tt, yy]=ode78(@(t,y) TwBP_EL(t,y,u,sc_param_ad),[0 (tf-t0)/TU],y0,odeopt);

    plot3D(yy);

    S=SwFun(tt,yy,sc_param_ad);
    H=Hamil(tt,yy,sc_param_ad);

    figure
    subplot(3,1,1)
    plot(S)
    legend('S')
    grid on
    grid minor
    subplot(3,1,3)
    plot(H)
    legend('H')
    grid on
    grid minor
    subplot(3,1,2)
    plot(yy(:,7))
    legend('m')
    grid on
    grid minor
end