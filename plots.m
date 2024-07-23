function plots(llt,t0,sc_param,targ)

% testing branch
% dimensional input (but costates)
% adimensional output

    odeopt=odeset('RelTol',1e-13,'AbsTol',1e-13);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    l0_g=llt(1:7);
    tf=llt(8);

    x0=[SEL2_ND(t0); sc_param(3)];

    [x0_ad, sc_param_ad]=ADIM(x0,sc_param);

    y0=[x0_ad; l0_g];

    [tt, yy]=ode78(@(t,y) TwBP_EL(t,y,u,sc_param_ad),[0 (tf-t0)/TU],y0,odeopt);

    plot3D(t0,tt,yy,targ);

    S=SwFun(tt,yy,sc_param_ad);
    H=Hamil(tt,yy,sc_param_ad);

    figure
    subplot(4,1,1)
    plot(tt,S)
    legend('S')
    grid on
    grid minor

    subplot(4,1,2)
    plot(tt,1.*(S<0)+0)
    legend('u')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(tt,yy(:,7))
    legend('m')
    grid on
    grid minor

    subplot(4,1,4)
    plot(tt,H)
    legend('H')
    grid on
    grid minor
end