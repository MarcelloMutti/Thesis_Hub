function [tt,yy,S,H] = DispRes(lltf,t0,sc_param,targ,epsilon)

% testing branch
% dimensional input (but costates and final time)
% adimensional output

    odeopt=odeset('RelTol',1e-13,'AbsTol',1e-13);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    ll0=lltf(1:7);
    tf=lltf(8)*TU+t0;

    x0=[SEL2_ND(t0); sc_param(3)];

    [x0_ad, sc_param_ad]=ADIM(x0,sc_param);

    y0=[x0_ad; ll0];

    Phi0=eye(length(y0));
    vPhi0=reshape(Phi0,[length(y0)^2,1]);

    [tt, yy]=ode78(@(t,y) TwBP_EL(t,y,u,sc_param_ad),[0 (tf-t0)/TU],[y0; vPhi0],odeopt);

    ToF=(tf-t0)/86400;
    mf=yy(end,7)*sc_param(3);
    mp=(yy(1,7)-yy(end,7))*sc_param(3);

    ttd=tt*TU/86400;

    fprintf('Departure date: %s (%.1f MJD2000)\n',cspice_et2utc(t0,'C',3),et2MJD2000(t0))
    fprintf('Arrival date: %s (%.1f MJD2000)\n',cspice_et2utc(tf,'C',3),et2MJD2000(tf))
    fprintf('ToF: %.3f days, Final Mass: %.3f kg, Propellant Mass: %.3f kg\n\n',ToF,mf,mp)

    plot3D(t0,tt,yy,targ);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
    xtf_ad=ADIM([xtf; 1],sc_param);
    rvtf=xtf_ad(1:6);

    df=yy(end,1:6).'-rvtf;

    fprintf('Position error %.3e km\nVelocity error %.3e km/s\n',norm(df(1:3))*LU, norm(df(4:6))*LU/TU)

    S=SwFun(tt,yy,sc_param_ad,epsilon);
    H=Hamil(tt,yy,sc_param_ad);

    figure
    subplot(4,1,1)
    plot(ttd,S)
    xlim([ttd(1) ttd(end)])
    legend('S','Location','best')
    grid on
    grid minor

    subplot(4,1,2)
    plot(ttd,1.*(S<0)+0)
    xlim([ttd(1) ttd(end)])
    legend('u','Location','best')
    grid on
    grid minor
    
    subplot(4,1,3)
    plot(ttd,yy(:,7))
    xlim([ttd(1) ttd(end)])
    legend('m','Location','best')
    grid on
    grid minor

    subplot(4,1,4)
    plot(ttd,H)
    xlim([ttd(1) ttd(end)])
    legend('H','Location','best')
    grid on
    grid minor
    xlabel('t [days]')
end