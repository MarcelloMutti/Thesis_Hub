function [prob] = DispRes(prob,output)

% testing branch
% dimensional input (but costates and final time)
% adimensional output

    if nargin<2
        output=1;
    end

    % Set the default font size for axes labels, tick labels and legend
    set(0,'DefaultAxesFontSize',15);
    set(0,'DefaultLegendFontSize',15);
    
    % Set the default linewidth
    set(0,'DefaultLineLineWidth',2);
    
    % Set the default text and legend interpreter
    set(groot,'defaultTextInterpreter','latex');
    set(groot,'defaultLegendInterpreter','latex');

    g0=9.80665e-3;      % [km/s^2]

    m0=prob.m0;
    t0=prob.t0;
    targ=prob.targ;
    epsilon=prob.epsilon;

    Pmin=prob.Plim(1);
    Pmax=prob.Plim(2);

    odeopt=odeset('RelTol',1e-13,'AbsTol',1e-13);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1
    MU=prob.m0;

    u=1;

    tf=prob.tf;

%     x0=[SEL2_ND(t0); m0];
% 
%     x0_ad=ADIM(x0,m0);

    y0=prob.y0;

    Phi0=eye(length(y0));
    vPhi0=reshape(Phi0,[length(y0)^2,1]);

    [tt, yy]=ode78(@(t,y) TwBP_EL(t,y,u),[0 (tf-t0)/TU],[y0; vPhi0],odeopt);

    ToF=(tf-t0)/86400;          % [d]
    mf=yy(end,7)*m0;            % [kg]
    mp=(yy(1,7)-yy(end,7))*m0;  % [kg]

    ttd=tt*TU/86400; % [d]

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
    xtf_ad=ADIM([xtf; 1],m0);
    rvtf=xtf_ad(1:6);

    df=yy(end,1:6).'-rvtf;

    r=sqrt(sum(yy(:,1:3).^2,2));

    [TTcc,~,P]=MARGO_param(r);

    TT=(MU*LU)/(TU^2)*TTcc(1,:)*1e6; % [mN]
    II=LU/TU*TTcc(2,:)/g0;

    S=SwFun(tt,yy,epsilon);
    H=Hamil(tt,yy);

    prob.mf=mf;
    prob.mp=mp;
    prob.S=S;
    prob.H=H;
    prob.tt_ad=tt;
    prob.tt=ttd;
    prob.yy=yy;

    if output==1

        plot3D(t0,tt,yy,targ);
    
        fprintf('Departure date: %s (%.1f MJD2000)\n',cspice_et2utc(t0,'C',3),et2MJD2000(t0))
        fprintf('Arrival date: %s (%.1f MJD2000)\n',cspice_et2utc(tf,'C',3),et2MJD2000(tf))
        fprintf('ToF: %.3f days, Final Mass: %.3f kg, Propellant Mass: %.3f kg\n\n',ToF,mf,mp)
    
        fprintf('Position error %.3e km\nVelocity error %.3e km/s\n',norm(df(1:3))*LU, norm(df(4:6))*LU/TU)

        fprintf('Max relative Hamiltonian variation %.3e%%\n\n',abs((max(H)-min(H))/max(H))*100)
    

        figure

        %-throttle---------------------------------------------------------
        subplot(4,2,1)
        plot(ttd,1.*(S<0)+0)
        xlim([ttd(1) ttd(end)])
        ylim([0 1.1])
        ylabel('$u$')
        grid on
        grid minor

        %-switching-function-----------------------------------------------
        subplot(4,2,3)
        plot(ttd,S)
        xlim([ttd(1) ttd(end)])
        ylabel('$S$')
        grid on
        grid minor
            
        %-mass-------------------------------------------------------------
        subplot(4,2,5)
        plot(ttd,yy(:,7)*m0)
%         hold on
%         plot([ttd(1) ttd(end)],[yy(1,7) yy(end,7)]*m0,'--','color',[.8 .8 .8],'LineWidth',0.1)
        xlim([ttd(1) ttd(end)])
        ylabel('$m\,[kg]$')
        grid on
        grid minor

        %-input-power------------------------------------------------------
        subplot(4,2,2)
        plot(ttd,P)
        if max(P)>Pmax
            hold on
            plot([ttd(1) ttd(end)],[Pmax,Pmax])
        end
        if min(P)<Pmin
            hold on
            plot([ttd(1) ttd(end)],[Pmin,Pmin])
        end
        xlim([ttd(1) ttd(end)])
        ylabel('$P_{in}\,[W]$')
        grid on
        grid minor

        %-thrust-----------------------------------------------------------
        subplot(4,2,4)
        plot(ttd,II)
        xlim([ttd(1) ttd(end)])
        ylabel('$I_{sp}\,[s]$')
        grid on
        grid minor
    
        %-specific-impulse-------------------------------------------------
        subplot(4,2,6)
        plot(ttd,TT)
        xlim([ttd(1) ttd(end)])
        ylabel('$T_{max}\,[mN]$')
        grid on
        grid minor

        %-hamiltonian------------------------------------------------------
        subplot(4,2,8)
        plot(ttd,H)
        xlim([ttd(1) ttd(end)])
        ylabel('H')
        grid on
        grid minor
        xlabel('$t [days]$')
    end

    % User plot setting removal
    set(0,'DefaultAxesFontSize','remove');
    set(0,'DefaultLegendFontSize','remove');
    set(0,'DefaultLineLineWidth','remove');
    set(groot,'defaultTextInterpreter','remove');
    set(groot,'defaultLegendInterpreter','remove');
end