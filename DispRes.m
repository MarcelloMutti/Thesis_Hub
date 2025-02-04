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

    Pmin=prob.Plim(1);
    Pmax=prob.Plim(2);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1
    MU=prob.m0;

    tf_ad=prob.tf_ad;
    tf=prob.tf;

    y0=prob.y0;

    Phi0=eye(length(y0));
    vPhi0=reshape(Phi0,[length(y0)^2,1]);

    if prob.isFO==1

        [tt,zz]=FO_ode87(prob,[0 tf_ad],[y0; vPhi0]);

    else

        [tt,zz]=TO_ode87(prob,[0 tf_ad],[y0; vPhi0]);

    end

    ToF=(tf-t0)/86400;          % [d]
    mf=zz(end,7)*m0;            % [kg]
    mp=(zz(1,7)-zz(end,7))*m0;  % [kg]

    ttd=tt*TU/86400; % [d]

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
    xtf_ad=ADIM([xtf; 1],m0);
    rvtf=xtf_ad(1:6);

    df=zz(end,1:6).'-rvtf;
    % df2=FO_ZFP(prob.y0(8:14),prob);

    r=sqrt(sum(zz(:,1:3).^2,2));

    [TTcc,~,P]=MARGO_param(r);

    TT=(MU*LU)/(TU^2)*TTcc(1,:)*1e6; % [mN]
    II=LU/TU*TTcc(2,:)/g0;           % [s]

    ep=prob.epsilon;
    Se=SwFun(tt,zz,prob.isFO);

    if ep>0

        u=1.*(Se<-ep)+(ep-Se)./(2*ep).*(abs(Se)<=ep)+0;
        % u_d=u;
        % ttd=ttd;
        % zz=zz;

    else
        
        u=zeros(size(Se));
        ut=1.*(Se(1)<ep)+0;
        u(1)=ut;

        for j=2:length(Se)
            if ttd(j)==ttd(j-1) && isapprox(Se(j),ep,AbsoluteTolerance=1e-10)
                u(j)=1-ut;
                ut=u(j);
            else
                u(j)=ut;
            end
        end

    end

    % % only for ep=0
    % if ep==0
    %     s=1;
    %     u0=u(s);
    % 
    %     ttd_d=[];   % used for bang-bang display
    %     u_d=[];     % used for bang-bang display
    %     zz_d=[];
    % 
    %     idu=0;
    % 
    %     while ~isempty(idu)
    %         idu=find(u(s:end)~=u0,1,'first')+s-1;
    % 
    %         if ~isempty(idu)
    % 
    %             u_d=[u_d; u(s:idu-1); u(idu-u0)];
    %             ttd_d=[ttd_d; ttd(s:idu-1); ttd(idu-1+u0)];
    %             zz_d=[zz_d; zz(s:idu-1,:); zz(idu-1+u0,:)];
    % 
    %             u0=u(idu);
    %             s=idu;
    % 
    %         else
    % 
    %             u_d=[u_d; u(s:end)];
    %             ttd_d=[ttd_d; ttd(s:end)];
    %             zz_d=[zz_d; zz(s:end,:)];
    % 
    %         end
    % 
    %     end
    % else
    %     u_d=u;
    %     ttd_d=ttd;
    %     zz_d=zz;
    % end

    H=Hamil(tt,zz,u,prob);

    prob.mf=mf;
    prob.mp=mp;
    prob.S=Se;
    prob.H=H;
    prob.tt_ad=tt;
    prob.tt=ttd;
    prob.zz=zz;

    if output==1

        plot3D(t0,ttd,zz,targ,u);
    
        fprintf('Departure date: %s (%.1f MJD2000)\n',cspice_et2utc(t0,'C',0),et2MJD2000(t0))
        fprintf('Arrival date: %s (%.1f MJD2000)\n',cspice_et2utc(tf,'C',0),et2MJD2000(tf))
        fprintf('ToF: %.3f days, Final Mass: %.3f kg, Propellant Mass: %.3f kg\n\n',ToF,mf,mp)
    
        fprintf('Position error %.3e km\nVelocity error %.3e km/s\n',norm(df(1:3))*LU, norm(df(4:6))*LU/TU)

        fprintf('Max relative Hamiltonian variation %.3e%%\n\n',abs((max(H)-min(H))/max(H))*100)
    

        figure

        %-throttle---------------------------------------------------------
        subplot(4,2,3)
        plot(ttd,u)
        axis tight
        ylim([0 1.1])
        ylabel('$u$')
        grid on
        grid minor

        %-switching-function-----------------------------------------------
        subplot(4,2,5)
        plot(ttd,Se)
        axis tight
        ylabel('$S$')
        grid on
        grid minor
            
        %-mass-------------------------------------------------------------
        subplot(4,2,7)
        plot(ttd,zz(:,7)*m0)
        axis tight
        ylabel('$m\,[kg]$')
        grid on
        grid minor

        %-input-power------------------------------------------------------
        subplot(4,2,4)
        plot(ttd,P)
        if max(P)>Pmax
            hold on
            plot([ttd(1) ttd(end)],[Pmax,Pmax],'r')
        end
        if min(P)<Pmin
            hold on
            plot([ttd(1) ttd(end)],[Pmin,Pmin],'r')
        end
        axis tight
        ylabel('$P_{in}\,[W]$')
        grid on
        grid minor

        %-thrust-----------------------------------------------------------
        subplot(4,2,6)
        plot(ttd,II)
        axis tight
        ylabel('$I_{sp}\,[s]$')
        grid on
        grid minor
    
        %-specific-impulse-------------------------------------------------
        subplot(4,2,8)
        plot(ttd,TT)
        axis tight
        ylabel('$T_{max}\,[mN]$')
        grid on
        grid minor

        %-hamiltonian------------------------------------------------------
        subplot(4,2,1)
        plot(ttd,H)
        axis tight
        ylabel('H')
        grid on
        grid minor
        
        subplot(4,2,2)
        plot(ttd,H)
        axis tight
        ylabel('H')
        grid on
        grid minor

        figure

        %-throttle---------------------------------------------------------
        subplot(3,2,1)
        plot(ttd,u)
        axis tight
        ylim([0 1.1])
        ylabel('$u$')
        grid on
        grid minor

        %-switching-function-----------------------------------------------
        subplot(3,2,3)
        plot(ttd,Se)
        axis tight
        ylabel('$S$')
        grid on
        grid minor
            
        %-mass-------------------------------------------------------------
        subplot(3,2,5)
        plot(ttd,zz(:,7)*m0)
        axis tight
        ylabel('$m\,[kg]$')
        grid on
        grid minor

        %-input-power------------------------------------------------------
        subplot(3,2,2)
        plot(ttd,P)
        if max(P)>Pmax
            hold on
            plot([ttd(1) ttd(end)],[Pmax,Pmax],'--r','LineWidth',1)
        end
        if min(P)<Pmin
            hold on
            plot([ttd(1) ttd(end)],[Pmin,Pmin],'--r','LineWidth',1)
        end
        axis tight
        ylabel('$P_{in}\,[W]$')
        grid on
        grid minor

        %-thrust-----------------------------------------------------------
        subplot(3,2,4)
        plot(ttd,II)
        axis tight
        ylabel('$I_{sp}\,[s]$')
        grid on
        grid minor
    
        %-specific-impulse-------------------------------------------------
        subplot(3,2,6)
        plot(ttd,TT)
        axis tight
        ylabel('$T_{max}\,[mN]$')
        grid on
        grid minor
        
    end



    % User plot setting removal
    set(0,'DefaultAxesFontSize','remove');
    set(0,'DefaultLegendFontSize','remove');
    set(0,'DefaultLineLineWidth','remove');
    set(groot,'defaultTextInterpreter','remove');
    set(groot,'defaultLegendInterpreter','remove');
end

function plot3D(t0,tt,zz,targ,u)
    
    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    ttd=tt.'.*86400+t0;

    xx_SEL2=cspice_spkezr('392',ttd,'ECLIPJ2000','NONE','Sun');
    rr_SEL2=xx_SEL2(1:3,:)./LU;

    rrt=cspice_spkpos(targ,ttd,'ECLIPJ2000','NONE','Sun')./LU;

%     rr_on=zz(S<0,1:3);
%     rr_off=zz(S>=0,1:3);

    % r=sqrt(sum(zz(:,1:3).^2,2));
    % Tc=MARGO_param(r);
    % T=Tc(1,:);
    % c=Tc(2,:);

    L=length(u);
    
    figure
    plot3(rr_SEL2(1,1),rr_SEL2(2,1),rr_SEL2(3,1),'ob')
    hold on
    plot3(rrt(1,end),rrt(2,end),rrt(3,end),'kx')
    plot3(0,0,0,'+k')
    plot3(rr_SEL2(1,:),rr_SEL2(2,:),rr_SEL2(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
    plot3(rrt(1,:),rrt(2,:),rrt(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
    for i=1:L-1
        line([zz(i,1), zz(i+1,1)],[zz(i,2), zz(i+1,2)],[zz(i,3), zz(i+1,3)],'color',[u(i) 0 1-u(i)],'linewidth',u(i)+1);
    end
    % for i=1:L
    %     if u(i)~=0
    %         th=-u(i)*zz(i,11:13)/norm(zz(i,11:13))*T(i)/zz(i,7);
    %         quiver3(zz(i,1),zz(i,2),zz(i,3),th(1),th(2),th(3),'color',[u(i) 0 1-u(i)])
    %     end
    % end
    view([55, 55])    
    xlim([-1.5 1.5])
    ylim([-1.5 1.5])

%     sp=get(gca,'DataAspectRatio');
%     if sp(3)==1
%           set(gca,'DataAspectRatio',[1 1 1/max(sp(1:2))])
%     else
%           set(gca,'DataAspectRatio',[1 1 sp(3)])
%     end

    grid on
    grid minor
 
    xlabel('$x [AU]$')
    ylabel('$y [AU]$')
    zlabel('$z [AU]$')
    legend('SEL2','AST','Sun')

end