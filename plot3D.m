function plot3D(t0,tt,yy,targ)
    
    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    ttd=tt.'.*TU+t0;

    xx_SEL2=SEL2_ND(ttd);
    rr_SEL2=xx_SEL2(1:3,:)./LU;

    rrt=cspice_spkpos(targ,ttd,'ECLIPJ2000','NONE','Sun')./LU;
    
    figure
    plot3(yy(:,1),yy(:,2),yy(:,3),'r')
    view([55, 55])
    hold on
    plot3(yy(1,1),yy(1,2),yy(1,3),'ob')
    plot3(yy(end,1),yy(end,2),yy(end,3),'kx')
    plot3(0,0,0,'+k')
    plot3(rr_SEL2(1,:),rr_SEL2(2,:),rr_SEL2(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
    plot3(rrt(1,:),rrt(2,:),rrt(3,:),'color',[.8 .8 .8],'LineWidth',0.1)
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
    legend('','SEL2','AST','Sun')

end

