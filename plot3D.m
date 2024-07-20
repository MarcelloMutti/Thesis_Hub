function plot3D(yy)
    
    figure
    plot3(yy(:,1),yy(:,2),yy(:,3),'r')
    hold on
    plot3(yy(1,1),yy(1,2),yy(1,3),'ob')
    plot3(yy(end,1),yy(end,2),yy(end,3),'kx')
    plot3(0,0,0,'+k')
    grid on
    grid minor
    xlabel('x')
    ylabel('y')
    zlabel('z')
    legend('','SEL2','AST','Sun')

end

