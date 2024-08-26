function dy = TBdyn(t,y)
    
     mu_S=cspice_bodvrd('Sun','GM',1);

    dy=zeros(6,1);
    dy(1:3)=y(4:6);
    dy(4:6)=-mu_S*y(1:3)/norm(y(1:3))^3;

end