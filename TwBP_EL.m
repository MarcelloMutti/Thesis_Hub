function [dy] = TwBP_EL(~, y, u, sc_param)
% v1 (master)

% adimensional dynamics

    rr=y(1:3);
    vv=y(4:6);
    m=y(7);

    llr=y(8:10);
    llv=y(11:13);

    r=norm(rr);
    lv=norm(llv);

    gr=-rr/r^3;
    hv=zeros(3,1);

    G=3*(rr*rr.')/r^5-eye(3)/r^3;
    H=zeros(3,3);

    T=sc_param(1);
    c=sc_param(2);

    dy=zeros(14,1);

    dy(1:3)=vv;
    dy(4:6)=gr+hv-u*T/m*llv/lv;
    dy(7)=-u*T/c;

    dy(8:10)=-G.'*llv;
    dy(11:13)=-H.'*llv-llr;
    dy(14)=-u*lv*T/m^2;

end

