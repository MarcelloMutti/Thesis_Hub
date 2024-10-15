function [ll0] = ACT(prob)
% v1 (master)

% dimensional input
    
    t0=prob.t0;
    m0=prob.m0;

    xx_SEL2=cspice_spkezr('392',t0,'ECLIPJ2000','NONE','Sun');
    
    x0d=[xx_SEL2; m0];
    
    x0=ADIM(x0d,m0);

    rr=x0(1:3);
    vv=x0(4:6);
    m=x0(7);

    r=norm(rr);
    v=norm(vv);
    
    a=rad2deg(unifrnd(-10, 10));
    da=rad2deg(unifrnd(-5, 5));
    b=rad2deg(unifrnd(-1, 1));
    db=rad2deg(unifrnd(-0.1, 0.1));

    S=unifrnd(-1, 0)+(prob.isFO-1);
    dS=unifrnd(-0.01, 0.01);
    
    lm=1;
%     lm=unifrnd(0,1);

%     T=sc_param_ad(1);
%     c=sc_param_ad(2);

    Tc=MARGO_param(r);
    T=Tc(1);
    c=Tc(2);
    
    u=1;

    uup=[cos(a)*cos(b);
        sin(a)*cos(b);
        sin(b)];

    duup=[-sin(a)*da*cos(b)-cos(a)*sin(b)*db;
         cos(a)*da*cos(b)-sin(a)*sin(b)*db;
         cos(b)*db];

    lv=-m/c*(S+lm-prob.isFO);

    dm=-u*T/c;
    dlm=-u*lv*T/m^2;

    dlv=-dm/c*(S+lm)-m/c*(dS+dlm);

    vu=vv/v;

    hh=cross(rr,vv);
    h=norm(hh);
    hu=hh/h;

    R=[vu, cross(hu,vu), hu];
    uu=R*uup;

    llv=-lv*uu;
%     llv=lv*uu;

    gr=-rr/r^3;
    hv=zeros(3,1);

%     G=3*(rr*rr.')/r^5-eye(3)/r^3;
    H=zeros(3,3);

%     dvv=gr+hv-u*T/m*llv/lv; % why tho??
    dvv=gr+hv;
    dv=dot(dvv,vv)/v;
    dvu=dvv/v-dv*vv/v^2;

    dhh=cross(rr,dvv);
    dh=dot(dhh,hh)/h;
    dhu=dhh/h-dh*hh/h^2;

    dR=[dvu, cross(dhu,vu)+cross(hu,dvu), dhu];
    duu=R*duup+dR*uup;

    dllv=-dlv*uu-lv*duu;
%     dllv=dlv*uu+lv*duu;
    llr=-dllv-H.'*llv;

    ll0=zeros(7,1);
    ll0(1:3)=llr;
    ll0(4:6)=llv;
    ll0(7)=lm;

end