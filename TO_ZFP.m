function [df] = TO_ZFP(llt,t0,sc_param,targ)

% testing branch
% dimensional input (but costates)
% adimensional output

    odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    l0_g=llt(1:7);
    tf=llt(8);  % dimensional

    x0=[SEL2_ND(t0); sc_param(3)];
%     x0=[cspice_spkezr('Earth',t0,'ECLIPJ2000','NONE','Sun'); sc_param(3)]; % testing

    [x0_ad, sc_param_ad]=ADIM(x0,sc_param);

    y0=[x0_ad; l0_g];

    [~, yy]=ode78(@(t,y) TwBP_EL(t,y,u,sc_param_ad),[0 (tf-t0)/TU],y0,odeopt);

    yf=yy(end,:).';

    rvf=yf(1:6);
    llrf=yf(8:10);
    llvf=yf(11:13);
    lmf=yf(14);

    ff=TwBP_EL((tf-t0)/TU,yf,u,sc_param_ad);
    ffx=ff(1:7);

    Hf=dot([llrf; llvf; lmf],ffx);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
%     xtf=cspice_spkezr('Venus',tf,'ECLIPJ2000','NONE','Sun'); % testing

    %--seems correct up to here--

%     rrtf=xtf(1:3);
%     vvtf=xtf(4:6);
%     aatf=-rrtf./norm(rrtf)^3;

    xtf_ad=ADIM([xtf; 1],sc_param);
    rvtf=xtf_ad(1:6);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);
    aatf=-rrtf./norm(rrtf)^3;

    Om=Hf+1-dot([llrf; llvf],[vvtf; aatf]);

    df=zeros(8,1);
    df(1:6)=rvf-rvtf;
    df(7)=lmf;
    df(8)=Om;   

end

