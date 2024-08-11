function [df,J] = TO_ZFP(llt,t0,sc_param,targ)

% testing branch
% adimensional guess input llt
% dimensional input t0,sc_param
% adimensional output df,J

    odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    l0_g=llt(1:7);
    tf_ad=llt(8);  % adimensional
    tf=tf_ad*TU+t0;

    x0=[SEL2_ND(t0); sc_param(3)];

    [x0_ad, sc_param_ad]=ADIM(x0,sc_param);

    T=sc_param_ad(1);

    y0=[x0_ad; l0_g];

    Phi0=eye(length(y0));
    vPhi0=reshape(Phi0,[length(y0)^2,1]);

    [~, yy]=ode78(@(t,y) TwBP_EL(t,y,u,sc_param_ad),[0 tf_ad],[y0; vPhi0],odeopt);

    yf=yy(end,:).';

    rrf=yf(1:3);
    vvf=yf(4:6);
    mf=yf(7);

    llrf=yf(8:10);
    llvf=yf(11:13);
    lmf=yf(14);

    rf=norm(rrf);
    lvf=norm(llvf);

    vPhif=yf(15:210);
    Phif=reshape(vPhif,[14,14]);

    ff=TwBP_EL(tf_ad,yf,u,sc_param_ad);
    ffx=ff(1:7);

    dmf=ffx(7);

    dllrf=ff(8:10);
    dllvf=ff(11:13);
    dlmf=ff(14);

    Hf=dot([llrf; llvf; lmf],ffx);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

    xtf_ad=ADIM([xtf; 1],sc_param);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    rtf=norm(rrtf);

    aatf=-rrtf./norm(rrtf)^3;
    daatf=3*dot(rrtf,vvtf)*rrtf/rtf^5-vvtf/rtf^3;

    Om=Hf+1-dot([llrf; llvf],[vvtf; aatf]);

    df=zeros(8,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;
    df(8)=Om;

    %-Jacobian-------------------------------------------------------------
    J=zeros(8,8);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

    Dff=zeros([7,7]);
    Dff(1:3,:)=Phif(4:6,8:14);
    Dff(4:6,:)=(3*(rrf*rrf.')./rf^5-eye(3)./rf^3)*Phif(1:3,8:14)+...
        -u*T/mf*(Phif(11:13,8:14)./lvf+...
        -(llvf*llvf.')*Phif(11:13,8:14)./lvf^3+...
        -llvf*Phif(7,8:14)./(lvf*mf));
   
    J(8,1:7)=[llrf; llvf; lmf].'*Dff+ffx.'*Phif(8:14,8:14)+...
        -vvtf.'*Phif(8:10,8:14)-aatf.'*Phif(11:13,8:14);

    
    %-time-dependency------------------------------------------------------
    J(1:6,8)=ffx(1:6)-[vvtf; aatf];   % dpsi/dtf
    J(7,8)=dlmf;

    dffx=zeros(7,1);
    dffx(1:3)=ffx(4:6);
    dffx(4:6)=3*dot(rrf,vvf)*rrf/rf^5-vvf/rf^3+...
        -u*T/(mf*lvf)*(dllvf-dot(dllvf,llvf)*llvf/lvf^2-dmf*llvf/mf);

    J(8,8)=dot([llrf; llvf; lmf],dffx)+dot(ffx,[dllrf; dllvf; dlmf])+...
        -dot(llrf,aatf)-dot(vvtf,dllrf)-dot(llvf,daatf)-dot(aatf,dllvf);

end

