function J = TO_jacobian(yf,prob)

    u=1;

    m0=prob.m0;
    tf=prob.tf;
    tf_ad=prob.tf_ad;
    targ=prob.targ;

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

    ff=TwBP_EL(tf_ad,yf,u);
    ffx=ff(1:7);

    dmf=ffx(7);

    dllrf=ff(8:10);
    dllvf=ff(11:13);
    dlmf=ff(14);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

    xtf_ad=ADIM([xtf; 1],m0);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    rtf=norm(rrtf);

    aatf=-rrtf./norm(rrtf)^3;
    daatf=3*dot(rrtf,vvtf)*rrtf/rtf^5-vvtf/rtf^3;

    [Tcf,dTcf]=MARGO_param(rf);

    T=Tcf(1);
    c=Tcf(2);

    Tp=dTcf(1); % dT/dr
    cp=dTcf(2); % dc/dr

    %-Jacobian-------------------------------------------------------------
    J=zeros(8,8);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

    Dff=zeros([7,7]);
    Dff(1:3,:)=Phif(4:6,8:14);
    Dff(4:6,:)=(3*(rrf*rrf.')./rf^5-eye(3)./rf^3)*Phif(1:3,8:14)+...
        -u*T/(lvf*mf)*((eye(3)-(llvf*llvf.')./lvf^2)*Phif(11:13,8:14)-llvf*Phif(7,8:14)./mf)+...
        -u/mf*llvf/lvf*Tp*rrf.'*Phif(1:3,8:14)./rf;
    Dff(7,:)=-u*(c*Tp-T*cp)./c^2*rrf.'*Phif(1:3,8:14)./rf;
   
    J(8,1:7)=[llrf; llvf; lmf].'*Dff+ffx.'*Phif(8:14,8:14)+...
        -vvtf.'*Phif(8:10,8:14)-aatf.'*Phif(11:13,8:14);

    
    %-time-dependency------------------------------------------------------
    J(1:6,8)=ffx(1:6)-[vvtf; aatf];   % dpsi/dtf
    J(7,8)=dlmf;

    dffx=zeros(7,1);
    dffx(1:3)=ffx(4:6);
    dffx(4:6)=3*dot(rrf,vvf)*rrf/rf^5-vvf/rf^3+...
        -u*T/(mf*lvf)*(dllvf-dot(dllvf,llvf)*llvf/lvf^2-dmf*llvf/mf)+...
        -u*llvf/lvf*Tp/mf*dot(rrf,vvf)/rf;
    dffx(7)=-u*(c*Tp-T*cp)./c^2*dot(rrf,vvf)/rf;

    J(8,8)=dot([llrf; llvf; lmf],dffx)+dot(ffx,[dllrf; dllvf; dlmf])+...
        -dot(llrf,aatf)-dot(vvtf,dllrf)-dot(llvf,daatf)-dot(aatf,dllvf);

end