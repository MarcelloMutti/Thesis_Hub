function df = TO_gamma(yf,prob)

%     LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
%     TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    u=1;

    m0=prob.m0;
    tf=prob.tf;
    tf_ad=prob.tf_ad;
    targ=prob.targ;

    rrf=yf(1:3);
    vvf=yf(4:6);

    llrf=yf(8:10);
    llvf=yf(11:13);
    lmf=yf(14);
    ff=TwBP_EL(tf_ad,yf,u);
    ffx=ff(1:7);

    Hf=dot([llrf; llvf; lmf],ffx);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

    xtf_ad=ADIM([xtf; 1],m0);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    aatf=-rrtf./norm(rrtf)^3;

    Om=Hf+1-dot([llrf; llvf],[vvtf; aatf]);

    df=zeros(8,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;
    df(8)=Om;
end