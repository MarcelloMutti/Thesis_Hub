function [df,J,prob] = FO_ZFP(ll0,prob)

% testing branch
% adimensional guess input llt
% dimensional input t0,sc_param
% adimensional output df,J

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    m0=prob.m0;
    t0=prob.t0;

%     ll0=ll(1:7);
%     tf_ad=ll(8);  % adimensional
%     tf=tf_ad*TU+t0;

    tf_ad=prob.tf_ad;
%     tf=prob.tf;

    xx_SEL2=cspice_spkezr('392',t0,'ECLIPJ2000','NONE','Sun');

    xx0=[xx_SEL2; m0];

    xx0_ad=ADIM(xx0,m0);

    yy0=[xx0_ad; ll0];

    Phi0=eye(length(yy0));
    vPhi0=reshape(Phi0,[length(yy0)^2,1]);

%     [~,zz] =FO_ode78(prob,[0 tf_ad],[yy0; vPhi0]);
    [~,zz] =FO_ode87(prob,[0 tf_ad],[yy0; vPhi0]);

    zzf=zz(end,:).';

    df=FO_gamma(zzf,prob);
    J=FO_jacobian(zzf);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

function df = FO_gamma(zf,prob)

    m0=prob.m0;
    tf=prob.tf;
    targ=prob.targ;

    rrf=zf(1:3);
    vvf=zf(4:6);

    lmf=zf(14);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

    xtf_ad=ADIM([xtf; 1],m0);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

    df=zeros(7,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;

end

function J = FO_jacobian(zf)

    vPhif=zf(15:210);
    Phif=reshape(vPhif,[14,14]);

    %-Jacobian-initialization----------------------------------------------
    J=zeros(7,7);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

end

