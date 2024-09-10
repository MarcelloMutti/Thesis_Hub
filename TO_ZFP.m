function [df,J,prob] = TO_ZFP(llt,prob)

% testing branch
% adimensional guess input llt
% dimensional input t0,sc_param
% adimensional output df,J

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

    m0=prob.m0;
    t0=prob.t0;

    ll0=llt(1:7);
    tf_ad=llt(8);  % adimensional
    tf=tf_ad*TU+t0;

    prob.tf_ad=tf_ad;
    prob.tf=tf;

    xx_SEL2=cspice_spkezr('392',t0,'ECLIPJ2000','NONE','Sun');

    xx0=[xx_SEL2; m0];

    xx0_ad=ADIM(xx0,m0);

    yy0=[xx0_ad; ll0];

    Phi0=eye(length(yy0));
    vPhi0=reshape(Phi0,[length(yy0)^2,1]);

    [~,zz] =ode78_cust(prob,[0 tf_ad],[yy0; vPhi0]);


    zzf=zz(end,:).';

    rf=norm(zzf(1:3));
    [~,~,Sp]=MARGO_param(rf);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    df=TO_gamma(zzf,prob,Ptype);
    J=TO_jacobian(zzf,prob,Ptype);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

