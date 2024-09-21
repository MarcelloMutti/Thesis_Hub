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

    rf=norm(zzf(1:3));
    [~,~,Sp]=MARGO_param(rf);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    if SwFun(tf_ad,zzf,prob.isFO)+prob.epsilon<0
        utype='on';
    elseif SwFun(tf_ad,zzf,prob.isFO)-prob.epsilon>0 || (SwFun(tf_ad,zzf,prob.isFO)-prob.epsilon>=0 && prob.epsilon==0)
        utype='off';
    elseif abs(SwFun(tf_ad,zzf,prob.isFO))-prob.epsilon<=0 && prob.epsilon~=0
        utype='med';
    end

    df=FO_gamma(zzf,prob,Ptype,utype);
    J=FO_jacobian(zzf,prob,Ptype,utype);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

function df = FO_gamma(zf,prob,Ptype,utype)

%     LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
%     TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

%     u=1;

    m0=prob.m0;
    tf=prob.tf;
%     tf_ad=prob.tf_ad;
    targ=prob.targ;

    rrf=zf(1:3);
    vvf=zf(4:6);

%     llrf=zf(8:10);
%     llvf=zf(11:13);
    lmf=zf(14);

%     ff=TO_2BP_SEP(tf_ad,zf,Ptype);
%     ffx=ff(1:7);

%     Hf=dot([llrf; llvf; lmf],ffx);

    xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

    xtf_ad=ADIM([xtf; 1],m0);

    rrtf=xtf_ad(1:3);
    vvtf=xtf_ad(4:6);

%     aatf=-rrtf./norm(rrtf)^3;

%     Om=Hf+1-dot([llrf; llvf],[vvtf; aatf]);

    df=zeros(7,1);
    df(1:6)=[rrf; vvf]-[rrtf; vvtf];
    df(7)=lmf;
%     df(8)=Om;
end

function J = FO_jacobian(zf,prob,Ptype,utype)

%     u=1;

%     m0=prob.m0;
%     tf=prob.tf;
%     tf_ad=prob.tf_ad;
%     targ=prob.targ;
% 
%     rrf=zf(1:3);
%     vvf=zf(4:6);
%     mf=zf(7);
% 
%     llvf=zf(11:13);

%     rf=norm(rrf);
%     lvf=norm(llvf);

%     G=3*(rrf*rrf.')/rf^5-eye(3)/rf^3;

    vPhif=zf(15:210);
    Phif=reshape(vPhif,[14,14]);

%     ff=TO_2BP_SEP(tf_ad,zf,Ptype);
%     ffx=ff(1:7);

%     dmf=ffx(7);
% 
%     dllrf=ff(8:10);
%     dllvf=ff(11:13);
%     dlmf=ff(14);

%     xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');

%     xtf_ad=ADIM([xtf; 1],m0);

%     rrtf=xtf_ad(1:3);
%     vvtf=xtf_ad(4:6);

%     rtf=norm(rrtf);

%     aatf=-rrtf./norm(rrtf)^3;
%     daatf=3*dot(rrtf,vvtf)*rrtf/rtf^5-vvtf/rtf^3;

%     [Tcf,dTcf]=MARGO_param(rf);

%     T=Tcf(1);
%     c=Tcf(2);
% 
%     Tp=dTcf(1); % dT/dr
%     cp=dTcf(2); % dc/dr

    %-Jacobian-------------------------------------------------------------
    J=zeros(7,7);
    
    %-costate-depenndency--------------------------------------------------
    J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
    J(7,1:7)=Phif(14,8:14);         % phi_lm,ll

%     Dffx=zeros([7,7]);
% 
%     Dffx(1:3,:)=Phif(4:6,8:14);
%     Dffx(4:6,:)=G*Phif(1:3,8:14)-u*llvf/lvf*(Tp*rrf.'*Phif(1:3,8:14)/rf-T*Phif(7,8:14)/mf)/mf-u*T/mf*(eye(3)-(llvf*llvf.')/lvf^2)*Phif(11:13,8:14)/lvf;
%     Dffx(7,:)=-u/c^2*(c*Tp-cp*T)*rrf.'*Phif(1:3,8:14)/rf;
% 
%     J(8,1:7)=zf(8:14).'*Dffx+(ffx-[vvtf; aatf; 0]).'*Phif(8:14,8:14);
% 
%     
%     %-time-dependency------------------------------------------------------
%     J(1:6,8)=ffx(1:6)-[vvtf; aatf];   % dpsi/dtf
%     J(7,8)=dlmf;
% 
%     dffx=zeros(7,1);
% 
%     dffx(1:3)=ffx(4:6);
%     dffx(4:6)=3*(rrf.'*vvf)*rrf/rf^5-vvf/rf^3-u*llvf/lvf*(Tp/mf*rrf.'*vvf/rf-T*dmf/mf^2)-u*T/mf*(dllvf-llvf*(llvf.'*dllvf)/lvf^2)/lvf;
%     dffx(7)=-u/c^2*(c*Tp-cp*T)*(rrf.'*vvf)/rf;
% 
%     J(8,8)=zf(8:14).'*(dffx-[aatf; daatf; 0])+(ffx-[vvtf; aatf; 0]).'*[dllrf; dllvf; dlmf];

end

