function [df,J,prob] = TO_ZFP(llt,prob)

% testing branch
% adimensional guess input llt
% dimensional input t0,sc_param
% adimensional output df,J



    odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);

    LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
    TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

%     u=1;

    m0=prob.m0;
    t0=prob.t0;
%     targ=prob.targ;

    ll0=llt(1:7);
    tf_ad=llt(8);  % adimensional
    tf=tf_ad*TU+t0;

    prob.tf_ad=tf_ad;
    prob.tf=tf;

    xx0=[SEL2_ND(t0); m0];

    xx0_ad=ADIM(xx0,m0);

%     T=sc_param_ad(1);

    yy0=[xx0_ad; ll0];

    Phi0=eye(length(yy0));
    vPhi0=reshape(Phi0,[length(yy0)^2,1]);

%     [~, yy]=ode78(@(t,y) TwBP_EL(t,y),[0 tf_ad],[yy0; vPhi0],odeopt);
    [~,yy] =ode78_cust(prob,[0 tf_ad],[yy0; vPhi0]);

%     options.AbsTol=1e-6;
%     options.RelTol=1e-6;
%     [~, yy]=rk78_cust([0 tf_ad],[yy0; vPhi0],options);

    yyf=yy(end,:).';

%     rrf=yf(1:3);
%     vvf=yf(4:6);
%     mf=yf(7);
% 
%     llrf=yf(8:10);
%     llvf=yf(11:13);
%     lmf=yf(14);
% 
%     rf=norm(rrf);
%     lvf=norm(llvf);
% 
%     vPhif=yf(15:210);
%     Phif=reshape(vPhif,[14,14]);
% 
%     ff=TwBP_EL(tf_ad,yf,u);
%     ffx=ff(1:7);
% 
%     dmf=ffx(7);
% 
%     dllrf=ff(8:10);
%     dllvf=ff(11:13);
%     dlmf=ff(14);
% 
%     Hf=dot([llrf; llvf; lmf],ffx);
% 
%     xtf=cspice_spkezr(targ,tf,'ECLIPJ2000','NONE','Sun');
% 
%     xtf_ad=ADIM([xtf; 1],m0);
% 
%     rrtf=xtf_ad(1:3);
%     vvtf=xtf_ad(4:6);
% 
%     rtf=norm(rrtf);
% 
%     aatf=-rrtf./norm(rrtf)^3;
%     daatf=3*dot(rrtf,vvtf)*rrtf/rtf^5-vvtf/rtf^3;
% 
%     Om=Hf+1-dot([llrf; llvf],[vvtf; aatf]);
% 
% 
%     [Tcf,dTcf]=MARGO_param(rf);
% 
%     T=Tcf(1)*M;
%     c=Tcf(2)*M;
% 
%     Tp=dTcf(1)*M; % dT/dr
%     cp=dTcf(2)*M; % dc/dr
% 
%     df=zeros(8,1);
%     df(1:6)=[rrf; vvf]-[rrtf; vvtf];
%     df(7)=lmf;
%     df(8)=Om;
% 
%     %-Jacobian-------------------------------------------------------------
%     J=zeros(8,8);
%     
%     %-costate-depenndency--------------------------------------------------
%     J(1:6,1:7)=Phif(1:6,8:14);      % phi_r,ll phi_v,ll
%     J(7,1:7)=Phif(14,8:14);         % phi_lm,ll
% 
%     Dff=zeros([7,7]);
%     Dff(1:3,:)=Phif(4:6,8:14);
%     Dff(4:6,:)=(3*(rrf*rrf.')./rf^5-eye(3)./rf^3)*Phif(1:3,8:14)+...
%         -u*T/(lvf*mf)*((eye(3)-(llvf*llvf.')./lvf^2)*Phif(11:13,8:14)-llvf*Phif(7,8:14)./mf)+...
%         -u/mf*llvf/lvf*Tp*rrf.'*Phif(1:3,8:14)./rf;
%     Dff(7,:)=-u*(c*Tp-T*cp)./c^2*rrf.'*Phif(1:3,8:14)./rf;
%    
%     J(8,1:7)=[llrf; llvf; lmf].'*Dff+ffx.'*Phif(8:14,8:14)+...
%         -vvtf.'*Phif(8:10,8:14)-aatf.'*Phif(11:13,8:14);
% 
%     
%     %-time-dependency------------------------------------------------------
%     J(1:6,8)=ffx(1:6)-[vvtf; aatf];   % dpsi/dtf
%     J(7,8)=dlmf;
% 
%     dffx=zeros(7,1);
%     dffx(1:3)=ffx(4:6);
%     dffx(4:6)=3*dot(rrf,vvf)*rrf/rf^5-vvf/rf^3+...
%         -u*T/(mf*lvf)*(dllvf-dot(dllvf,llvf)*llvf/lvf^2-dmf*llvf/mf)+...
%         -u*llvf/lvf*Tp/mf*dot(rrf,vvf)/rf;
%     dffx(7)=-u*(c*Tp-T*cp)./c^2*dot(rrf,vvf)/rf;
% 
%     J(8,8)=dot([llrf; llvf; lmf],dffx)+dot(ffx,[dllrf; dllvf; dlmf])+...
%         -dot(llrf,aatf)-dot(vvtf,dllrf)-dot(llvf,daatf)-dot(aatf,dllvf);
    
    rf=norm(yyf(1:3));
    [~,~,Sp]=MARGO_param(rf);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    df=TO_gamma(yyf,prob,Ptype);
    J=TO_jacobian(yyf,prob,Ptype);

    prob.y0=yy0;
    prob.gamma=df;
    prob.jac=J;

end

