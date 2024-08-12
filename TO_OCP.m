% v1 (branch)

clear all; close all; clc;

% cspice_furnsh('kernels\naif0012.tls'); %LSK leapseconds
% cspice_furnsh('kernels\gm_de432.tpc'); %PCK gravitational parameters
% cspice_furnsh('kernels\de432s.bsp');   %SPK solar system eph
% 
% cspice_furnsh('kernels\3702319_2014YD.bsp');    % 2014 YD ephemerides
% cspice_furnsh('kernels\3550232_2010UE51.bsp');  % 2010 UE51 ephemerides
% cspice_furnsh('kernels\3568303_2011MD.bsp')     % 2011 MD ephemerides
% cspice_furnsh('kernels\3054374_2000SG344.bsp'); % 2000 SG334 ephemerides
% cspice_furnsh('kernels\20478784_2012UV136.bsp');% 2012 UV136 ephemerides
%
% cspice_furnsh('kernels\pck00010.tpc'); %PCK planetary/lunar const, orientations, rotations
% cspice_furnsh('de431.bsp');

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

% %-empeheris check in selected time window----------------------------------
% 
% t_wo=cspice_str2et('2022-12-31 12:00 TDB');
% t_wc=cspice_str2et('2027-12-31 12:00 TDB');
% 
% tt=linspace(t_wo,t_wc,500);
% 
% xx_SEL2=SEL2_ND(tt);
% rr_2014YD=cspice_spkpos('3702319',tt,'ECLIPJ2000','NONE','Sun');
% rr_2010UE51=cspice_spkpos('3550232',tt,'ECLIPJ2000','NONE','Sun');
% rr_2011MD=cspice_spkpos('3568303',tt,'ECLIPJ2000','NONE','Sun');
% rr_2000SG344=cspice_spkpos('3054374',tt,'ECLIPJ2000','NONE','Sun');
% rr_2012UV136=cspice_spkpos('20478784',tt,'ECLIPJ2000','NONE','Sun');
% 
% figure
% plot3(xx_SEL2(1,:),xx_SEL2(2,:),xx_SEL2(3,:))
% hold on
% plot3(rr_2014YD(1,:),rr_2014YD(2,:),rr_2014YD(3,:))
% plot3(rr_2010UE51(1,:),rr_2010UE51(2,:),rr_2010UE51(3,:))
% plot3(rr_2011MD(1,:),rr_2011MD(2,:),rr_2011MD(3,:))
% plot3(rr_2000SG344(1,:),rr_2000SG344(2,:),rr_2000SG344(3,:))
% plot3(rr_2012UV136(1,:),rr_2012UV136(2,:),rr_2012UV136(3,:))
% legend('SEL2','2014 YD','2010 UE51','2011 MD','2000 SG344','2012 UV136')

%-shooting problem set up--------------------------------------------------

t_wo=cspice_str2et('2022-12-31 12:00 UTC'); % [8400.5 mjd2000]
t_wc=cspice_str2et('2024-12-31 12:00 UTC'); % [10226.5 mjd2000]

m0=22.3; % [kg]
t0=MJD20002et(8400); % [mjd2000]
targ='3550232';

ToF_g=500; %[d]

% 3702319  2014 YD    [683]
% 3550232  2010 UE51  [390]
% 3568303  2011 MD    [490]
% 3054374  2000 SG344 [535]
% 20478784 2012 UV136 [335]

%--------------------------------------------------------------------------

SC_param=[MARGO_param(1); m0]; % [kg*km/s^2, km/s, kg]

l0_g=ACT(t0,SC_param);

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

tf_g=ToF_g*86400/TU;

fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true);

% % Single solution
% [lltf_TO,~,~,~,jacFD]=fsolve(@(llt) TO_ZFP(llt,t0,SC_param,targ),[l0_g; tf_g],fsopt);
% 
% [df,jacAN]=TO_ZFP(lltf_TO,t0,SC_param,targ);
% 
% DispRes(lltf_TO,t0,SC_param,targ,0);

% % First solution through Tcont
% lltf=TO_TCONT(t0,SC_param,targ);
% 
% [df,jacAN]=TO_ZFP(lltf,t0,SC_param,targ);
% 
% DispRes(lltf,t0,SC_param,targ,0);

[tt0,lltf_M,mp,Hf,mS]=TO_t0CONT([t_wo t_wc],SC_param,targ);

figure
subplot(2,2,1)
plot(et2MJD2000(tt0),lltf_M(8,:)*TU/86400)
grid on
grid minor
xlim([min(et2MJD2000(tt0)) max(et2MJD2000(tt0))])
legend('tf','Location','best')

subplot(2,2,3)
plot(et2MJD2000(tt0),mp)
grid on
grid minor
xlim([min(et2MJD2000(tt0)) max(et2MJD2000(tt0))])
legend('mf','Location','best')

subplot(2,2,2)
plot(et2MJD2000(tt0),mS)
grid on
grid minor
xlim([min(et2MJD2000(tt0)) max(et2MJD2000(tt0))])
legend('max S','Location','best')

subplot(2,2,4)
plot(et2MJD2000(tt0),Hf)
grid on
grid minor
xlim([min(et2MJD2000(tt0)) max(et2MJD2000(tt0))])
legend('Hf','Location','best')

% %-AUX solution attempt-----------------------------------------------------
% 
% % works with gamma=1 @ eps=1, but to be generalized in continuation process
% 
% epsilon=1;
% 
% ll_AUX=fsolve(@(ll) AUX_ZFP(ll,[t0 tf_TO],SC_param,targ,epsilon),2*epsilon*ll_TO,fsopt);
% 
% DispRes(ll_AUX,[t0 tf_TO],SC_param,targ,epsilon);

% cspice_kclear();