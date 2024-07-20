% v1 (branch)

clear all; close all; clc;

cspice_furnsh('kernels\naif0012.tls'); %LSK leapseconds
cspice_furnsh('kernels\gm_de432.tpc'); %PCK gravitational parameters
cspice_furnsh('kernels\de432s.bsp');   %SPK solar system eph

% cspice_furnsh('kernels\3702319_2014YD.bsp');    % 2014 YD ephemerides
% cspice_furnsh('kernels\3550232_2010UE51.bsp');  % 2010 UE51 ephemerides
% cspice_furnsh('kernels\3568303_2011MD.bsp')     % 2011 MD ephemerides
cspice_furnsh('kernels\3054374_2000SG344.bsp'); % 2000 SG334 ephemerides
% cspice_furnsh('kernels\20478784_2012UV136.bsp');% 2012 UV136 ephemerides

% cspice_furnsh('kernels\pck00010.tpc'); %PCK planetary/lunar const, orientations, rotations
% cspice_furnsh('de431.bsp');

% cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

%-empeheris check in selected time window----------------------------------

t_wo=cspice_str2et('2022-12-31 12:00 TDB');
t_wc=cspice_str2et('2027-12-31 12:00 TDB');

tt=linspace(t_wo,t_wc,500);

xx_SEL2=SEL2_ND(tt);
% rr_2014YD=cspice_spkpos('3702319',tt,'ECLIPJ2000','NONE','Sun');
% rr_2010UE51=cspice_spkpos('3550232',tt,'ECLIPJ2000','NONE','Sun');
% rr_2011MD=cspice_spkpos('3568303',tt,'ECLIPJ2000','NONE','Sun');
rr_2000SG344=cspice_spkpos('3054374',tt,'ECLIPJ2000','NONE','Sun');
% rr_2012UV136=cspice_spkpos('20478784',tt,'ECLIPJ2000','NONE','Sun');

figure
plot3(xx_SEL2(1,:),xx_SEL2(2,:),xx_SEL2(3,:))
hold on
% plot3(rr_2014YD(1,:),rr_2014YD(2,:),rr_2014YD(3,:))
% plot3(rr_2010UE51(1,:),rr_2010UE51(2,:),rr_2010UE51(3,:))
% plot3(rr_2011MD(1,:),rr_2011MD(2,:),rr_2011MD(3,:))
plot3(rr_2000SG344(1,:),rr_2000SG344(2,:),rr_2000SG344(3,:))
% plot3(rr_2012UV136(1,:),rr_2012UV136(2,:),rr_2012UV136(3,:))
% legend('SEL2','2014 YD','2010 UE51','2011 MD','2000 SG344','2012 UV136')

%-shooting problem set up--------------------------------------------------

m0=22.3; % [kg]
t0=cspice_str2et('2023 JUL 19 12:00:00 UTC');

SC_param=[MARGO_param(1); m0]; % [kg*km/s^2, km/s, kg]

x0=[SEL2_ND(t0); m0];

[x0_ad, SC_param_ad]=ADIM(x0,SC_param);
l0_g=ACT(x0_ad,SC_param_ad);

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

tf=cspice_str2et('2024 SEP 19 23:31:00 UTC');

fsopt=optimoptions('fsolve','Display','iter-detailed');

lltf=fsolve(@(llt) TO_ZFP(llt,t0,SC_param),[l0_g; tf],fsopt);

disp([l0_g lltf(1:7)])
disp([tf lltf(8)])
disp(TO_ZFP(lltf,t0,SC_param))


plots(lltf,t0,SC_param);

cspice_kclear();