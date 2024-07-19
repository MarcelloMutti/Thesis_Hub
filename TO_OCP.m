clear all; close all; clc;

cspice_furnsh('kernels\naif0012.tls'); %LSK leapseconds
cspice_furnsh('kernels\gm_de432.tpc'); %PCK gravitational parameters
cspice_furnsh('kernels\de432s.bsp');   %SPK solar system eph
cspice_furnsh('kernels\20155782_2000SG344.bsp') % 2000SG334 ephemeris

% cspice_furnsh('kernels\pck00010.tpc'); %PCK planetary/lunar const, orientations, rotations
% cspice_furnsh('de431.bsp');

% cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

%-empeheris check in selected time window----------------------------------
t_wo=cspice_str2et('2024-01-01 12:00 TDB');
t_wc=cspice_str2et('2026-12-31 12:00 TDB');

tt=linspace(t_wo,t_wc,500);

rr_E=cspice_spkpos('Earth',tt,'ECLIPJ2000','NONE','Sun');
rr_T=cspice_spkpos('20155782',tt,'ECLIPJ2000','NONE','Sun');

figure
plot3(rr_E(1,:),rr_E(2,:),rr_E(3,:))
hold on
plot3(rr_T(1,:),rr_T(2,:),rr_T(3,:))

%-shooting problem set up--------------------------------------------------

m0=22.3; % [kg]
t0=cspice_str2et('2024 JUN 10 12:00:00');


SC_param=[MARGO_param(1); m0]; % [kg*km/s^2, km/s, kg]

x0=[SEL2_ND(t0); m0];

[x0_ad, SC_param_ad]=ADIM(x0,SC_param);
l0_g=ACT(x0_ad,SC_param_ad);

y0_g=[x0_ad; l0_g];

cspice_kclear();