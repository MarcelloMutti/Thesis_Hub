clear all; close all; clc;

% cspice_furnsh('kernels\naif0012.tls'); %LSK leapseconds
% cspice_furnsh('kernels\gm_de440.tpc'); %PCK gravitational parameters
% cspice_furnsh('kernels\de440.bsp');    %SPK solar system eph
% 
% cspice_furnsh('kernels\3702319_2014YD.bsp');    % 2014 YD ephemerides
% cspice_furnsh('kernels\3550232_2010UE51.bsp');  % 2010 UE51 ephemerides
% cspice_furnsh('kernels\3568303_2011MD.bsp')     % 2011 MD ephemerides
% cspice_furnsh('kernels\3054374_2000SG344.bsp'); % 2000 SG334 ephemerides
% cspice_furnsh('kernels\20478784_2012UV136.bsp');% 2012 UV136 ephemerides

%-shooting-problem-set-up--------------------------------------------------

str_wo='2022-12-31 12:00 UTC';
str_wc='2024-12-31 12:00 UTC';

targ='3702319';
m0=22.6;  % [kg]
Pmax=120; % [W]
Pmin=20;  % [W]
isFO=1;

% 3054374  2000 SG344 [538]
% 3550232  2010 UE51  [392]
% 3568303  2011 MD    [504]
% 20478784 2012 UV136 [428]
% 3702319  2014 YD    [682]

% t0=MJD20002et(8400); % [mjd2000]
% ToF_g=550; %[d]

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

t_wo=cspice_str2et(str_wo); % [8400.5 mjd2000]
t_wc=cspice_str2et(str_wc); % [9131.5 mjd2000]

FO_prob=struct_assembly(targ,[t_wo t_wc],m0,[Pmin Pmax],isFO);

% FO_SEP_sym_dyn(FO_prob); % to KEEP
addpath('FO_dyn')

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

% % Single solution
% FO_prob.t0=t_wo+0*86400;
% Tof=1095;
% % 391.845
% tf_ad=Tof*86400/TU;
% 
% FO_prob.tf_ad=tf_ad;
% FO_prob.tf=FO_prob.t0+tf_ad*TU;
% 
% xx0d=[cspice_spkezr('392',FO_prob.t0,'ECLIPJ2000','NONE','Sun'); FO_prob.m0];
% xx0=ADIM(xx0d, FO_prob.m0);
% 
% ll0=ACT(FO_prob);
% % ll0=-0.5*rand([7,1]).*rand([7,1]);
% % ll0=[6.703257711968917e-01
% %     -5.127666636081690e-01
% %      1.449497544940961e-02
% %      5.919420522718783e-01
% %      4.741306434873093e-02
% %     -3.671240124429167e-03
% %      1.000000000000000e+00];
% 
% FO_prob.y0=[xx0; ll0];
% 
% z0=[xx0; ll0; reshape(eye(14),[14*14,1])];
% 
% FO_prob.epsilon=1.0;
% 
% % [tt,zz]=FO_ode78(FO_prob,[0 tf_ad],z0);
% 
% fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'OptimalityTolerance',1e-8,'FunctionTolerance',1e-8,'MaxIterations',2e2);
% 
% [ll_FO,~,ex_flag]=fsolve(@(ll) FO_ZFP(ll,FO_prob),ll0,fsopt);
% 
% [~,~,FO_prob]=FO_ZFP(ll_FO,FO_prob);
% FO_prob=DispRes(FO_prob);

load('2014TOs_cp.mat');

EO_prob=FO_CONT(FO_prob,TO_prob);


