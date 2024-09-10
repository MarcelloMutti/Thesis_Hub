% v1 (branch)

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

targ='3054374';
m0=22.6;  % [kg]
Pmax=120; % [W]
Pmin=20;  % [W]
epsilon=0;

% 3054374  2000 SG344 [538] (+2)
% 3550232  2010 UE51  [392] (+2)
% 3568303  2011 MD    [504] (+2)
% 20478784 2012 UV136 [428] (+14)
% 3702319  2014 YD    [682] (+30)

% t0=MJD20002et(8400); % [mjd2000]
% ToF_g=550; %[d]

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

t_wo=cspice_str2et(str_wo); % [8400.5 mjd2000]
t_wc=cspice_str2et(str_wc); % [9131.5 mjd2000]

prob=struct_assembly(targ,[t_wo t_wc],m0,[Pmin Pmax],epsilon);

% TO_SEP_sym_dyn(prob); % to KEEP

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

% t0 Continuation
prob=TO_t0CONT(prob);

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