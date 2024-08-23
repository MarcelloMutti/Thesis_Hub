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
%
% cspice_furnsh('kernels\pck00010.tpc'); %PCK planetary/lunar const, orientations, rotations
% cspice_furnsh('de431.bsp');

%-shooting-problem-set-up--------------------------------------------------

str_wo='2022-12-31 12:00 UTC';
str_wc='2024-12-31 12:00 UTC';

targ='3054374';
m0=22.3; % [kg]
Pmax=120; % [W]
Pmin=20;  % [W]
epsilon=0;

% 3702319  2014 YD    [654]
% 3550232  2010 UE51  [390]
% 3568303  2011 MD    [500]
% 3054374  2000 SG344 [536]
% 20478784 2012 UV136 [411]

% t0=MJD20002et(8400); % [mjd2000]
% ToF_g=550; %[d]

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

t_wo=cspice_str2et(str_wo); % [8400.5 mjd2000]
t_wc=cspice_str2et(str_wc); % [10226.5 mjd2000]

% t_wo=MJD20002et(8400.5);
% t_wc=MJD20002et(8410.5);

prob=struct_assembly(targ,[t_wo t_wc],m0,[Pmin Pmax],epsilon);

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

% % Single solution
% 
% SC_param=[MARGO_param(1); m0]; % [kg*km/s^2, km/s, kg]
% 
% l0_g=ACT(prob); 
% 
% tf_g=ToF_g*86400/TU;
% 
% fsopt=optimoptions('fsolve','Display','iter-detailed','SpecifyObjectiveGradient',true,'FunctionTolerance',1e-14);
% 
% [lltf_TO,~,~,~,jacFD]=fsolve(@(llt) TO_ZFP(llt,prob),[l0_g; tf_g],fsopt);
% 
% [~,~,prob]=TO_ZFP(lltf_TO,prob);
% 
% prob=DispRes(lltf_TO,prob);

% t0 Continuation
prob=TO_t0CONT(prob);

Y0=[prob.y0];
S0=zeros(size(prob));

for i=1:length(prob)
    S0(i)=prob(i).S(1);
end

figure
plot(et2MJD2000([prob.t0]),[prob.tf_ad]*TU/86400,'linewidth',2)
grid on
grid minor
axis tight
ylim([100 1100])
title('tf')

figure
subplot(3,2,1)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(1:3,:).^2)))
grid on
grid minor
axis tight

subplot(3,2,3)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(4:6,:).^2)))
grid on
grid minor
axis tight

subplot(3,2,5)
plot(et2MJD2000([prob.t0]),Y0(7,:))
grid on
grid minor
axis tight

subplot(3,2,2)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(8:10,:).^2)))
grid on
grid minor
axis tight

subplot(3,2,4)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(11:13,:).^2)))
grid on
grid minor
axis tight

subplot(3,2,6)
plot(et2MJD2000([prob.t0]),Y0(14,:))
grid on
grid minor
axis tight

figure
plot(et2MJD2000([prob.t0]),S0)
grid on
grid minor
axis tight
title('S0')

% figure
% plot(et2MJD2000([prob.t0]),max([prob.S],2))
% grid on
% grid minor

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