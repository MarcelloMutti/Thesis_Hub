Y0=[prob.y0];
S0=zeros(size(prob));

for i=1:length(prob)
    S0(i)=prob(i).S(1);
end

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

figure
plot(et2MJD2000([prob.t0]),[prob.tf_ad]*TU/86400,'linewidth',2)
grid on
grid minor
axis tight
% ylim([100 1100])
title('tf')

figure
subplot(3,2,1)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(1:3,:).^2)))
grid on
grid minor
axis tight
ylabel('r')

subplot(3,2,3)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(4:6,:).^2)))
grid on
grid minor
axis tight
ylabel('v')

subplot(3,2,5)
plot(et2MJD2000([prob.t0]),Y0(7,:))
grid on
grid minor
axis tight
ylabel('m')

subplot(3,2,2)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(8:10,:).^2)))
grid on
grid minor
axis tight
ylabel('lr')

subplot(3,2,4)
plot(et2MJD2000([prob.t0]),sqrt(sum(Y0(11:13,:).^2)))
grid on
grid minor
axis tight
ylabel('lv')

subplot(3,2,6)
plot(et2MJD2000([prob.t0]),Y0(14,:))
grid on
grid minor
axis tight
ylabel('lm')

figure
plot(et2MJD2000([prob.t0]),S0)
grid on
grid minor
axis tight
ylabel('S0')

figure
subplot(3,2,1)
plot(et2MJD2000([prob.t0]),Y0(8,:))
grid on
grid minor
axis tight
ylabel('lrx')

subplot(3,2,3)
plot(et2MJD2000([prob.t0]),Y0(9,:))
grid on
grid minor
axis tight
ylabel('lry')

subplot(3,2,5)
plot(et2MJD2000([prob.t0]),Y0(10,:))
grid on
grid minor
axis tight
ylabel('lrz')

subplot(3,2,2)
plot(et2MJD2000([prob.t0]),Y0(11,:))
grid on
grid minor
axis tight
ylabel('lvx')

subplot(3,2,4)
plot(et2MJD2000([prob.t0]),Y0(12,:))
grid on
grid minor
axis tight
ylabel('lvy')

subplot(3,2,6)
plot(et2MJD2000([prob.t0]),Y0(13,:))
grid on
grid minor
axis tight
ylabel('lvz')

figure
plot(et2MJD2000([prob.t0]),Y0(14,:))
grid on
grid minor
axis tight
ylabel('lm0')

DF=[prob.gamma];

DR=sqrt(sum(DF(1:3,:).^2))*LU;
DV=sqrt(sum(DF(4:6,:).^2))*LU/TU;

figure
subplot(1,2,1)
plot(et2MJD2000([prob.t0]),DR)
grid on
grid minor
axis tight
title('r error')

subplot(1,2,2)
plot(et2MJD2000([prob.t0]),DV)
grid on
grid minor
axis tight
title('v error')

[tP,tF]=plomb([prob.tf_ad],([prob.t0]));
[vP,vF]=plomb(sqrt(sum(Y0(4:6,:).^2)),([prob.t0]));

figure
subplot(1,2,1)
plot(tF,tP)
grid on
grid minor
axis tight
title('tf spectra')

subplot(1,2,2)
plot(vF,vP)
grid on
grid minor
axis tight
title('v spectra')