clc;

Y0=[prob.y0];
S0=zeros(size(prob));

for i=1:length(prob)
    S0(i)=prob(i).S(1);
end

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

isPareto=true(1,length(prob));

for i=1:length(prob)
    for j=1:length(prob)
        if i~=j
            if (prob(j).tf_ad <= prob(i).tf_ad && prob(j).mp < prob(i).mp) || (prob(j).tf_ad < prob(i).tf_ad && prob(j).mp <= prob(i).mp)
                isPareto(i) = false;
                break;
            end
        end
    end    
end

figure
plot(et2MJD2000([prob.t0]),[prob.tf_ad]*TU/86400,'linewidth',2)
hold on
plot(et2MJD2000([prob(isPareto).t0]),[prob(isPareto).tf_ad]*TU/86400,'xk','linewidth',2)
grid on
grid minor
axis tight
ylim([100 1100])
title('tf')

figure
plot([prob.tf_ad]*TU/86400,[prob.mp],'LineWidth',2)
hold on
plot([prob(isPareto).tf_ad]*TU/86400,[prob(isPareto).mp],'xk','LineWidth',2)
grid on
grid minor
axis tight
xlabel('tf')
ylabel('mp')

figure
plot(et2MJD2000([prob.t0]),[prob.tf_ad]/max([prob.tf_ad]),'r','linewidth',2)
hold on
plot(et2MJD2000([prob.t0]),[prob.mp]/max([prob.mp]),'b','linewidth',2)
grid on
grid minor
axis tight
xlabel('tf')
legend('tf','mp','location','best')

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

rrM=zeros(size(prob));
rrm=rrM;
for i=1:length(prob)
    rrM(i)=max(sqrt(sum(prob(i).zz(:,1:3).^2,2)));
    rrm(i)=min(sqrt(sum(prob(i).zz(:,1:3).^2,2)));
end

P=@(r) dot(r.^(0:4),[840.11  -1754.3 1625.01 -739.87  134.45]);
rmlim=fzero(@(x) P(x)-120,1);

figure
plot(et2MJD2000([prob.t0]),rrm)
hold on
plot(et2MJD2000([prob.t0]),rmlim*ones(size(prob)),'r--')
plot(et2MJD2000([prob.t0]),rrM)
plot(et2MJD2000([prob.t0]),1.7047*ones(size(prob)),'r--')
grid on
grid minor
axis tight
title('r bounds')