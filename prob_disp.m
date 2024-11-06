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
plot(et2MJD2000([prob.t0]),[prob.tf_ad]*TU/86400,'.','linewidth',2)
hold on
plot(et2MJD2000([prob(isPareto).t0]),[prob(isPareto).tf_ad]*TU/86400,'xk','linewidth',2)
grid on
grid minor
axis tight
ylim([100 1100])
xlabel('$$t_0\,[MJD2000]$$','Interpreter','latex')
ylabel('$$ToF\,[d]$$','Interpreter','latex')

figure
plot([prob.tf_ad]*TU/86400,[prob.mp],'.','LineWidth',2)
hold on
plot([prob(isPareto).tf_ad]*TU/86400,[prob(isPareto).mp],'xk','LineWidth',2)
grid on
grid minor
axis tight
xlabel('tf')
ylabel('mp')

figure
plot(et2MJD2000([prob.t0]),[prob.tf_ad]/max([prob.tf_ad]),'r.','linewidth',2)
hold on
plot(et2MJD2000([prob.t0]),[prob.mp]/max([prob.mp]),'b.','linewidth',2)
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

% tt=prob(1).tt;
% zz=prob(1).zz;
% 
% dt=tt(2:end)-tt(1:end-1);
% % dz=sqrt(sum((zz(2:end,:)-zz(1:end-1,:)).^2,2));
% dz=max(abs(zz(2:end,:)-zz(1:end-1,:)),[],2);
% 
% figure
% plot(dt,dz)
% grid on
% grid minor
% axis tight
% xlabel('step size')
% ylabel('step error')

xlin=et2MJD2000(linspace(min([prob.t0]),max([prob.t0]),length(prob)/10));
ylin=TU/86400*(linspace(min([prob.tf_ad]),max([prob.tf_ad]),length(prob)/10));

K=boundary([prob.t0].',[prob.tf_ad].',0.95); % specify manually

% K=[];
% for i=1:length(TO_ref)
%     di=find([prob.t0]==TO_ref(i).t0);
%     if i==1 || i==length(TO_ref)
%         K=[K, di];
%     else
%         K=[K, di(1), di(end)];
%     end
% end
% [~,oid]=sort(K);
% K=K(oid);

[X,Y]=meshgrid(xlin,ylin);
Z=griddata(et2MJD2000([prob.t0]),TU/86400*([prob.tf_ad]),[prob.mp],X,Y,'linear');
[in,on]=inpolygon(X,Y,et2MJD2000([prob(K).t0]),TU/86400*([prob(K).tf_ad]));
Z(~in)=NaN;

figure
contourf(X,Y,Z,min([prob.mp]):0.3:max([prob.mp]),'fill','on')
c=colorbar;
hold on
contour(X,Y,Z,[2.8 2.8],'k--','LineWidth',2)
plot(et2MJD2000([TO_ref.t0]),TU/86400*([TO_ref.tf_ad]),'r','LineWidth',2)
plot(et2MJD2000([prob(isPareto).t0]),[prob(isPareto).tf_ad]*TU/86400,'xk','linewidth',2)
xlabel('$$t_0\, [MJD2000]$$','Interpreter','latex')
ylabel('$$ToF\, [days]$$','Interpreter','latex')
c.Label.String='mp';


prob_TO=prob(strcmp({prob.sts},'TO'));
prob_SK=prob(strcmp({prob.sts},'skp'));

prob_TOEO=[prob_TO prob_SK];
[~,oid]=sort([prob_TOEO.t0]);
prob_TOEO=prob_TOEO(oid);

y0TOEO=[prob_TOEO.y0];
l0TOEO=y0TOEO(8:14,:);

prob_TOTO=TO_ref(1:(length(prob_TO)+length(prob_SK)));

MS=[];
for i=1:length(prob_TOTO)
    MS(i)=max(prob_TOTO(i).S);
end

y0TOTO=[prob_TOTO.y0];
l0TOTO=y0TOTO(8:14,:);

gamma=l0TOEO./l0TOTO;

llm0=y0TOEO(14,:);

figure
subplot(2,2,1)
plot(et2MJD2000([prob_TOEO.t0]),gamma.')
grid on
grid minor
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\gamma=\lambda^{EO}_0/\lambda^{TO}_0$$','Interpreter','latex')

subplot(2,2,3)
plot(et2MJD2000([prob_TOEO.t0]),llm0)
grid on
grid minor
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\lambda^{EO}_{m,0}$$','Interpreter','latex')

subplot(2,2,2)
plot(et2MJD2000([prob_TOEO.t0]),mean(gamma)./llm0)
grid on
grid minor
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\frac{\gamma}{\lambda^{EO}_{m,0}}$$','Interpreter','latex')

subplot(2,2,4)
plot(et2MJD2000([prob_TOEO.t0]),-1./MS)
grid on
grid minor
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$-1/\max(S)$$','Interpreter','latex')

figure
subplot(3,1,1)
semilogy(et2MJD2000([prob_TOEO.t0]),mean(gamma),'-ob')
hold on
semilogy(et2MJD2000([prob_TOEO.t0]),-1./MS,'--xr')
semilogy(et2MJD2000([prob_TOEO.t0]),-2./MS,'--k')
% semilogy(et2MJD2000([prob_TOEO.t0]),-1./MS*mean(mean(gamma)./(-1./MS)),'--k')
grid on
grid minor
xlim(et2MJD2000([min([prob_TOEO.t0]) max([prob_TOEO.t0])]))
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\gamma,\,-1/\max(S)$$','Interpreter','latex')

subplot(3,1,2)
semilogy(et2MJD2000([prob_TOEO.t0]),mean(gamma)./(-1./MS),'-ob')
grid on
grid minor
xlim(et2MJD2000([min([prob_TOEO.t0]) max([prob_TOEO.t0])]))
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\frac{\gamma}{-1/\max(S)}$$','Interpreter','latex')

subplot(3,1,3)
semilogy(et2MJD2000([prob_TOEO.t0]),mean(gamma)+1./MS,'-ob')
grid on
grid minor
xlim(et2MJD2000([min([prob_TOEO.t0]) max([prob_TOEO.t0])]))
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\gamma-(-1/\max(S))$$','Interpreter','latex')

figure
semilogy(et2MJD2000([prob_TOEO.t0]),mean(gamma),'-ob')
hold on
semilogy(et2MJD2000([prob_TOEO.t0]),-1./MS,'--xr')
grid on
grid minor
xlim(et2MJD2000([min([prob_TOEO.t0]) max([prob_TOEO.t0])]))
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\gamma,\,-1/\max(S)$$','Interpreter','latex')
legend('$$\gamma$$','$$\frac{-1}{\max(S)}$$','interpreter','latex');

figure
semilogy(et2MJD2000([prob_TOEO.t0]),-MS.*mean(gamma),'-ob')
grid on
grid minor
xlim(et2MJD2000([min([prob_TOEO.t0]) max([prob_TOEO.t0])]))
xlabel('$$t_0 \, [MDJ2000]$$','Interpreter','latex')
ylabel('$$\frac{\gamma}{-1/\max(S)}$$','Interpreter','latex')