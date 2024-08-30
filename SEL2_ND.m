function [xx_SEL2] = SEL2_ND(t)
% v1 (master) SEL2 approximation

% dimensional output

    t_wo=cspice_str2et('2022-12-31 12:00 UTC');

    mu_S=cspice_bodvrd('Sun','GM',1);
    mu_E=cspice_bodvrd('Earth','GM',1);
    
    % eph
    xx_E=cspice_spkezr('Earth',t,'ECLIPJ2000','NONE','Sun');

%     % prop
%     if length(t)==1
%         if t==t_wo
%             xx_E=cspice_spkezr('Earth',t,'ECLIPJ2000','NONE','Sun');
%         else
%             odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);
%             [~,yy]=ode113(@(t,y) TBdyn(t,y),[t_wo t],cspice_spkezr('Earth',t_wo,'ECLIPJ2000','NONE','Sun'),odeopt);
%             xx_E=yy(end,:).';
%         end
%     else
%         odeopt=odeset('RelTol',1e-12,'AbsTol',1e-12);
%         [~,yy]=ode113(@(t,y) TBdyn(t,y),t,cspice_spkezr('Earth',t_wo,'ECLIPJ2000','NONE','Sun'),odeopt);
%         xx_E=yy.';
%     end

%     xx_SEL2=xx_E.*(1+(mu_E/(3*mu_S))^(1/3));

    xx_SEL2=cspice_spkezr('392',t,'ECLIPJ2000','NONE','Sun');
end

