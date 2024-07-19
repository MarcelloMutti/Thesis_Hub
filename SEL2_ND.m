function [xx_SEL2] = SEL2_ND(t)
% v1 (master) SEL2 approximation
% dimensional output

    mu_S=cspice_bodvrd('Sun','GM',1);
    mu_E=cspice_bodvrd('Earth','GM',1);

    xx_E=cspice_spkezr('Earth',t,'ECLIPJ2000','NONE','Sun');

    xx_SEL2=xx_E.*(1+(mu_E/(3*mu_S))^(1/3));
end

