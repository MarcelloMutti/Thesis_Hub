function [tt,zz]=TO_ode78(prob,tspan,z0)

    options=odeset('AbsTol',1e-12,'RelTol',1e-12,'Events',@(t,y) crossings(t,y,prob));

    r0=norm(z0(1:3));
    [~,~,Sp]=MARGO_param(r0);

    if Sp<prob.Plim(2)
        Ptype='med';
    elseif Sp>=prob.Plim(2)
        Ptype='max';
    end

    tt=[];
    zz=[];
    complete=0;

    while ~complete
        [tti,zzi,te,ze,ie]=ode78(@(t,z) TO_2BP_SEP(t,z,Ptype),tspan,z0,options);
    
        if ~isempty(ie)

            tt=[tt; tti(1:end)];   % to try correcting for duplicates
            zz=[zz; zzi(1:end,:)];

            ti=te(1);
            z0=ze(1,:).';
            tspan(1)=ti;

            switch ie(1)

                case 1*~strcmp(Ptype,'max')  % Pmed -> Pmax

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=TO_2BP_SEP(te,ze.','med');
                    dze_p=TO_2BP_SEP(te,ze.','max');

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve);

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='max';

                case 2*~strcmp(Ptype,'med')  % Pmax -> Pmed

                    re=ze(1:3).';
                    ve=ze(4:6).';

                    dze_m=TO_2BP_SEP(te,ze.','max');
                    dze_p=TO_2BP_SEP(te,ze.','med');

                    Psi=eye(14)+(dze_p(1:14)-dze_m(1:14))*[re.', zeros(1,11)]./(re.'*ve);

                    Phie_m=reshape(ze(15:210),[14,14]);
                    Phie_p=Psi*Phie_m;
                    z0(15:210)=reshape(Phie_p,[1,14*14]);

                    Ptype='med';

                case 3  % u_on -> u_off
                case 4  % u_off -> u_on
                otherwise
            end
        else
            tt=[tt; tti];
            zz=[zz; zzi];
            complete=1;
        end
    end
    
end

function [value,isterminal,direction]=crossings(t,y,prob)

    r=norm(y(1:3));
    [~,~,Sp]=MARGO_param(r);
    f1=Sp-prob.Plim(2);
    f2=SwFun(t,y,prob.epsilon);

    value=[f1;  % Power switching function                      
           f1;
           f2;  % throttle switching function (for TO)
           f2];

    isterminal=ones(4,1);
    
    direction=[+1;  % from Pmed to Pmax [ie=1]
               -1;  % from Pmax to Pmed [ie=2]
               +1;  % on to off         [ie=3]
               -1]; % off to on         [ie=4]

end