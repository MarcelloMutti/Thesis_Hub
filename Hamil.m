function [H] = Hamil(tt,zz,u,prob)
% adimesional inputs

    H=zeros(length(tt),1);

    for i=1:length(H)

        r=norm(zz(i,1:3));
        [~,~,Sp]=MARGO_param(r);
    
        if Sp<prob.Plim(2)
            Ptype='med';
        elseif Sp>=prob.Plim(2)
            Ptype='max';
        end

        ep=prob.epsilon;
        
        % Se=SwFun(tt(i),zz(i,:),prob.isFO);
        % 
        % if ep~=0
        %     if Se<-ep
        %         utype='on';
        %     elseif Se>ep
        %         utype='off';
        %     else
        %         utype='med';
        %     end
        % else
        %     if Se<0
        %         utype='on';
        %     else
        %         utype='off';
        %     end
        % end

        if u(i)==1
            utype='on';
        elseif u(i)==0
            utype='off';
        else
            utype='med';
        end

        if prob.isFO==1

            ff=FO_2BP_SEP(tt(i),zz(i,:).',Ptype,utype,ep);

        else

            ff=TO_2BP_SEP(tt(i),zz(i,:).',Ptype);

        end

        % if ep>0
        % 
        %     u=1.*(Se<-ep)+(ep-Se)./(2*ep).*(abs(Se)<=ep)+0;
        % 
        % else
        % 
        %     u=1.*(Se<0)+0;
        % 
        % end

        Tc=MARGO_param(r);
        T=Tc(1);
        c=Tc(2);

        ffx=ff(1:7);
        ll=zz(i,8:14).';
        H(i)=dot(ll,ffx)+prob.isFO*(T/c*(u(i)-u(i)*ep*(1-u(i))));
    end
    
end

