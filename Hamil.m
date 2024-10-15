function [H] = Hamil(tt,zz,prob)
% adimesional inputs

%     u=1;
    H=zeros(length(tt),1);

    for i=1:length(H)

        r=norm(zz(i,1:3));
        [~,~,Sp]=MARGO_param(r);
    
        if Sp<prob.Plim(2)
            Ptype='med';
        elseif Sp>=prob.Plim(2)
            Ptype='max';
        end

        Sw=SwFun(tt(i),zz(i,:),prob.isFO);

        if prob.isFO==1

            ff=FO_2BP_SEP(tt(i),zz(i,:).',Ptype,prob.epsilon);
            u=1/2*(1-tanh(Sw/prob.epsilon));

        else

            ff=TO_2BP_SEP(tt(i),zz(i,:).',Ptype);
            u=1;

        end

        Tc=MARGO_param(r);
        T=Tc(1);
        c=Tc(2);

        ffx=ff(1:7);
        ll=zz(i,8:14).';
        H(i)=dot(ll,ffx)+prob.isFO*(T/c*(u-u*prob.epsilon*(1-u)));
    end
    
end

