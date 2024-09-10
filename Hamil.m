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

        ff=TwBP_EL(tt(i),zz(i,:).',Ptype);
        ffx=ff(1:7);
        ll=zz(i,8:14).';
        H(i)=dot(ll,ffx);
    end
    
end

