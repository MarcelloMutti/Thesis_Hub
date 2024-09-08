function [H] = Hamil(tt,yy,prob)
% adimesional inputs

%     u=1;
    H=zeros(length(tt),1);

    for i=1:length(H)

        r=norm(yy(i,1:3));
        [~,~,Sp]=MARGO_param(r);
    
        if Sp<prob.Plim(2)
            Ptype='med';
        elseif Sp>=prob.Plim(2)
            Ptype='max';
        end

        ff=TwBP_EL(tt(i),yy(i,:).',Ptype);
        ffx=ff(1:7);
        ll=yy(i,8:14).';
        H(i)=dot(ll,ffx);
    end
    
end

