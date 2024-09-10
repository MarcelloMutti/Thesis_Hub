function [S] = SwFun(tt,zz,epsilon)
% adimesional inputs

%     c=sc_param(2);
    
    S=zeros(length(tt),1);
    
    for i=1:length(tt)
        
        if length(tt)==1
            r=norm(zz(1:3));
            m=zz(7);
            llv=zz(11:13);
            lm=zz(14);
        else
            r=norm(zz(i,1:3));
            m=zz(i,7);
            llv=zz(i,11:13);
            lm=zz(i,14);
        end

        Tc=MARGO_param(r);
    
        c=Tc(2);
    
        lv=norm(llv);
        
        S(i)=-m.*lv./c-lm+epsilon;
    
    end

end

