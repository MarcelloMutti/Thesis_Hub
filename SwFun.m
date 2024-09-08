function [S] = SwFun(tt,yy,epsilon)
% adimesional inputs

%     c=sc_param(2);
    
    S=zeros(length(tt),1);
    
    for i=1:length(tt)
        
        if length(tt)==1
            r=norm(yy(1:3));
            m=yy(7);
            llv=yy(11:13);
            lm=yy(14);
        else
            r=norm(yy(i,1:3));
            m=yy(i,7);
            llv=yy(i,11:13);
            lm=yy(i,14);
        end

        Tc=MARGO_param(r);
    
        c=Tc(2);
    
%         m=yy(i,7);
%         llv=yy(i,11:13);
    
        lv=norm(llv);
%     
%         lm=yy(i,14);
        
        S(i)=-m.*lv./c-lm+epsilon;
    
    end

%     r=sqrt(sum(prob(end).yy(:,1:3).^2,2));
% 
%     Tc=MARGO_param(r);
%     c=Tc(2);
%      
%     m=yy(:,7);

end

