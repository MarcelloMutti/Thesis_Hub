function [S] = SwFun(tt,yy,epsilon)
% adimesional inputs

%     c=sc_param(2);
    
    S=zeros(length(tt),1);
    
    for i=1:length(tt)

    r=norm(yy(i,1:3));
    Tc=MARGO_param(r);

    c=Tc(2);

    m=yy(i,7);
    llv=yy(i,11:13);

    lv=norm(llv);

    lm=yy(i,14);
    
    S(i)=-m.*lv./c-lm+epsilon;
    
    end
end

