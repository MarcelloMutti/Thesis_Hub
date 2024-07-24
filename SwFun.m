function [S] = SwFun(~,yy,sc_param,epsilon)
% adimesional inputs

    c=sc_param(2);

    m=yy(:,7);
    llv=yy(:,11:13);
    lv=sqrt(sum(llv.^2,2));
    lm=yy(:,14);
    
    S=-m.*lv./c-lm+epsilon;
    
end

