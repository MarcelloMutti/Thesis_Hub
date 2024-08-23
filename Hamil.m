function [H] = Hamil(tt,yy)
% adimesional inputs

    u=1;
    H=zeros(length(tt),1);

    for i=1:length(H)
        ff=TwBP_EL(tt(i),yy(i,:).',u);
        ffx=ff(1:7);
        ll=yy(i,8:14).';
        H(i)=dot(ll,ffx);
    end
    
end

