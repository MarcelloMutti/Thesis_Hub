function [H] = Hamil(tt,yy,sc_param)
% adimesional inputs

    u=1;
    H=zeros(size(yy,1),1);

    for i=1:length(H)
        ff=TwBP_EL(tt(i),yy(1,:).',u,sc_param);
        ffx=ff(1:7);
        ll=yy(i,8:14).';
        H(i)=dot(ffx,ll);
    end
end

