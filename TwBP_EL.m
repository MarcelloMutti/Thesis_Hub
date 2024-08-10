function [dydPhi] = TwBP_EL(~, y, u, sc_param)
% v1 (master)

% adimensional dynamics

    rr=y(1:3);
    vv=y(4:6);
    m=y(7);

    ll=y(8:14);
    llr=y(8:10);
    llv=y(11:13);

    r=norm(rr);
    lv=norm(llv);

    vPhi=y(15:210);
    Phi=reshape(vPhi,[14,14]);

    gr=-rr/r^3;
    hv=zeros(3,1);

    G=3*(rr*rr.')/r^5-eye(3)/r^3;
    H=zeros(3,3);

    T=sc_param(1);
    c=sc_param(2);

    x=rr(1);
    y=rr(2);
    z=rr(3);

    lvx=llv(1);
    lvy=llv(2);
    lvz=llv(3);

    DyF=zeros(14,14);

    DyF(1:3,4:6)=eye(3);
    DyF(4:6,1:3)=G;
    DyF(4:6,7)=llv/lv*u*T/m^2;
    DyF(4:6,11:13)=-u*T/m*(eye(3)./lv-(llv*llv.')./lv^3);

%     DyF(8:10,1:3)=[[-(3*(2*lvx*x*abs(x)^2 - 5*lvx*x^2*real(x) + 2*lvx*x*abs(y)^2 + lvy*y*abs(x)^2 + 2*lvx*x*abs(z)^2 + lvy*y*abs(y)^2 + lvz*z*abs(x)^2 + lvy*y*abs(z)^2 + lvz*z*abs(y)^2 + lvz*z*abs(z)^2 + lvx*real(x)*abs(x)^2 + lvx*real(x)*abs(y)^2 + lvx*real(x)*abs(z)^2 - 5*lvy*x*y*real(x) - 5*lvz*x*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             -(3*(lvy*x*abs(x)^2 - 5*lvx*x^2*real(y) + lvy*x*abs(y)^2 + lvy*x*abs(z)^2 + lvx*real(y)*abs(x)^2 + lvx*real(y)*abs(y)^2 + lvx*real(y)*abs(z)^2 - 5*lvy*x*y*real(y) - 5*lvz*x*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             -(3*(lvz*x*abs(x)^2 - 5*lvx*x^2*real(z) + lvz*x*abs(y)^2 + lvz*x*abs(z)^2 + lvx*real(z)*abs(x)^2 + lvx*real(z)*abs(y)^2 + lvx*real(z)*abs(z)^2 - 5*lvy*x*y*real(z) - 5*lvz*x*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]
%                    [                                                                                                            -(3*(lvx*y*abs(x)^2 - 5*lvy*y^2*real(x) + lvx*y*abs(y)^2 + lvx*y*abs(z)^2 + lvy*real(x)*abs(x)^2 + lvy*real(x)*abs(y)^2 + lvy*real(x)*abs(z)^2 - 5*lvx*x*y*real(x) - 5*lvz*y*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), -(3*(lvx*x*abs(x)^2 - 5*lvy*y^2*real(y) + lvx*x*abs(y)^2 + 2*lvy*y*abs(x)^2 + lvx*x*abs(z)^2 + 2*lvy*y*abs(y)^2 + lvz*z*abs(x)^2 + 2*lvy*y*abs(z)^2 + lvz*z*abs(y)^2 + lvz*z*abs(z)^2 + lvy*real(y)*abs(x)^2 + lvy*real(y)*abs(y)^2 + lvy*real(y)*abs(z)^2 - 5*lvx*x*y*real(y) - 5*lvz*y*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             -(3*(lvz*y*abs(x)^2 - 5*lvy*y^2*real(z) + lvz*y*abs(y)^2 + lvz*y*abs(z)^2 + lvy*real(z)*abs(x)^2 + lvy*real(z)*abs(y)^2 + lvy*real(z)*abs(z)^2 - 5*lvx*x*y*real(z) - 5*lvz*y*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]
%                    [                                                                                                            -(3*(lvx*z*abs(x)^2 - 5*lvz*z^2*real(x) + lvx*z*abs(y)^2 + lvx*z*abs(z)^2 + lvz*real(x)*abs(x)^2 + lvz*real(x)*abs(y)^2 + lvz*real(x)*abs(z)^2 - 5*lvx*x*z*real(x) - 5*lvy*y*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             -(3*(lvy*z*abs(x)^2 - 5*lvz*z^2*real(y) + lvy*z*abs(y)^2 + lvy*z*abs(z)^2 + lvz*real(y)*abs(x)^2 + lvz*real(y)*abs(y)^2 + lvz*real(y)*abs(z)^2 - 5*lvx*x*z*real(y) - 5*lvy*y*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), -(3*(lvx*x*abs(x)^2 - 5*lvz*z^2*real(z) + lvx*x*abs(y)^2 + lvy*y*abs(x)^2 + lvx*x*abs(z)^2 + lvy*y*abs(y)^2 + 2*lvz*z*abs(x)^2 + lvy*y*abs(z)^2 + 2*lvz*z*abs(y)^2 + 2*lvz*z*abs(z)^2 + lvz*real(z)*abs(x)^2 + lvz*real(z)*abs(y)^2 + lvz*real(z)*abs(z)^2 - 5*lvx*x*z*real(z) - 5*lvy*y*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]];
    
    DyF(8:10,1:3)=-[[(3*(2*lvx*x*abs(x)^2 - 5*lvx*x^2*real(x) + 2*lvx*x*abs(y)^2 + 2*lvx*x*abs(z)^2 + lvy*y*abs(x)^2 + lvy*y*abs(y)^2 + lvy*y*abs(z)^2 + lvz*z*abs(x)^2 + lvz*z*abs(y)^2 + lvz*z*abs(z)^2 + lvx*real(x)*abs(x)^2 + lvx*real(x)*abs(y)^2 + lvx*real(x)*abs(z)^2 - 5*lvy*x*y*real(x) - 5*lvz*x*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             (3*(lvy*x*abs(x)^2 - 5*lvx*x^2*real(y) + lvy*x*abs(y)^2 + lvy*x*abs(z)^2 + lvx*real(y)*abs(x)^2 + lvx*real(y)*abs(y)^2 + lvx*real(y)*abs(z)^2 - 5*lvy*x*y*real(y) - 5*lvz*x*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             (3*(lvz*x*abs(x)^2 - 5*lvx*x^2*real(z) + lvz*x*abs(y)^2 + lvz*x*abs(z)^2 + lvx*real(z)*abs(x)^2 + lvx*real(z)*abs(y)^2 + lvx*real(z)*abs(z)^2 - 5*lvy*x*y*real(z) - 5*lvz*x*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]
                    [                                                                                                            (3*(lvx*y*abs(x)^2 - 5*lvy*y^2*real(x) + lvx*y*abs(y)^2 + lvx*y*abs(z)^2 + lvy*real(x)*abs(x)^2 + lvy*real(x)*abs(y)^2 + lvy*real(x)*abs(z)^2 - 5*lvx*x*y*real(x) - 5*lvz*y*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), (3*(lvx*x*abs(x)^2 - 5*lvy*y^2*real(y) + lvx*x*abs(y)^2 + lvx*x*abs(z)^2 + 2*lvy*y*abs(x)^2 + 2*lvy*y*abs(y)^2 + 2*lvy*y*abs(z)^2 + lvz*z*abs(x)^2 + lvz*z*abs(y)^2 + lvz*z*abs(z)^2 + lvy*real(y)*abs(x)^2 + lvy*real(y)*abs(y)^2 + lvy*real(y)*abs(z)^2 - 5*lvx*x*y*real(y) - 5*lvz*y*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             (3*(lvz*y*abs(x)^2 - 5*lvy*y^2*real(z) + lvz*y*abs(y)^2 + lvz*y*abs(z)^2 + lvy*real(z)*abs(x)^2 + lvy*real(z)*abs(y)^2 + lvy*real(z)*abs(z)^2 - 5*lvx*x*y*real(z) - 5*lvz*y*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]
                    [                                                                                                            (3*(lvx*z*abs(x)^2 - 5*lvz*z^2*real(x) + lvx*z*abs(y)^2 + lvx*z*abs(z)^2 + lvz*real(x)*abs(x)^2 + lvz*real(x)*abs(y)^2 + lvz*real(x)*abs(z)^2 - 5*lvx*x*z*real(x) - 5*lvy*y*z*real(x)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2),                                                                                                             (3*(lvy*z*abs(x)^2 - 5*lvz*z^2*real(y) + lvy*z*abs(y)^2 + lvy*z*abs(z)^2 + lvz*real(y)*abs(x)^2 + lvz*real(y)*abs(y)^2 + lvz*real(y)*abs(z)^2 - 5*lvx*x*z*real(y) - 5*lvy*y*z*real(y)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2), (3*(lvx*x*abs(x)^2 - 5*lvz*z^2*real(z) + lvx*x*abs(y)^2 + lvx*x*abs(z)^2 + lvy*y*abs(x)^2 + lvy*y*abs(y)^2 + lvy*y*abs(z)^2 + 2*lvz*z*abs(x)^2 + 2*lvz*z*abs(y)^2 + 2*lvz*z*abs(z)^2 + lvz*real(z)*abs(x)^2 + lvz*real(z)*abs(y)^2 + lvz*real(z)*abs(z)^2 - 5*lvx*x*z*real(z) - 5*lvy*y*z*real(z)))/(abs(x)^2 + abs(y)^2 + abs(z)^2)^(7/2)]];
    
    DyF(8:10,11:13)=-G.';
    DyF(11:13,8:10)=-eye(3);
    DyF(14,7)=2*lv*u*T/m^3;
    DyF(14,11:13)=-llv.'*u*T/(lv*m^2);

    dPhi=DyF*Phi;

    vdPhi=reshape(dPhi,[14*14,1]);

    dydPhi=zeros(210,1);

    dydPhi(1:3)=vv;
    dydPhi(4:6)=gr+hv-u*T/m*llv/lv;
    dydPhi(7)=-u*T/c;

    dydPhi(8:10)=-G.'*llv;
    dydPhi(11:13)=-H.'*llv-llr;
    dydPhi(14)=-u*lv*T/m^2;
    
    dydPhi(15:210)=vdPhi;

end

