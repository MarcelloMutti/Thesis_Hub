function TO_SEP_sym_dyn(prob)

    y=sym('y',[14,1]);
    assume(y,'real');
    
    rr=y(1:3);
    vv=y(4:6);
    m=y(7);
    llr=y(8:10);
    llv=y(11:13);
    lm=y(14);
    
    U=sym('U',[1,4]); % to be substituted
    assume(U,'real');
    
    LU=U(1);
    TU=U(2);
    MU=U(3);
    cf=U(4);% =1e-6
    
    syms g0; % to be substituted
    assume(g0,'real');
    
    MP=sym('MP',[3,5]);
    assume(MP,'real');
    
    ap=MP(1,:);
    bp=MP(2,:);
    cp=MP(3,:);
    
    %-MED-Dyn-computed-inside-F-and-DyF------------------------------------
    
    P=dot(norm(rr).^(0:4),cp);  % [W]
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    c=g0*I*TU/LU;         % [-]
    
    g=-rr/norm(rr)^3;
    h=zeros(3,1);
    
    % insert S and u dependency in EO/FO
    u=1;
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx);
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_med=[ffx; ffl];
    A_med=jacobian(FF_med,y);

    %-MAX-Dyn-computed-inside-F-and-DyF------------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    c=g0*I*TU/LU;         % [-]
    
    g=-rr/norm(rr)^3;
    h=zeros(3,1);
    
    % insert S and u dependency in EO/FO
    u=1;
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx);
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_max=[ffx; ffl];
    A_max=jacobian(FF_max,y);

    dir_name='TO_dyn';

    if ~exist(dir_name,'dir')
        mkdir(dir_name)
        addpath(dir_name)
    end

    matlabFunction(FF_med,A_med,'vars',{y,MP,U,g0},'file',fullfile(dir_name,'TO_med_dz'));
    matlabFunction(FF_max,A_max,'vars',{y,MP,U,g0},'file',fullfile(dir_name,'TO_max_dz'));

end