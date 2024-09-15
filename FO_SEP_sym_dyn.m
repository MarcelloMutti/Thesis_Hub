function FO_SEP_sym_dyn(prob)

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

    syms epsilon
    assume(epsilon,'real');
    
    %-MED-Dyn-computed-inside-F-and-DyF------------------------------------
    
    P=dot(norm(rr).^(0:4),cp);  % [W]
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    c=g0*I*TU/LU;         % [-]
    
    g=-rr/norm(rr)^3;
    h=zeros(3,1);

    %-thrust-on------------------------------------------------------------
    u=1;
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_med_uon=[ffx; ffl];
    A_med_uon=jacobian(FF_med_uon,y);

    %-thrust-off-----------------------------------------------------------
    u=0;

    ffx=[vv;
     g+h-u*llv/norm(llv)*T/m;
     -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_med_uoff=[ffx; ffl];
    A_med_uoff=jacobian(FF_med_uoff,y);

    %-thrust-med-----------------------------------------------------------
    Se=-norm(llv)*c/m-lm+1;
    u=(epsilon-Se)/(2*epsilon);
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_med_umed=[ffx; ffl];
    A_med_umed=jacobian(FF_med_umed,y);



    %-MAX-Dyn-computed-inside-F-and-DyF------------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    c=g0*I*TU/LU;         % [-]
    
    g=-rr/norm(rr)^3;
    h=zeros(3,1);

    %-thrust-on------------------------------------------------------------
    u=1;
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_max_uon=[ffx; ffl];
    A_max_uon=jacobian(FF_max_uon,y);

    %-thrust-off-----------------------------------------------------------
    u=0;

    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_max_uoff=[ffx; ffl];
    A_max_uoff=jacobian(FF_max_uoff,y);

    %-thrust-med-----------------------------------------------------------
    Se=-norm(llv)*c/m-lm+1;
    u=(epsilon-Se)/(2*epsilon);
    
    ffx=[vv;
         g+h-u*llv/norm(llv)*T/m;
         -u*T/c];
    
    Hamil=dot([llr; llv; lm],ffx)+T/c*(u-epsilon*u*(1-u));
    
    ffl=-jacobian(Hamil,[rr; vv; m]).';
    
    FF_max_umed=[ffx; ffl];
    A_max_umed=jacobian(FF_max_umed,y);
    

    dir_name='FO_dyn';

    if ~exist(dir_name,'dir')
        mkdir(dir_name)
        addpath(dir_name)
    end

    matlabFunction(FF_med_uon, A_med_uon, 'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_uon_dz' ));
    matlabFunction(FF_max_uon, A_max_uon, 'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uon_dz' ));
    matlabFunction(FF_med_uoff,A_med_uoff,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_uoff_dz'));
    matlabFunction(FF_max_uoff,A_max_uoff,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uoff_dz'));
    matlabFunction(FF_med_umed,A_med_umed,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_umed_dz'));
    matlabFunction(FF_max_umed,A_max_umed,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_umed_dz'));

end