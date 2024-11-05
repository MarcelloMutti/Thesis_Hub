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

    syms u
    assume('real' & u<1 & u>0)

    syms epsilon
    assume(epsilon,'real');
    
    %-MED-Dyn-computed-inside-F-and-DyF------------------------------------
    
    P=dot(norm(rr).^(0:4),cp);  % [W]
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]
    
    dP=dot(norm(rr).^(0:3),(1:4).*cp(2:end));
    dT=dot(P.^(0:3),(1:4).*ap(2:end));
    dI=dot(P.^(0:3),(1:4).*bp(2:end));    
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]

    tr=cf*TU^2/(MU*LU)*dT*dP*rr/norm(rr);
    ir=1/TU*dI*dP*rr/norm(rr);
    c=I*TU^2/LU*g0;

    %-med-on---------------------------------------------------------------
    u=1;
        
    FF_med_uon=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3+u*norm(llv)/m*tr+((lm-1+epsilon)*u-epsilon*u^2)/(c)*(tr-T/I*ir);
        -llr;
        -u*norm(llv)*T/m^2];

    A_med_uon=jacobian(FF_med_uon,y);

    %-med-off--------------------------------------------------------------
    u=0;
        
    FF_med_uoff=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3+u*norm(llv)/m*tr+((lm-1+epsilon)*u-epsilon*u^2)/(c)*(tr-T/I*ir);
        -llr;
        -u*norm(llv)*T/m^2];

    A_med_uoff=jacobian(FF_med_uoff,y);

    %-med-med--------------------------------------------------------------
    S=-norm(llv)*c/m-lm+1;
    u=(epsilon-S)/(2*epsilon);
        
    FF_med_umed=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3+u*norm(llv)/m*tr+((lm-1+epsilon)*u-epsilon*u^2)/(c)*(tr-T/I*ir);
        -llr;
        -u*norm(llv)*T/m^2];

    A_med_umed=jacobian(FF_med_umed,y);

    %-MAX-Dyn-computed-inside-F-and-DyF------------------------------------

    P=sym(prob.Plim(2));
    T=dot(P.^(0:4),ap);         % [mN]
    I=dot(P.^(0:4),bp);         % [s]   
    
    T=cf*T*TU^2/(MU*LU);  % [-]
    I=I/TU;               % [-]
    c=I*TU^2/LU*g0;

    %-max-on---------------------------------------------------------------
    u=1;
        
    FF_max_uon=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A_max_uon=jacobian(FF_max_uon,y);

    %-max-off--------------------------------------------------------------
    u=0;
        
    FF_max_uoff=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A_max_uoff=jacobian(FF_max_uoff,y);

    %-med-med--------------------------------------------------------------
    S=-norm(llv)*c/m-lm+1;
    u=(epsilon-S)/(2*epsilon);
        
    FF_max_umed=[vv;
        -rr/norm(rr)^3-u*T/m*llv/norm(llv);
        -u*T/(c);
        -3*dot(rr,llv)*rr/norm(rr)^5+llv/norm(rr)^3;
        -llr;
        -u*norm(llv)*T/m^2];

    A_max_umed=jacobian(FF_max_umed,y);

    
    dir_name='FO_dyn';

    if ~exist(dir_name,'dir')
        mkdir(dir_name)
        addpath(dir_name)
    end

%     FF_meds=subs(FF_med,u,uS);
%     FF_maxs=subs(FF_max,u,uS);
%     A_meds =subs(A_med ,u,uS);
%     A_maxs =subs(A_max ,u,uS);

    matlabFunction(FF_med_uon, A_med_uon ,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_uon_dz' ));
    matlabFunction(FF_med_uoff,A_med_uoff,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_uoff_dz'));
    matlabFunction(FF_med_umed,A_med_umed,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_med_umed_dz'));
    matlabFunction(FF_max_uon, A_max_uon ,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uon_dz' ));
    matlabFunction(FF_max_uoff,A_max_uoff,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_uoff_dz'));
    matlabFunction(FF_max_umed,A_max_umed,'vars',{y,MP,U,g0,epsilon},'file',fullfile(dir_name,'FO_max_umed_dz'));
    
   
end