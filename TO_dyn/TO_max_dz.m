function [FF,A] = TO_max_dz(in1,in2,in3,g0)
%TO_max_dz
%    [FF,A] = TO_max_dz(IN1,IN2,IN3,G0)

%    This function was generated by the Symbolic Math Toolbox version 24.2.
%    02-Dec-2024 17:07:42

MP1_1 = in2(1);
MP1_2 = in2(4);
MP1_3 = in2(7);
MP1_4 = in2(10);
MP1_5 = in2(13);
MP2_1 = in2(2);
MP2_2 = in2(5);
MP2_3 = in2(8);
MP2_4 = in2(11);
MP2_5 = in2(14);
U1 = in3(:,1);
U2 = in3(:,2);
U3 = in3(:,3);
U4 = in3(:,4);
y1 = in1(1,:);
y2 = in1(2,:);
y3 = in1(3,:);
y4 = in1(4,:);
y5 = in1(5,:);
y6 = in1(6,:);
y7 = in1(7,:);
y8 = in1(8,:);
y9 = in1(9,:);
y10 = in1(10,:);
y11 = in1(11,:);
y12 = in1(12,:);
y13 = in1(13,:);
t2 = abs(y1);
t3 = abs(y2);
t4 = abs(y3);
t5 = abs(y11);
t6 = abs(y12);
t7 = abs(y13);
t8 = sign(y1);
t9 = sign(y2);
t10 = sign(y3);
t11 = sign(y11);
t12 = sign(y12);
t13 = sign(y13);
t14 = U2.^2;
t21 = y1.*y11.*3.0;
t22 = y2.*y12.*3.0;
t23 = y3.*y13.*3.0;
t24 = 1.0./U1;
t25 = 1.0./U3;
t26 = 1.0./y7;
t28 = MP1_2.*1.2e+2;
t29 = MP1_3.*1.44e+4;
t30 = MP1_4.*1.728e+6;
t31 = MP1_5.*2.0736e+8;
t15 = t2.^2;
t16 = t3.^2;
t17 = t4.^2;
t18 = t5.^2;
t19 = t6.^2;
t20 = t7.^2;
t27 = t26.^2;
t34 = t21+t22+t23;
t50 = MP1_1+t28+t29+t30+t31;
t32 = t15+t16+t17;
t33 = t18+t19+t20;
t35 = sqrt(t33);
t36 = 1.0./t32.^(3.0./2.0);
t37 = 1.0./t32.^(5.0./2.0);
t38 = 1.0./t32.^(7.0./2.0);
t39 = 1.0./t35;
t41 = -t36;
t42 = t37.*y1.*y2.*3.0;
t43 = t37.*y1.*y3.*3.0;
t44 = t37.*y2.*y3.*3.0;
t48 = t34.*t37;
t40 = t39.^3;
t45 = -t42;
t46 = -t43;
t47 = -t44;
t49 = -t48;
t51 = U4.*t14.*t24.*t25.*t26.*t39.*t50;
t52 = -t51;
FF = [y4;y5;y6;t41.*y1+t52.*y11;t41.*y2+t52.*y12;t41.*y3+t52.*y13;-(U2.*U4.*t25.*t50)./(g0.*(MP2_1+MP2_2.*1.2e+2+MP2_3.*1.44e+4+MP2_4.*1.728e+6+MP2_5.*2.0736e+8));t36.*y11+t49.*y1;t36.*y12+t49.*y2;t36.*y13+t49.*y3;-y8;-y9;-y10;-U4.*t14.*t24.*t25.*t27.*t35.*t50];
if nargout > 1
    mt1 = [0.0,0.0,0.0,t41+t2.*t8.*t37.*y1.*3.0,t2.*t8.*t37.*y2.*3.0,t2.*t8.*t37.*y3.*3.0,0.0,t49-t37.*y1.*y11.*3.0-t2.*t8.*t37.*y11.*3.0+t2.*t8.*t34.*t38.*y1.*5.0,t37.*y2.*y11.*-3.0-t2.*t8.*t37.*y12.*3.0+t2.*t8.*t34.*t38.*y2.*5.0,t37.*y3.*y11.*-3.0-t2.*t8.*t37.*y13.*3.0+t2.*t8.*t34.*t38.*y3.*5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t9.*t37.*y1.*3.0,t41+t3.*t9.*t37.*y2.*3.0,t3.*t9.*t37.*y3.*3.0,0.0,t37.*y1.*y12.*-3.0-t3.*t9.*t37.*y11.*3.0+t3.*t9.*t34.*t38.*y1.*5.0,t49-t37.*y2.*y12.*3.0-t3.*t9.*t37.*y12.*3.0+t3.*t9.*t34.*t38.*y2.*5.0,t37.*y3.*y12.*-3.0-t3.*t9.*t37.*y13.*3.0+t3.*t9.*t34.*t38.*y3.*5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t10.*t37.*y1.*3.0,t4.*t10.*t37.*y2.*3.0];
    mt2 = [t41+t4.*t10.*t37.*y3.*3.0,0.0,t37.*y1.*y13.*-3.0-t4.*t10.*t37.*y11.*3.0+t4.*t10.*t34.*t38.*y1.*5.0,t37.*y2.*y13.*-3.0-t4.*t10.*t37.*y12.*3.0+t4.*t10.*t34.*t38.*y2.*5.0,t49-t37.*y3.*y13.*3.0-t4.*t10.*t37.*y13.*3.0+t4.*t10.*t34.*t38.*y3.*5.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t14.*t24.*t25.*t27.*t39.*t50.*y11,U4.*t14.*t24.*t25.*t27.*t39.*t50.*y12,U4.*t14.*t24.*t25.*t27.*t39.*t50.*y13,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t14.*t24.*t25.*t26.^3.*t35.*t50.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0];
    mt3 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t52+U4.*t5.*t11.*t14.*t24.*t25.*t26.*t40.*t50.*y11,U4.*t5.*t11.*t14.*t24.*t25.*t26.*t40.*t50.*y12,U4.*t5.*t11.*t14.*t24.*t25.*t26.*t40.*t50.*y13,0.0,t36-t37.*y1.^2.*3.0,t45,t46,0.0,0.0,0.0,-U4.*t5.*t11.*t14.*t24.*t25.*t27.*t39.*t50,0.0,0.0,0.0,U4.*t6.*t12.*t14.*t24.*t25.*t26.*t40.*t50.*y11,t52+U4.*t6.*t12.*t14.*t24.*t25.*t26.*t40.*t50.*y12,U4.*t6.*t12.*t14.*t24.*t25.*t26.*t40.*t50.*y13,0.0,t45,t36-t37.*y2.^2.*3.0,t47,0.0,0.0,0.0,-U4.*t6.*t12.*t14.*t24.*t25.*t27.*t39.*t50,0.0,0.0,0.0,U4.*t7.*t13.*t14.*t24.*t25.*t26.*t40.*t50.*y11,U4.*t7.*t13.*t14.*t24.*t25.*t26.*t40.*t50.*y12,t52+U4.*t7.*t13.*t14.*t24.*t25.*t26.*t40.*t50.*y13,0.0,t46,t47,t36-t37.*y3.^2.*3.0,0.0,0.0,0.0];
    mt4 = [-U4.*t7.*t13.*t14.*t24.*t25.*t27.*t39.*t50,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    A = reshape([mt1,mt2,mt3,mt4],14,14);
end
end
