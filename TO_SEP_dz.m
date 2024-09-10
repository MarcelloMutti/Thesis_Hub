function [FF_med,A_med,FF_max,A_max] = TO_SEP_dz(in1,in2,in3,g0)
%TO_SEP_dz
%    [FF_med,A_med,FF_max,A_max] = TO_SEP_dz(IN1,IN2,IN3,G0)

%    This function was generated by the Symbolic Math Toolbox version 9.2.
%    10-Sep-2024 17:30:12

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
MP3_1 = in2(3);
MP3_2 = in2(6);
MP3_3 = in2(9);
MP3_4 = in2(12);
MP3_5 = in2(15);
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
y14 = in1(14,:);
t2 = abs(y1);
t3 = abs(y2);
t4 = abs(y3);
t5 = abs(y11);
t6 = abs(y12);
t7 = abs(y13);
t8 = dirac(y1);
t9 = dirac(y2);
t10 = dirac(y3);
t11 = sign(y1);
t12 = sign(y2);
t13 = sign(y3);
t14 = sign(y11);
t15 = sign(y12);
t16 = sign(y13);
t17 = U2.^2;
t18 = y11.^2;
t19 = y12.^2;
t20 = y13.^2;
t30 = 1.0./U1;
t31 = 1.0./U3;
t32 = 1.0./g0;
t33 = -y8;
t34 = -y9;
t35 = -y10;
t36 = 1.0./y7;
t39 = MP1_2.*1.2e+2;
t46 = MP1_3.*1.44e+4;
t50 = MP1_4.*1.728e+6;
t51 = MP1_5.*2.0736e+8;
t21 = t2.^2;
t22 = t3.^2;
t23 = t4.^2;
t24 = t5.^2;
t25 = t6.^2;
t26 = t7.^2;
t27 = t11.^2;
t28 = t12.^2;
t29 = t13.^2;
t37 = t36.^2;
t38 = t36.^3;
t40 = MP3_3.*t2.*t8.*4.0;
t41 = MP3_3.*t3.*t9.*4.0;
t42 = MP3_3.*t4.*t10.*4.0;
t43 = MP3_3.*t2.*t11.*2.0;
t44 = MP3_3.*t3.*t12.*2.0;
t45 = MP3_3.*t4.*t13.*2.0;
t55 = MP3_5.*t2.*t3.*t11.*t12.*8.0;
t56 = MP3_5.*t2.*t4.*t11.*t13.*8.0;
t57 = MP3_5.*t3.*t4.*t12.*t13.*8.0;
t176 = MP1_1+t39+t46+t50+t51;
t47 = MP3_3.*t27.*2.0;
t48 = MP3_3.*t28.*2.0;
t49 = MP3_3.*t29.*2.0;
t52 = MP3_5.*t21.*t27.*8.0;
t53 = MP3_5.*t22.*t28.*8.0;
t54 = MP3_5.*t23.*t29.*8.0;
t58 = t21+t22+t23;
t59 = t24+t25+t26;
t60 = t58.^2;
t61 = MP3_3.*t58;
t63 = sqrt(t58);
t69 = 1.0./sqrt(t59);
t73 = MP3_5.*t27.*t58.*4.0;
t74 = MP3_5.*t28.*t58.*4.0;
t75 = MP3_5.*t29.*t58.*4.0;
t79 = MP3_5.*t2.*t11.*t58.*4.0;
t80 = MP3_5.*t3.*t12.*t58.*4.0;
t81 = MP3_5.*t4.*t13.*t58.*4.0;
t82 = MP3_5.*t2.*t8.*t58.*8.0;
t83 = MP3_5.*t3.*t9.*t58.*8.0;
t84 = MP3_5.*t4.*t10.*t58.*8.0;
t62 = MP3_5.*t60;
t64 = t63.^3;
t65 = 1.0./t63;
t70 = t69.^3;
t71 = MP3_2.*t63;
t92 = MP3_4.*t27.*t63.*3.0;
t93 = MP3_4.*t28.*t63.*3.0;
t94 = MP3_4.*t29.*t63.*3.0;
t101 = MP3_4.*t2.*t8.*t63.*6.0;
t102 = MP3_4.*t3.*t9.*t63.*6.0;
t103 = MP3_4.*t4.*t10.*t63.*6.0;
t107 = MP3_4.*t2.*t11.*t63.*3.0;
t108 = MP3_4.*t3.*t12.*t63.*3.0;
t109 = MP3_4.*t4.*t13.*t63.*3.0;
t180 = U4.*t17.*t30.*t31.*t36.*t69.*t176;
t66 = 1.0./t64;
t67 = t65.^5;
t68 = t65.^7;
t72 = MP3_4.*t64;
t89 = MP3_2.*t27.*t65;
t90 = MP3_2.*t28.*t65;
t91 = MP3_2.*t29.*t65;
t104 = MP3_2.*t2.*t11.*t65;
t105 = MP3_2.*t3.*t12.*t65;
t106 = MP3_2.*t4.*t13.*t65;
t125 = MP3_2.*t2.*t8.*t65.*2.0;
t126 = MP3_2.*t3.*t9.*t65.*2.0;
t127 = MP3_2.*t4.*t10.*t65.*2.0;
t146 = MP3_4.*t21.*t27.*t65.*3.0;
t147 = MP3_4.*t22.*t28.*t65.*3.0;
t148 = MP3_4.*t23.*t29.*t65.*3.0;
t149 = MP3_4.*t2.*t3.*t11.*t12.*t65.*3.0;
t150 = MP3_4.*t2.*t4.*t11.*t13.*t65.*3.0;
t151 = MP3_4.*t3.*t4.*t12.*t13.*t65.*3.0;
t181 = -t180;
t76 = t66.*y1;
t77 = t66.*y2;
t78 = t66.*y3;
t85 = -t66;
t95 = t2.*t11.*t67.*3.0;
t96 = t3.*t12.*t67.*3.0;
t97 = t2.*t11.*t67.*6.0;
t98 = t4.*t13.*t67.*3.0;
t99 = t3.*t12.*t67.*6.0;
t100 = t4.*t13.*t67.*6.0;
t119 = t27.*t67.*y1.*3.0;
t120 = t28.*t67.*y2.*3.0;
t121 = t29.*t67.*y3.*3.0;
t128 = t2.*t8.*t67.*y1.*6.0;
t129 = t3.*t9.*t67.*y2.*6.0;
t130 = t4.*t10.*t67.*y3.*6.0;
t131 = t2.*t11.*t67.*y1.*-3.0;
t132 = t2.*t11.*t67.*y2.*-3.0;
t133 = t3.*t12.*t67.*y1.*-3.0;
t134 = t2.*t11.*t67.*y3.*-3.0;
t135 = t3.*t12.*t67.*y2.*-3.0;
t136 = t4.*t13.*t67.*y1.*-3.0;
t137 = t3.*t12.*t67.*y3.*-3.0;
t138 = t4.*t13.*t67.*y2.*-3.0;
t139 = t4.*t13.*t67.*y3.*-3.0;
t140 = MP3_2.*t21.*t27.*t66;
t141 = MP3_2.*t22.*t28.*t66;
t142 = MP3_2.*t23.*t29.*t66;
t143 = MP3_2.*t2.*t3.*t11.*t12.*t66;
t144 = MP3_2.*t2.*t4.*t11.*t13.*t66;
t145 = MP3_2.*t3.*t4.*t12.*t13.*t66;
t158 = t21.*t27.*t68.*y1.*1.5e+1;
t159 = t22.*t28.*t68.*y2.*1.5e+1;
t160 = t23.*t29.*t68.*y3.*1.5e+1;
t161 = t2.*t3.*t11.*t12.*t68.*y1.*1.5e+1;
t162 = t2.*t3.*t11.*t12.*t68.*y2.*1.5e+1;
t163 = t2.*t4.*t11.*t13.*t68.*y1.*1.5e+1;
t164 = t2.*t3.*t11.*t12.*t68.*y3.*1.5e+1;
t165 = t2.*t4.*t11.*t13.*t68.*y2.*1.5e+1;
t166 = t3.*t4.*t12.*t13.*t68.*y1.*1.5e+1;
t167 = t2.*t4.*t11.*t13.*t68.*y3.*1.5e+1;
t168 = t3.*t4.*t12.*t13.*t68.*y2.*1.5e+1;
t169 = t3.*t4.*t12.*t13.*t68.*y3.*1.5e+1;
t182 = MP3_1+t61+t62+t71+t72;
t197 = t43+t79+t104+t107;
t198 = t44+t80+t105+t108;
t199 = t45+t81+t106+t109;
t86 = -t76;
t87 = -t77;
t88 = -t78;
t110 = t95.*y1;
t111 = t95.*y2;
t112 = t96.*y1;
t113 = t95.*y3;
t114 = t96.*y2;
t115 = t98.*y1;
t116 = t96.*y3;
t117 = t98.*y2;
t118 = t98.*y3;
t122 = -t95;
t123 = -t96;
t124 = -t98;
t152 = MP3_2.*t21.*t27.*t85;
t153 = MP3_2.*t22.*t28.*t85;
t154 = MP3_2.*t23.*t29.*t85;
t155 = MP3_2.*t2.*t3.*t11.*t12.*t85;
t156 = MP3_2.*t2.*t4.*t11.*t13.*t85;
t157 = MP3_2.*t3.*t4.*t12.*t13.*t85;
t170 = t166.*y11;
t171 = t165.*y12;
t172 = t164.*y13;
t173 = -t158;
t174 = -t159;
t175 = -t160;
t177 = t66+t131;
t178 = t66+t135;
t179 = t66+t139;
t186 = t182.^2;
t187 = t182.^3;
t189 = MP1_2.*t182;
t190 = MP2_2.*t182;
t206 = t197.^2;
t207 = t198.^2;
t208 = t199.^2;
t209 = MP1_2.*t197;
t210 = MP1_2.*t198;
t211 = MP2_2.*t197;
t212 = MP1_2.*t199;
t213 = MP2_2.*t198;
t214 = MP2_2.*t199;
t236 = MP1_3.*t182.*t197.*2.0;
t237 = MP1_3.*t182.*t198.*2.0;
t238 = MP2_3.*t182.*t197.*2.0;
t239 = MP1_3.*t182.*t199.*2.0;
t240 = MP2_3.*t182.*t198.*2.0;
t241 = MP2_3.*t182.*t199.*2.0;
t257 = MP1_3.*t197.*t198.*2.0;
t258 = MP1_3.*t197.*t199.*2.0;
t259 = MP2_3.*t197.*t198.*2.0;
t260 = MP1_3.*t198.*t199.*2.0;
t261 = MP2_3.*t197.*t199.*2.0;
t262 = MP2_3.*t198.*t199.*2.0;
t272 = MP1_4.*t182.*t197.*t198.*6.0;
t273 = MP1_4.*t182.*t197.*t199.*6.0;
t274 = MP2_4.*t182.*t197.*t198.*6.0;
t275 = MP1_4.*t182.*t198.*t199.*6.0;
t276 = MP2_4.*t182.*t197.*t199.*6.0;
t277 = MP2_4.*t182.*t198.*t199.*6.0;
t183 = t55+t149+t155;
t184 = t56+t150+t156;
t185 = t57+t151+t157;
t188 = t186.^2;
t200 = MP1_3.*t186;
t201 = MP1_4.*t187;
t203 = MP2_3.*t186;
t204 = MP2_4.*t187;
t215 = MP1_3.*t206.*2.0;
t216 = MP1_3.*t207.*2.0;
t217 = MP1_3.*t208.*2.0;
t242 = MP1_4.*t186.*t197.*3.0;
t243 = MP1_4.*t182.*t206.*6.0;
t244 = MP1_5.*t187.*t197.*4.0;
t245 = MP1_4.*t186.*t198.*3.0;
t246 = MP2_4.*t186.*t197.*3.0;
t247 = MP1_4.*t182.*t207.*6.0;
t248 = MP1_5.*t187.*t198.*4.0;
t249 = MP2_5.*t187.*t197.*4.0;
t250 = MP1_4.*t186.*t199.*3.0;
t251 = MP2_4.*t186.*t198.*3.0;
t252 = MP1_4.*t182.*t208.*6.0;
t253 = MP1_5.*t187.*t199.*4.0;
t254 = MP2_5.*t187.*t198.*4.0;
t255 = MP2_4.*t186.*t199.*3.0;
t256 = MP2_5.*t187.*t199.*4.0;
t263 = MP1_5.*t186.*t206.*1.2e+1;
t264 = MP1_5.*t186.*t207.*1.2e+1;
t265 = MP1_5.*t186.*t208.*1.2e+1;
t266 = t40+t47+t52+t73+t82+t89+t92+t101+t125+t146+t152;
t267 = t41+t48+t53+t74+t83+t90+t93+t102+t126+t147+t153;
t268 = t42+t49+t54+t75+t84+t91+t94+t103+t127+t148+t154;
t278 = MP1_5.*t186.*t197.*t198.*1.2e+1;
t279 = MP1_5.*t186.*t197.*t199.*1.2e+1;
t280 = MP2_5.*t186.*t197.*t198.*1.2e+1;
t281 = MP1_5.*t186.*t198.*t199.*1.2e+1;
t282 = MP2_5.*t186.*t197.*t199.*1.2e+1;
t283 = MP2_5.*t186.*t198.*t199.*1.2e+1;
t191 = MP1_2.*t183;
t192 = MP1_2.*t184;
t193 = MP2_2.*t183;
t194 = MP1_2.*t185;
t195 = MP2_2.*t184;
t196 = MP2_2.*t185;
t202 = MP1_5.*t188;
t205 = MP2_5.*t188;
t218 = MP1_3.*t182.*t183.*2.0;
t219 = MP1_3.*t182.*t184.*2.0;
t220 = MP2_3.*t182.*t183.*2.0;
t221 = MP1_3.*t182.*t185.*2.0;
t222 = MP2_3.*t182.*t184.*2.0;
t223 = MP2_3.*t182.*t185.*2.0;
t224 = MP1_4.*t183.*t186.*3.0;
t225 = MP1_5.*t183.*t187.*4.0;
t226 = MP1_4.*t184.*t186.*3.0;
t227 = MP1_5.*t184.*t187.*4.0;
t228 = MP2_4.*t183.*t186.*3.0;
t229 = MP1_4.*t185.*t186.*3.0;
t230 = MP2_5.*t183.*t187.*4.0;
t231 = MP1_5.*t185.*t187.*4.0;
t232 = MP2_4.*t184.*t186.*3.0;
t233 = MP2_5.*t184.*t187.*4.0;
t234 = MP2_4.*t185.*t186.*3.0;
t235 = MP2_5.*t185.*t187.*4.0;
t269 = MP1_2.*t266;
t270 = MP1_2.*t267;
t271 = MP1_2.*t268;
t284 = MP1_3.*t182.*t266.*2.0;
t285 = MP1_3.*t182.*t267.*2.0;
t286 = MP1_3.*t182.*t268.*2.0;
t287 = MP1_4.*t186.*t266.*3.0;
t288 = MP1_5.*t187.*t266.*4.0;
t289 = MP1_4.*t186.*t267.*3.0;
t290 = MP1_5.*t187.*t267.*4.0;
t291 = MP1_4.*t186.*t268.*3.0;
t292 = MP1_5.*t187.*t268.*4.0;
t300 = t209+t236+t242+t244;
t301 = t210+t237+t245+t248;
t302 = t211+t238+t246+t249;
t303 = t212+t239+t250+t253;
t304 = t213+t240+t251+t254;
t305 = t214+t241+t255+t256;
t293 = MP1_1+t189+t200+t201+t202;
t294 = MP2_1+t190+t203+t204+t205;
t306 = U4.*t17.*t30.*t31.*t36.*t69.*t300;
t307 = U4.*t17.*t30.*t31.*t36.*t69.*t301;
t308 = U4.*t17.*t30.*t31.*t36.*t69.*t303;
t318 = U4.*t17.*t18.*t30.*t31.*t37.*t69.*t300;
t319 = U4.*t17.*t19.*t30.*t31.*t37.*t69.*t300;
t320 = U4.*t17.*t20.*t30.*t31.*t37.*t69.*t300;
t321 = U4.*t17.*t18.*t30.*t31.*t37.*t69.*t301;
t322 = U4.*t17.*t19.*t30.*t31.*t37.*t69.*t301;
t323 = U4.*t17.*t20.*t30.*t31.*t37.*t69.*t301;
t324 = U4.*t17.*t18.*t30.*t31.*t37.*t69.*t303;
t325 = U4.*t17.*t19.*t30.*t31.*t37.*t69.*t303;
t326 = U4.*t17.*t20.*t30.*t31.*t37.*t69.*t303;
t351 = t191+t218+t224+t225+t257+t272+t278;
t352 = t192+t219+t226+t227+t258+t273+t279;
t353 = t193+t220+t228+t230+t259+t274+t280;
t354 = t194+t221+t229+t231+t260+t275+t281;
t355 = t195+t222+t232+t233+t261+t276+t282;
t356 = t196+t223+t234+t235+t262+t277+t283;
t402 = t215+t243+t263+t269+t284+t287+t288;
t403 = t216+t247+t264+t270+t285+t289+t290;
t404 = t217+t252+t265+t271+t286+t291+t292;
t295 = 1.0./t294;
t298 = U4.*t17.*t30.*t31.*t36.*t69.*t293;
t309 = t306.*y11;
t310 = t306.*y12;
t311 = t306.*y13;
t312 = t307.*y11;
t313 = t307.*y12;
t314 = t307.*y13;
t315 = t308.*y11;
t316 = t308.*y12;
t317 = t308.*y13;
t333 = -t318;
t334 = -t319;
t335 = -t320;
t336 = -t321;
t337 = -t322;
t338 = -t323;
t339 = -t324;
t340 = -t325;
t341 = -t326;
t360 = U4.*t17.*t30.*t31.*t36.*t69.*t351.*y11;
t361 = U4.*t17.*t30.*t31.*t36.*t69.*t351.*y12;
t362 = U4.*t17.*t30.*t31.*t36.*t69.*t351.*y13;
t363 = U4.*t17.*t30.*t31.*t36.*t69.*t352.*y11;
t364 = U4.*t17.*t30.*t31.*t36.*t69.*t352.*y12;
t365 = U4.*t17.*t30.*t31.*t36.*t69.*t352.*y13;
t366 = U4.*t17.*t30.*t31.*t36.*t69.*t354.*y11;
t367 = U4.*t17.*t30.*t31.*t36.*t69.*t354.*y12;
t368 = U4.*t17.*t30.*t31.*t36.*t69.*t354.*y13;
t296 = t295.^2;
t297 = t295.^3;
t299 = -t298;
t327 = -t310;
t328 = -t311;
t329 = -t312;
t330 = -t314;
t331 = -t315;
t332 = -t316;
t348 = U2.*U4.*t31.*t32.*t295.*t300;
t349 = U2.*U4.*t31.*t32.*t295.*t301;
t350 = U2.*U4.*t31.*t32.*t295.*t303;
t369 = t164+t362;
t370 = t165+t364;
t371 = t166+t366;
t375 = t122+t162+t361;
t376 = t123+t161+t360;
t377 = t122+t167+t365;
t378 = t124+t163+t363;
t379 = t123+t169+t368;
t380 = t124+t168+t367;
t399 = U2.*U4.*t31.*t32.*t295.*t351.*y14;
t400 = U2.*U4.*t31.*t32.*t295.*t352.*y14;
t401 = U2.*U4.*t31.*t32.*t295.*t354.*y14;
t408 = t333+t334+t335;
t409 = t336+t337+t338;
t410 = t339+t340+t341;
t342 = t111+t327;
t343 = t113+t328;
t344 = t112+t329;
t345 = t116+t330;
t346 = t115+t331;
t347 = t117+t332;
t357 = U2.*U4.*t31.*t32.*t293.*t296.*t302;
t358 = U2.*U4.*t31.*t32.*t293.*t296.*t304;
t359 = U2.*U4.*t31.*t32.*t293.*t296.*t305;
FF_med = [y4;y5;y6;t86+t299.*y11;t87+t299.*y12;t88+t299.*y13;-U2.*U4.*t31.*t32.*t293.*t295;y11.*(t177+t309)-t342.*y12-t343.*y13+t348.*y14-t357.*y14;y12.*(t178+t313)-t344.*y11-t345.*y13+t349.*y14-t358.*y14;y13.*(t179+t317)-t346.*y11-t347.*y12+t350.*y14-t359.*y14;t33;t34;t35;-U4.*t17.*t18.*t30.*t31.*t37.*t69.*t293-U4.*t17.*t19.*t30.*t31.*t37.*t69.*t293-U4.*t17.*t20.*t30.*t31.*t37.*t69.*t293];
if nargout > 1
    t372 = t369.*y13;
    t373 = t370.*y12;
    t374 = t371.*y11;
    t381 = t376.*y11;
    t382 = t375.*y12;
    t383 = t378.*y11;
    t384 = t377.*y13;
    t385 = t380.*y12;
    t386 = t379.*y13;
    t387 = U2.*U4.*t31.*t32.*t296.*t300.*t304.*y14;
    t388 = U2.*U4.*t31.*t32.*t296.*t301.*t302.*y14;
    t389 = U2.*U4.*t31.*t32.*t296.*t300.*t305.*y14;
    t390 = U2.*U4.*t31.*t32.*t296.*t302.*t303.*y14;
    t391 = U2.*U4.*t31.*t32.*t296.*t301.*t305.*y14;
    t392 = U2.*U4.*t31.*t32.*t296.*t303.*t304.*y14;
    t405 = U2.*U4.*t31.*t32.*t293.*t297.*t302.*t304.*y14.*2.0;
    t406 = U2.*U4.*t31.*t32.*t293.*t297.*t302.*t305.*y14.*2.0;
    t407 = U2.*U4.*t31.*t32.*t293.*t297.*t304.*t305.*y14.*2.0;
    t411 = U2.*U4.*t31.*t32.*t293.*t296.*t353.*y14;
    t412 = U2.*U4.*t31.*t32.*t293.*t296.*t355.*y14;
    t413 = U2.*U4.*t31.*t32.*t293.*t296.*t356.*y14;
    t393 = -t387;
    t394 = -t388;
    t395 = -t389;
    t396 = -t390;
    t397 = -t391;
    t398 = -t392;
    t414 = -t411;
    t415 = -t412;
    t416 = -t413;
    t417 = t372+t381+t382+t393+t394+t399+t405+t414;
    t418 = t373+t383+t384+t395+t396+t400+t406+t415;
    t419 = t374+t385+t386+t397+t398+t401+t407+t416;
    mt1 = [0.0,0.0,0.0,t85+t110-t309,t342,t343,-t348+t357,-y12.*(t27.*t67.*y2.*3.0+t2.*t8.*t67.*y2.*6.0-t21.*t27.*t68.*y2.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t402.*y12)-y13.*(t27.*t67.*y3.*3.0+t2.*t8.*t67.*y3.*6.0-t21.*t27.*t68.*y3.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t402.*y13)-y11.*(t97+t119+t128+t173-U4.*t17.*t30.*t31.*t36.*t69.*t402.*y11)+U2.*U4.*t31.*t32.*t295.*t402.*y14+U2.*U4.*t31.*t32.*t293.*t297.*t302.^2.*y14.*2.0-U2.*U4.*t31.*t32.*t296.*t300.*t302.*y14.*2.0-U2.*U4.*t31.*t32.*t293.*t296.*y14.*(MP2_3.*t206.*2.0+MP2_2.*t266+MP2_4.*t182.*t206.*6.0+MP2_5.*t186.*t206.*1.2e+1+MP2_3.*t182.*t266.*2.0+MP2_4.*t186.*t266.*3.0+MP2_5.*t187.*t266.*4.0),t417,t418,0.0,0.0,0.0,t408,0.0,0.0,0.0,t344];
    mt2 = [t85+t114-t313,t345,-t349+t358,t417,-y11.*(t28.*t67.*y1.*3.0+t3.*t9.*t67.*y1.*6.0-t22.*t28.*t68.*y1.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t403.*y11)-y13.*(t28.*t67.*y3.*3.0+t3.*t9.*t67.*y3.*6.0-t22.*t28.*t68.*y3.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t403.*y13)-y12.*(t99+t120+t129+t174-U4.*t17.*t30.*t31.*t36.*t69.*t403.*y12)+U2.*U4.*t31.*t32.*t295.*t403.*y14+U2.*U4.*t31.*t32.*t293.*t297.*t304.^2.*y14.*2.0-U2.*U4.*t31.*t32.*t296.*t301.*t304.*y14.*2.0-U2.*U4.*t31.*t32.*t293.*t296.*y14.*(MP2_3.*t207.*2.0+MP2_2.*t267+MP2_4.*t182.*t207.*6.0+MP2_5.*t186.*t207.*1.2e+1+MP2_3.*t182.*t267.*2.0+MP2_4.*t186.*t267.*3.0+MP2_5.*t187.*t267.*4.0),t419,0.0,0.0,0.0,t409,0.0,0.0,0.0,t346,t347,t85+t118-t317];
    mt3 = [-t350+t359,t418,t419,-y11.*(t29.*t67.*y1.*3.0+t4.*t10.*t67.*y1.*6.0-t23.*t29.*t68.*y1.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t404.*y11)-y12.*(t29.*t67.*y2.*3.0+t4.*t10.*t67.*y2.*6.0-t23.*t29.*t68.*y2.*1.5e+1-U4.*t17.*t30.*t31.*t36.*t69.*t404.*y12)-y13.*(t100+t121+t130+t175-U4.*t17.*t30.*t31.*t36.*t69.*t404.*y13)+U2.*U4.*t31.*t32.*t295.*t404.*y14+U2.*U4.*t31.*t32.*t293.*t297.*t305.^2.*y14.*2.0-U2.*U4.*t31.*t32.*t296.*t303.*t305.*y14.*2.0-U2.*U4.*t31.*t32.*t293.*t296.*y14.*(MP2_3.*t208.*2.0+MP2_2.*t268+MP2_4.*t182.*t208.*6.0+MP2_5.*t186.*t208.*1.2e+1+MP2_3.*t182.*t268.*2.0+MP2_4.*t186.*t268.*3.0+MP2_5.*t187.*t268.*4.0),0.0,0.0,0.0,t410,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0];
    mt4 = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t293.*y11,U4.*t17.*t30.*t31.*t37.*t69.*t293.*y12,U4.*t17.*t30.*t31.*t37.*t69.*t293.*y13,0.0,t408,t409,t410,0.0,0.0,0.0,U4.*t17.*t18.*t30.*t31.*t38.*t69.*t293.*2.0+U4.*t17.*t19.*t30.*t31.*t38.*t69.*t293.*2.0+U4.*t17.*t20.*t30.*t31.*t38.*t69.*t293.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t299+U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t293.*y11,U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t293.*y12,U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t293.*y13,0.0];
    mt5 = [t177+t309+y11.*(t306-U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t300.*y11)-U4.*t5.*t14.*t17.*t19.*t30.*t31.*t36.*t70.*t300-U4.*t5.*t14.*t17.*t20.*t30.*t31.*t36.*t70.*t300,t133+t312+y11.*(t307-U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t301.*y11)-U4.*t5.*t14.*t17.*t19.*t30.*t31.*t36.*t70.*t301-U4.*t5.*t14.*t17.*t20.*t30.*t31.*t36.*t70.*t301,t136+t315+y11.*(t308-U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t303.*y11)-U4.*t5.*t14.*t17.*t19.*t30.*t31.*t36.*t70.*t303-U4.*t5.*t14.*t17.*t20.*t30.*t31.*t36.*t70.*t303,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t293.*y11.*-2.0+U4.*t5.*t14.*t17.*t18.*t30.*t31.*t37.*t70.*t293+U4.*t5.*t14.*t17.*t19.*t30.*t31.*t37.*t70.*t293+U4.*t5.*t14.*t17.*t20.*t30.*t31.*t37.*t70.*t293,0.0,0.0,0.0,U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t293.*y11];
    mt6 = [t299+U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t293.*y12,U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t293.*y13,0.0,t132+t310+y12.*(t306-U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t300.*y12)-U4.*t6.*t15.*t17.*t18.*t30.*t31.*t36.*t70.*t300-U4.*t6.*t15.*t17.*t20.*t30.*t31.*t36.*t70.*t300,t178+t313+y12.*(t307-U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t301.*y12)-U4.*t6.*t15.*t17.*t18.*t30.*t31.*t36.*t70.*t301-U4.*t6.*t15.*t17.*t20.*t30.*t31.*t36.*t70.*t301,t138+t316+y12.*(t308-U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t303.*y12)-U4.*t6.*t15.*t17.*t18.*t30.*t31.*t36.*t70.*t303-U4.*t6.*t15.*t17.*t20.*t30.*t31.*t36.*t70.*t303,0.0,0.0,0.0];
    mt7 = [U4.*t17.*t30.*t31.*t37.*t69.*t293.*y12.*-2.0+U4.*t6.*t15.*t17.*t18.*t30.*t31.*t37.*t70.*t293+U4.*t6.*t15.*t17.*t19.*t30.*t31.*t37.*t70.*t293+U4.*t6.*t15.*t17.*t20.*t30.*t31.*t37.*t70.*t293,0.0,0.0,0.0,U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t293.*y11,U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t293.*y12,t299+U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t293.*y13,0.0,t134+t311+y13.*(t306-U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t300.*y13)-U4.*t7.*t16.*t17.*t18.*t30.*t31.*t36.*t70.*t300-U4.*t7.*t16.*t17.*t19.*t30.*t31.*t36.*t70.*t300,t137+t314+y13.*(t307-U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t301.*y13)-U4.*t7.*t16.*t17.*t18.*t30.*t31.*t36.*t70.*t301-U4.*t7.*t16.*t17.*t19.*t30.*t31.*t36.*t70.*t301];
    mt8 = [t179+t317+y13.*(t308-U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t303.*y13)-U4.*t7.*t16.*t17.*t18.*t30.*t31.*t36.*t70.*t303-U4.*t7.*t16.*t17.*t19.*t30.*t31.*t36.*t70.*t303,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t293.*y13.*-2.0+U4.*t7.*t16.*t17.*t18.*t30.*t31.*t37.*t70.*t293+U4.*t7.*t16.*t17.*t19.*t30.*t31.*t37.*t70.*t293+U4.*t7.*t16.*t17.*t20.*t30.*t31.*t37.*t70.*t293,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t348-t357,t349-t358,t350-t359,0.0,0.0,0.0,0.0];
    A_med = reshape([mt1,mt2,mt3,mt4,mt5,mt6,mt7,mt8],14,14);
end
if nargout > 2
    FF_max = [y4;y5;y6;t86+t181.*y11;t87+t181.*y12;t88+t181.*y13;-(U2.*U4.*t31.*t32.*t176)./(MP2_1+MP2_2.*1.2e+2+MP2_3.*1.44e+4+MP2_4.*1.728e+6+MP2_5.*2.0736e+8);t132.*y12+t134.*y13+t177.*y11;t133.*y11+t137.*y13+t178.*y12;t136.*y11+t138.*y12+t179.*y13;t33;t34;t35;-U4.*t17.*t18.*t30.*t31.*t37.*t69.*t176-U4.*t17.*t19.*t30.*t31.*t37.*t69.*t176-U4.*t17.*t20.*t30.*t31.*t37.*t69.*t176];
end
if nargout > 3
    mt9 = [0.0,0.0,0.0,t85+t110,t111,t113,0.0,-y11.*(t97+t119+t128+t173)-t27.*t67.*y2.*y12.*3.0-t27.*t67.*y3.*y13.*3.0-t2.*t8.*t67.*y2.*y12.*6.0-t2.*t8.*t67.*y3.*y13.*6.0+t21.*t27.*t68.*y2.*y12.*1.5e+1+t21.*t27.*t68.*y3.*y13.*1.5e+1,t172+t161.*y11-y12.*(t95-t162)-t3.*t12.*t67.*y11.*3.0,t171+t163.*y11-y13.*(t95-t167)-t4.*t13.*t67.*y11.*3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t112,t85+t114,t116,0.0,t172+t162.*y12-y11.*(t96-t161)-t2.*t11.*t67.*y12.*3.0,-y12.*(t99+t120+t129+t174)-t28.*t67.*y1.*y11.*3.0-t28.*t67.*y3.*y13.*3.0-t3.*t9.*t67.*y1.*y11.*6.0-t3.*t9.*t67.*y3.*y13.*6.0+t22.*t28.*t68.*y1.*y11.*1.5e+1+t22.*t28.*t68.*y3.*y13.*1.5e+1];
    mt10 = [t170+t168.*y12-y13.*(t96-t169)-t4.*t13.*t67.*y12.*3.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t115,t117,t85+t118,0.0,t171+t167.*y13-y11.*(t98-t163)-t2.*t11.*t67.*y13.*3.0,t170+t169.*y13-y12.*(t98-t168)-t3.*t12.*t67.*y13.*3.0,-y13.*(t100+t121+t130+t175)-t29.*t67.*y1.*y11.*3.0-t29.*t67.*y2.*y12.*3.0-t4.*t10.*t67.*y1.*y11.*6.0-t4.*t10.*t67.*y2.*y12.*6.0+t23.*t29.*t68.*y1.*y11.*1.5e+1+t23.*t29.*t68.*y2.*y12.*1.5e+1,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t176.*y11];
    mt11 = [U4.*t17.*t30.*t31.*t37.*t69.*t176.*y12,U4.*t17.*t30.*t31.*t37.*t69.*t176.*y13,0.0,0.0,0.0,0.0,0.0,0.0,0.0,U4.*t17.*t18.*t30.*t31.*t38.*t69.*t176.*2.0+U4.*t17.*t19.*t30.*t31.*t38.*t69.*t176.*2.0+U4.*t17.*t20.*t30.*t31.*t38.*t69.*t176.*2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,t181+U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t176.*y11,U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t176.*y12,U4.*t5.*t14.*t17.*t30.*t31.*t36.*t70.*t176.*y13,0.0,t177,t133,t136,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t176.*y11.*-2.0+U4.*t5.*t14.*t17.*t18.*t30.*t31.*t37.*t70.*t176+U4.*t5.*t14.*t17.*t19.*t30.*t31.*t37.*t70.*t176+U4.*t5.*t14.*t17.*t20.*t30.*t31.*t37.*t70.*t176,0.0,0.0,0.0];
    mt12 = [U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t176.*y11,t181+U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t176.*y12,U4.*t6.*t15.*t17.*t30.*t31.*t36.*t70.*t176.*y13,0.0,t132,t178,t138,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t176.*y12.*-2.0+U4.*t6.*t15.*t17.*t18.*t30.*t31.*t37.*t70.*t176+U4.*t6.*t15.*t17.*t19.*t30.*t31.*t37.*t70.*t176+U4.*t6.*t15.*t17.*t20.*t30.*t31.*t37.*t70.*t176,0.0,0.0,0.0,U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t176.*y11,U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t176.*y12,t181+U4.*t7.*t16.*t17.*t30.*t31.*t36.*t70.*t176.*y13,0.0,t134,t137,t179,0.0,0.0,0.0,U4.*t17.*t30.*t31.*t37.*t69.*t176.*y13.*-2.0+U4.*t7.*t16.*t17.*t18.*t30.*t31.*t37.*t70.*t176+U4.*t7.*t16.*t17.*t19.*t30.*t31.*t37.*t70.*t176+U4.*t7.*t16.*t17.*t20.*t30.*t31.*t37.*t70.*t176,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];
    A_max = reshape([mt9,mt10,mt11,mt12],14,14);
end
