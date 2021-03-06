function dx2dF1 = fse_dx2dF1(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX2DF1
%    DX2DF1 = FSE_DX2DF1(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    23-Aug-2013 13:40:30

t81 = 1.0./a3;
t83 = F2.*a1;
t82 = F1-t83;
t88 = alpha.*t82;
t84 = beta-t88;
t85 = exp(t84);
t86 = t85+1.0;
t87 = 1.0./t86;
t89 = a2-1.0;
t90 = t82.*t87.*t89;
t91 = a7+1.0;
t100 = t81.*x1;
t92 = -t100+1.0;
t93 = 1.0./t91;
t94 = t92.^t93;
t95 = a3-x1;
t114 = alpha.*t95;
t96 = -beta-t114;
t97 = exp(t96);
t98 = t97+1.0;
t99 = 1.0./t98;
t101 = a3.*t91.*t94;
t102 = t101+t90;
t115 = alpha.*t102;
t103 = beta-t115;
t104 = exp(t103);
t105 = t104+1.0;
t106 = 1.0./t105;
t107 = t81.*t82.*t87.*t89.*t93;
t108 = t107+t94;
t109 = t108.^t91;
t110 = a3.*t106.*t109;
t111 = -a3+t110+x1;
t116 = t111.*t87.*t99;
t112 = -t116+t90;
t113 = t81.*t82.*x1;
t117 = alpha.*t112;
t118 = beta+t117;
t119 = exp(t118);
t120 = t119+1.0;
t121 = 1.0./t120;
t122 = t87.*t89;
t123 = 1.0./t86.^2;
t124 = alpha.*t123.*t82.*t85.*t89;
t125 = a3+t116-x1;
t139 = alpha.*t125;
t126 = beta-t139;
t127 = exp(t126);
t128 = t127+1.0;
t129 = 1.0./t128;
t130 = t81.*t87.*t89.*t93;
t131 = alpha.*t123.*t81.*t82.*t85.*t89.*t93;
t132 = t130+t131;
t133 = t108.^a7;
t134 = a3.*t106.*t132.*t133.*t91;
t135 = 1.0./t105.^2;
t136 = t122+t124;
t137 = a3.*alpha.*t104.*t109.*t135.*t136;
t138 = t134+t137;
t140 = t138.*t87.*t99;
t141 = alpha.*t111.*t123.*t85.*t99;
t142 = t87.*t95;
t143 = t116+t142;
t144 = a4.*t81.*x1.*(1.0./2.0);
t145 = t113+t144;
t146 = t122+t124-t140-t141;
t147 = t145.*t87;
t148 = t112.*t121;
t149 = a3.*t87;
t156 = t129.*t143;
t150 = t149-t156;
t151 = a4.*t150.*t81.*(1.0./2.0);
t152 = -t113+t147+t148+t151;
t153 = alpha.*t152;
t154 = beta+t153;
t155 = exp(t154);
t157 = t155+1.0;
t158 = alpha.*t123.*t85.*t95;
t159 = t140+t141+t158;
t160 = t129.*t159;
t161 = t140+t141;
t162 = 1.0./t128.^2;
t163 = alpha.*t127.*t143.*t161.*t162;
t164 = t160+t163-a3.*alpha.*t123.*t85;
t165 = a4.*t164.*t81.*(1.0./2.0);
t166 = 1.0./t120.^2;
t167 = alpha.*t112.*t119.*t146.*t166;
t168 = t100+t165+t167-t121.*t146-t81.*t87.*x1-alpha.*t123.*t145.*t85;
dx2dF1 = t167-t121.*t146-t168./t157+alpha.*t152.*t155.*1.0./t157.^2.*t168;
