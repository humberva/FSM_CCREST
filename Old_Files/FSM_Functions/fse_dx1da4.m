function dx1da4 = fse_dx1da4(F1,F2,a1,a2,a3,a4,a7,alpha,beta,x1)
%FSE_DX1DA4
%    DX1DA4 = FSE_DX1DA4(F1,F2,A1,A2,A3,A4,A7,ALPHA,BETA,X1)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    23-Aug-2013 13:40:54

t1034 = F2.*a1;
t1028 = F1-t1034;
t1035 = alpha.*t1028;
t1029 = beta-t1035;
t1030 = exp(t1029);
t1031 = t1030+1.0;
t1032 = 1.0./t1031;
t1033 = a3-x1;
t1036 = 1.0./a3;
t1037 = a7+1.0;
t1047 = t1036.*x1;
t1038 = -t1047+1.0;
t1039 = 1.0./t1037;
t1040 = t1038.^t1039;
t1041 = a2-1.0;
t1062 = alpha.*t1033;
t1042 = -beta-t1062;
t1043 = exp(t1042);
t1044 = t1043+1.0;
t1045 = 1.0./t1044;
t1046 = t1032.*t1041.*t1028;
t1048 = a3.*t1040.*t1037;
t1049 = t1046+t1048;
t1063 = alpha.*t1049;
t1050 = beta-t1063;
t1051 = exp(t1050);
t1052 = t1051+1.0;
t1053 = 1.0./t1052;
t1054 = t1032.*t1041.*t1036.*t1028.*t1039;
t1055 = t1040+t1054;
t1056 = t1055.^t1037;
t1057 = a3.*t1053.*t1056;
t1058 = -a3+t1057+x1;
t1059 = t1032.*t1045.*t1058;
t1060 = t1036.*t1028.*x1;
t1061 = t1032.*t1033;
t1064 = t1061+t1059;
t1065 = a3+t1059-x1;
t1074 = alpha.*t1065;
t1066 = beta-t1074;
t1067 = exp(t1066);
t1068 = t1067+1.0;
t1069 = 1.0./t1068;
t1073 = a3.*t1032;
t1075 = t1064.*t1069;
t1070 = t1073-t1075;
t1071 = a4.*t1036.*x1.*(1.0./2.0);
t1072 = t1060+t1071;
t1083 = t1032.*t1072;
t1084 = a4.*t1070.*t1036.*(1.0./2.0);
t1076 = t1060-t1083-t1084+x1;
t1085 = alpha.*t1076;
t1077 = -beta-t1085;
t1078 = exp(t1077);
t1079 = t1070.*t1036.*(1.0./2.0);
t1080 = t1032.*t1036.*x1.*(1.0./2.0);
t1081 = t1080+t1079;
t1082 = t1032-1.0;
t1086 = t1078+1.0;
dx1da4 = (t1081.*t1082)./t1086+alpha.*t1081.*t1082.*t1076.*1.0./t1086.^2.*t1078;
