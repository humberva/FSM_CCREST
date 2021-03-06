function dLM1da2 = fse_dLM1da2(F1,F2,a1,alpha,beta)
%FSE_DLM1DA2
%    DLM1DA2 = FSE_DLM1DA2(F1,F2,A1,ALPHA,BETA,LAMBDA2)

%    This function was generated by the Symbolic Math Toolbox version 5.6.
%    10-Jan-2014 14:29:43

syms lambda2;

t123 = F1-F2.*a1;
dLM1da2 = -lambda2-t123./(exp(beta-alpha.*t123)+1.0);
