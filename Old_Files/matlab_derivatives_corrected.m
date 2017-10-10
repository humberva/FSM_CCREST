clear;

syms x1 x2 x3 x4 x5 a1 a2 a3 a4 a5 a6 a7 F1 F2 alpha beta;

%Water Balance Computations
%Soil Precipitation
cond1 = (F1-a1*F2); %P - aET   
approx_f1 = 1+exp(-alpha*(cond1)+beta);
M1 = ((F1-a1*F2)*(1-a2))/approx_f1;

%Infiltration
cond2 = (a3*(1+a7)*(1-(x1/a3))^(1/(1+a7))-M1);
approx_f2 = 1 + exp(-alpha*cond2+beta);
cond3 = (a3-x1);
approx_f3 = 1 + exp(-alpha*cond3-beta);
M2 = (a3 - x1 - (a3*((1-(x1/a3))^(1/(1+a7))-(M1/(a3*(1+a7))))^(1+a7))/approx_f2)*(1/(approx_f3*approx_f1));

%Excess Rainfall - ER
approx_f5 = 1 + exp(-alpha*(M1 - M2)+beta);
M3 = (M1-M2)/(approx_f5);

%Soil Moisture Updating - W0 (crest.m)/WA(report) Gaining term
cond8 = a3 - x1 - M2;
approx_f8 = 1 + exp(-alpha*cond8+beta);
M6 = ((x1 - a3)/approx_f1 + M2)/approx_f8 + a3/approx_f1;

%temX - Gain or loss of water in the soil tank
M5 = ((a4*x1)/(2*a3)-(a1*F2-F1)*x1/a3)/approx_f1+M6*a4/(2*a3)+(a1*F2-F1)*x1/a3;

%Excess Interflow - ERI
approx_f7 = 1 + exp(-alpha*(M3-M5)+beta);
M4 = (M5-M3)/approx_f7 + M3;

%Excess Overland - ERO
M7 = M3 - M4 + (a2*(F1-a1*F2))/approx_f1;

%Update States + Flow Routing
%Soil Moisture
approx_f11 = 1 + exp(-alpha*(x1-M5)-beta);
newx(1) = M6 + ((x1-M5)/approx_f11)*(1-1/approx_f1);

%Interflow Reservoir
newx(2) = x2*(1-a6)+M4;

%Overland Reservoirs
newx(3) = x3*(1-a5) + M7;
newx(4) = x4*(1-a5) + a5*x3;
newx(5) = x5*(1-a5) + a5*x4;

%EXPRESSIONS FOR DERIVATIVES (Sensitivities)

%Derivatives w.r.t Forcing?
dx1dF1 = diff(newx(1),F1);
matlabFunction(dx1dF1, 'file', 'fse_dx1dF1', 'vars', symvar(dx1dF1));

dx2dF1 = diff(newx(2),F1);
matlabFunction(dx2dF1, 'file', 'fse_dx2dF1', 'vars', symvar(dx2dF1));

dx3dF1 = diff(newx(3),F1);
matlabFunction(dx3dF1, 'file', 'fse_dx3dF1', 'vars', symvar(dx3dF1));

dx1dF2 = diff(newx(1),F2);
matlabFunction(dx1dF2, 'file', 'fse_dx1dF2', 'vars', symvar(dx1dF2));

dx2dF2 = diff(newx(2),F2);
matlabFunction(dx2dF2, 'file', 'fse_dx2dF2', 'vars', symvar(dx2dF2));

dx3dF2 = diff(newx(3),F2);
matlabFunction(dx3dF2, 'file', 'fse_dx3dF2', 'vars', symvar(dx3dF2));

%Derivatives w.r.t. initial conditions
dx1dx1 = diff(newx(1),x1);
matlabFunction(dx1dx1, 'file', 'fse_dx1dx1', 'vars', symvar(dx1dx1));

dx1dx2 = diff(newx(1),x2);
matlabFunction(dx1dx2, 'file', 'fse_dx1dx2', 'vars', symvar(dx1dx2));

dx1dx3 = diff(newx(1),x3);
matlabFunction(dx1dx3, 'file', 'fse_dx1dx3', 'vars', symvar(dx1dx3));

dx1dx4 = diff(newx(1),x4);
matlabFunction(dx1dx4, 'file', 'fse_dx1dx4', 'vars', symvar(dx1dx4));

dx1dx5 = diff(newx(1),x5);
matlabFunction(dx1dx5, 'file', 'fse_dx1dx5', 'vars', symvar(dx1dx5));

dx2dx1 = diff(newx(2),x1);
matlabFunction(dx2dx1, 'file', 'fse_dx2dx1', 'vars', symvar(dx2dx1));

dx2dx2 = diff(newx(2),x2);
matlabFunction(dx2dx2, 'file', 'fse_dx2dx2', 'vars', symvar(dx2dx2));

dx2dx3 = diff(newx(2),x3);
matlabFunction(dx2dx3, 'file', 'fse_dx2dx3', 'vars', symvar(dx2dx3));

dx2dx4 = diff(newx(2),x4);
matlabFunction(dx2dx4, 'file', 'fse_dx2dx4', 'vars', symvar(dx2dx4));

dx2dx5 = diff(newx(2),x5);
matlabFunction(dx2dx5, 'file', 'fse_dx2dx5', 'vars', symvar(dx2dx5));

dx3dx1 = diff(newx(3),x1);
matlabFunction(dx3dx1, 'file', 'fse_dx3dx1', 'vars', symvar(dx3dx1));

dx3dx2 = diff(newx(3),x2);
matlabFunction(dx3dx2, 'file', 'fse_dx3dx2', 'vars', symvar(dx3dx2));

dx3dx3 = diff(newx(3),x3);
matlabFunction(dx3dx3, 'file', 'fse_dx3dx3', 'vars', symvar(dx3dx3));

dx3dx4 = diff(newx(3),x4);
matlabFunction(dx3dx4, 'file', 'fse_dx3dx4', 'vars', symvar(dx3dx4));

dx3dx5 = diff(newx(3),x5);
matlabFunction(dx3dx5, 'file', 'fse_dx3dx5', 'vars', symvar(dx3dx5));

dx4dx1 = diff(newx(4),x1);
matlabFunction(dx4dx1, 'file', 'fse_dx4dx1', 'vars', symvar(dx4dx1));

dx4dx2 = diff(newx(4),x2);
matlabFunction(dx4dx2, 'file', 'fse_dx4dx2', 'vars', symvar(dx4dx2));

dx4dx3 = diff(newx(4),x3);
matlabFunction(dx4dx3, 'file', 'fse_dx4dx3', 'vars', symvar(dx4dx3));

dx4dx4 = diff(newx(4),x4);
matlabFunction(dx4dx4, 'file', 'fse_dx4dx4', 'vars', symvar(dx4dx4));

dx4dx5 = diff(newx(4),x5);
matlabFunction(dx4dx5, 'file', 'fse_dx4dx5', 'vars', symvar(dx4dx5));

dx5dx1 = diff(newx(5),x1);
matlabFunction(dx5dx1, 'file', 'fse_dx5dx1', 'vars', symvar(dx5dx1));

dx5dx2 = diff(newx(5),x2);
matlabFunction(dx5dx2, 'file', 'fse_dx5dx2', 'vars', symvar(dx5dx2));

dx5dx3 = diff(newx(5),x3);
matlabFunction(dx5dx3, 'file', 'fse_dx5dx3', 'vars', symvar(dx5dx3));

dx5dx4 = diff(newx(5),x4);
matlabFunction(dx5dx4, 'file', 'fse_dx5dx4', 'vars', symvar(dx5dx4));

dx5dx5 = diff(newx(5),x5);
matlabFunction(dx5dx5, 'file', 'fse_dx5dx5', 'vars', symvar(dx5dx5));

% Derivatives w.r.t. parameters
dx1da1 = diff(newx(1),a1);
matlabFunction(dx1da1, 'file', 'fse_dx1da1', 'vars', symvar(dx1da1));

dx1da2 = diff(newx(1),a2);
matlabFunction(dx1da2, 'file', 'fse_dx1da2', 'vars', symvar(dx1da2));

dx1da3 = diff(newx(1),a3);
matlabFunction(dx1da3, 'file', 'fse_dx1da3', 'vars', symvar(dx1da3));

dx1da4 = diff(newx(1),a4);
matlabFunction(dx1da4, 'file', 'fse_dx1da4', 'vars', symvar(dx1da4));

dx1da5 = diff(newx(1),a5);
matlabFunction(dx1da5, 'file', 'fse_dx1da5', 'vars', symvar(dx1da5));

dx1da6 = diff(newx(1),a6);
matlabFunction(dx1da6, 'file', 'fse_dx1da6', 'vars', symvar(dx1da6));

dx1da7 = diff(newx(1),a7);
matlabFunction(dx1da7, 'file', 'fse_dx1da7', 'vars', symvar(dx1da7));

dx2da1 = diff(newx(2),a1);
matlabFunction(dx2da1, 'file', 'fse_dx2da1', 'vars', symvar(dx2da1));

dx2da2 = diff(newx(2),a2);
matlabFunction(dx2da2, 'file', 'fse_dx2da2', 'vars', symvar(dx2da2));

dx2da3 = diff(newx(2),a3);
matlabFunction(dx2da3, 'file', 'fse_dx2da3', 'vars', symvar(dx2da3));

dx2da4 = diff(newx(2),a4);
matlabFunction(dx2da4, 'file', 'fse_dx2da4', 'vars', symvar(dx2da4));

dx2da5 = diff(newx(2),a5);
matlabFunction(dx2da5, 'file', 'fse_dx2da5', 'vars', symvar(dx2da5));

dx2da6 = diff(newx(2),a6);
matlabFunction(dx2da6, 'file', 'fse_dx2da6', 'vars', symvar(dx2da6));

dx2da7 = diff(newx(2),a7);
matlabFunction(dx2da7, 'file', 'fse_dx2da7', 'vars', symvar(dx2da7));

dx3da1 = diff(newx(3),a1);
matlabFunction(dx3da1, 'file', 'fse_dx3da1', 'vars', symvar(dx3da1));

dx3da2 = diff(newx(3),a2);
matlabFunction(dx3da2, 'file', 'fse_dx3da2', 'vars', symvar(dx3da2));

dx3da3 = diff(newx(3),a3);
matlabFunction(dx3da3, 'file', 'fse_dx3da3', 'vars', symvar(dx3da3));

dx3da4 = diff(newx(3),a4);
matlabFunction(dx3da4, 'file', 'fse_dx3da4', 'vars', symvar(dx3da4));

dx3da5 = diff(newx(3),a5);
matlabFunction(dx3da5, 'file', 'fse_dx3da5', 'vars', symvar(dx3da5));

dx3da6 = diff(newx(3),a6);
matlabFunction(dx3da6, 'file', 'fse_dx3da6', 'vars', symvar(dx3da6));

dx3da7 = diff(newx(3),a7);
matlabFunction(dx3da7, 'file', 'fse_dx3da7', 'vars', symvar(dx3da7));

dx4da1 = diff(newx(4),a1);
matlabFunction(dx4da1, 'file', 'fse_dx4da1', 'vars', symvar(dx4da1));

dx4da2 = diff(newx(4),a2);
matlabFunction(dx4da2, 'file', 'fse_dx4da2', 'vars', symvar(dx4da2));

dx4da3 = diff(newx(4),a3);
matlabFunction(dx4da3, 'file', 'fse_dx4da3', 'vars', symvar(dx4da3));

dx4da4 = diff(newx(4),a4);
matlabFunction(dx4da4, 'file', 'fse_dx4da4', 'vars', symvar(dx4da4));

dx4da5 = diff(newx(4),a5);
matlabFunction(dx4da5, 'file', 'fse_dx4da5', 'vars', symvar(dx4da5));

dx4da6 = diff(newx(4),a6);
matlabFunction(dx4da6, 'file', 'fse_dx4da6', 'vars', symvar(dx4da6));

dx4da7 = diff(newx(4),a7);
matlabFunction(dx4da7, 'file', 'fse_dx4da7', 'vars', symvar(dx4da7));

dx5da1 = diff(newx(5),a1);
matlabFunction(dx5da1, 'file', 'fse_dx5da1', 'vars', symvar(dx5da1));

dx5da2 = diff(newx(5),a2);
matlabFunction(dx5da2, 'file', 'fse_dx5da2', 'vars', symvar(dx5da2));

dx5da3 = diff(newx(5),a3);
matlabFunction(dx5da3, 'file', 'fse_dx5da3', 'vars', symvar(dx5da3));

dx5da4 = diff(newx(5),a4);
matlabFunction(dx5da4, 'file', 'fse_dx5da4', 'vars', symvar(dx5da4));

dx5da5 = diff(newx(5),a5);
matlabFunction(dx5da5, 'file', 'fse_dx5da5', 'vars', symvar(dx5da5));

dx5da6 = diff(newx(5),a6);
matlabFunction(dx5da6, 'file', 'fse_dx5da6', 'vars', symvar(dx5da6));

dx5da7 = diff(newx(5),a7);
matlabFunction(dx5da7, 'file', 'fse_dx5da7', 'vars', symvar(dx5da7));