function [cSM, aET, cERI, cERO, infiltration, misc] = ccrest(stepHours, precipIn, petIn, parameters, states)
%This function runs the Continuous Coupled Routing and Excess STorage (CCREST) model
%This version differs from its predecessor in that the equations match
%those in EF5 (ef5.ou.edu) as of January 2018.
%
%No routing equations are included herein; water budget calculations and
%runoff partition into surface and subsurface components alone.

%Initialize Variables                                 
precip = precipIn .* stepHours; % precipIn is mm/hr, precip is mm
pet = petIn .* stepHours; % petIn in mm/hr, pet is mm
  
%Sigmoid function parameters
alpha = 1000000;
beta = 1000;

adjPET = parameters.PKE.*pet;

cond1 = (precip-adjPET);   
approx_f1 = 1./(1+exp(-alpha.*(cond1)+beta));

%cond2 = ExcessET
cond2 = (adjPET - precip).*(states.SM./parameters.PWM).*(1-approx_f1);
approx_f2 = 1./(1+exp(-alpha.*(cond2)+beta));

% M1 = aET
aET = adjPET.*approx_f1 + cond2.*(1-approx_f2) + (states.SM.*approx_f2 + precip).*(1-approx_f1);

% Soil Precipitation: This is the precip that makes it to the soil.
% M2 = cPSoil
M2 = ((precip-adjPET).*(1-parameters.PIM)).*approx_f1;

% Infiltration
im = parameters.PWM.*(1+parameters.PB);
it = im.*(1-(1-states.SM./parameters.PWM).^(1./(1+parameters.PB)));

cond3 =(im - (M2+it));
approx_f3 = 1./(1 + exp(-alpha.*cond3+beta));
cond4 = (parameters.PWM - states.SM);
approx_f4 = 1./(1 + exp(-alpha.*cond4+beta));

% M3 = I
infiltration = ((parameters.PWM - states.SM).*(1-approx_f3) + parameters.PWM.*((1-(it./im)).^(1+parameters.PB)-(1-((it+M2)./im)).^(1+parameters.PB)).*approx_f3).*approx_f4.*approx_f1;

%Excess Rainfall/Infiltration - R
approx_f5 = 1./(1 + exp(-alpha.*(M2 - infiltration)+beta));
M4 = (M2-infiltration).*(approx_f5);

% Soil Moisture Updating - W0 (crest.m)/WA(report) Gaining term
% M5 = Gaining water
M5 = parameters.PWM.*(1-approx_f4 + approx_f4.*approx_f3) + (states.SM + infiltration).*approx_f4.*(1-approx_f3);

% M6 = Losing water
M6 = (states.SM - cond2).*(1-approx_f2);

% newx(1) = New Soil Moisture
cSM = M5.*approx_f1 + M6.*(1-approx_f1);

%Partition to interflow and overland runoff
termX = ((states.SM+M5)./(2.*parameters.PWM)).*parameters.PFC.*approx_f1;

cond6 = (M4-termX);
approx_f6 = 1./(1 + exp(-alpha.*cond6+beta));

%Interflow water
cERI = (termX.*approx_f6 + M4.*(1-approx_f6)).*approx_f1;

%Overland water
cERO = M4 - cERI + (precip - adjPET).*parameters.PIM.*approx_f1;

misc.it = it;
misc.im = im;
misc.approx_f3 = approx_f3;
misc.approx_f4 = approx_f4;