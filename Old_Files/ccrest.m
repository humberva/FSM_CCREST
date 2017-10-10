function [prediction] = ccrest(rundata)
%This function runs the Lumped Continuous Coupled Routing and Excess STorage (CCREST) model
%This implementation uses a series of 3 linear tanks for overland routing

%Initialize Variables

%"Unpack" rundata
Pin = rundata.PP;
PET = rundata.PET;
States_0 = rundata.States;
Parameters = rundata.Parameters;


%Simulation Settings
basin_area = rundata.barea; %basin area in km^2
tperiod = rundata.tperiod; %time period for simulation

%Subset data for simulation
st_idx = find(rundata.PET_period == tperiod(1));
en_idx = find(rundata.PET_period == tperiod(end));
PET = PET(st_idx:en_idx);

st_idx = find(rundata.PP_period == tperiod(1));
en_idx = find(rundata.PP_period == tperiod(end));
Pin = Pin(st_idx:en_idx);

%Forcing
F(1,:) = Pin; %Precipitation (mm)
F(2,:) = PET; %Potential Evapotranspiration (mm)

%Model Parameters
a(1) = Parameters(1); %ET Factor
a(2) = Parameters(2)/100; %Impervious Area Percentage
a(3) = Parameters(3); %Soil's Maximum Water Capacity (mm)
a(4) = Parameters(4); %Ksat in (mm)
a(5) = Parameters(5); %Overland Reservoir Depletion Rate
a(6) = Parameters(6); %Interflow Reservoir Depletion Rate
a(7) = Parameters(7); %Exponential in VIC model

%Model States
x_1_P = States_0.IWU; %Soil Saturation (%)
x = zeros(5,length(F));
x(1,1) = a(3).*(x_1_P/100); %Soil Moisture(mm)
x(2,1) = States_0.ISU; %Interflow Reservoir (mm)
x(3:5,1) = States_0.ISO; %Overland Reservoirs (mm)

%Approximation Parameters
alpha1 = rundata.alpha.alpha1;
beta1 = rundata.beta.beta1;
alpha2 = rundata.alpha.alpha2;
beta2 = rundata.beta.beta2;
alpha3 = rundata.alpha.alpha3;
beta3 = rundata.beta.beta3;
alpha5 = rundata.alpha.alpha5;
beta5 = rundata.beta.beta5;
% alpha6 = rundata.alpha.alpha6;
% beta6 = rundata.beta.beta6;
alpha7 = rundata.alpha.alpha7;
beta7 = rundata.beta.beta7;
alpha8 = rundata.alpha.alpha8;
beta8 = rundata.beta.beta8;
% alpha9 = rundata.alpha.alpha9;
% beta9 = rundata.beta.beta9;
% alpha10 = rundata.alpha.alpha10;
% beta10 = rundata.beta.beta10;
alpha11 = rundata.alpha.alpha11;
beta11 = rundata.beta.beta11;
% alpha12 = rundata.alpha.alpha12;
% beta12 = rundata.beta.beta12;

%Process equations
M1 = zeros(1,length(F));
M2 = zeros(1,length(F));
M3 = zeros(1,length(F));
M4 = zeros(1,length(F));
M5 = zeros(1,length(F));
M6 = zeros(1,length(F));
M7 = zeros(1,length(F));
M8 = zeros(1,length(F)); %Streamflow (cms)

NoNL = zeros(7,length(F));

%Enter Loop
for tstep = 1:length(F)
    %Water Balance Computations
    %Soil Precipitation
    cond1 = (F(1,tstep)-a(1)*F(2,tstep)); %P - aET   
    approx_f1 = 1+exp(-alpha1*(cond1)+beta1); NoNL(1,tstep) = 1/approx_f1;
    M1(tstep) = ((F(1,tstep)-a(1)*F(2,tstep))*(1-a(2)))/approx_f1;

    %Infiltration
    cond2 = (a(3)*(1+a(7))*(1-(x(1,tstep)/a(3)))^(1/(1+a(7)))-M1(tstep));
    approx_f2 = 1 + exp(-alpha2*cond2+beta2); NoNL(2,tstep) = 1/approx_f2;
    cond3 = (a(3)-x(1,tstep));
    approx_f3 = 1 + exp(-alpha3*cond3-beta3); NoNL(3,tstep) = 1/approx_f3;
    signo = ((1-(x(1,tstep)/a(3)))^(1/(1+a(7)))-(M1(tstep)/(a(3)*(1+a(7)))))/abs((1-(x(1,tstep)/a(3)))^(1/(1+a(7)))-(M1(tstep)/(a(3)*(1+a(7)))));
    M2(tstep) = (a(3) - x(1,tstep) - (a(3)*signo*abs((1-(x(1,tstep)/a(3)))^(1/(1+a(7)))-(M1(tstep)/(a(3)*(1+a(7)))))^(1+a(7)))/approx_f2)*(1/(approx_f3*approx_f1));
    
    %Excess Rainfall - ER
    approx_f5 = 1 + exp(-alpha5*(M1(tstep) - M2(tstep))+beta5); NoNL(4,tstep) = 1/approx_f5;
    M3(tstep) = (M1(tstep)-M2(tstep))/(approx_f5);
    
    %Soil Moisture Updating - W0 (crest.m)/WA(report) Gaining term
    cond8 = a(3) - x(1,tstep) - M2(tstep);
    approx_f8 = 1 + exp(-alpha8*cond8+beta8); NoNL(5,tstep) = 1/approx_f8;
%     approx_f10 = 1 + exp(alpha10*cond8-beta10);
%     M6(tstep) = ((x(1,tstep)+M2(tstep))/(approx_f8)+(a(3)/(approx_f10)))*(1/approx_f1);
    M6(tstep) = ((x(1,tstep) - a(3))/approx_f1 + M2(tstep))/approx_f8 + a(3)/approx_f1;
    
    %temX - Gain or loss of water in the soil tank
%     approx_f9 = 1 + exp(-alpha9*(a(1)*F(2,tstep)-F(1,tstep))-beta9);
%     M5(tstep) = ((a(4)*(x(1,tstep)+M6(tstep)))/(2*a(3)*approx_f1)) + (((a(1)*F(2,tstep)-F(1,tstep))*x(1,tstep))/(approx_f9*a(3)));
    M5(tstep) = ((a(4)*x(1,tstep))/(2*a(3))-(a(1)*F(2,tstep)-F(1,tstep))*x(1,tstep)/a(3))/approx_f1+M6(tstep)*a(4)/(2*a(3))+(a(1)*F(2,tstep)-F(1,tstep))*x(1,tstep)/a(3);
    
    %Excess Interflow - ERI
%     approx_f6 = 1 + exp(-alpha6*(M5(tstep)-M3(tstep))-beta6);
    approx_f7 = 1 + exp(-alpha7*(M3(tstep)-M5(tstep))+beta7);
    NoNL(6,tstep) = 1/approx_f7;
%     M4(tstep) = ((M3(tstep)/approx_f6) + (M5(tstep)/approx_f7))*(1/approx_f1);
    M4(tstep) = (M5(tstep)-M3(tstep))/approx_f7 + M3(tstep);
    
    %Excess Overland - ERO
    M7(tstep) = M3(tstep) - M4(tstep) + (a(2)*(F(1,tstep)-a(1)*F(2,tstep)))/approx_f1;
    
    %Streamflow - Q
    M8(tstep) = (a(6)*x(2,tstep)+a(5)*x(5,tstep))*basin_area*(1/3.6);
    
    if (M8(tstep) < 0)
        error('');
    end
    
    %Update States + Flow Routing
    %Soil Moisture
    approx_f11 = 1 + exp(-alpha11*(x(1,tstep)-M5(tstep))-beta11);
    NoNL(7,tstep) = 1/approx_f11;
    x(1,tstep+1) = M6(tstep) + ((x(1,tstep)-M5(tstep))/approx_f11)*(1-1/approx_f1);
    
    %Interflow Reservoir
    x(2,tstep+1) = x(2,tstep)*(1-a(6))+M4(tstep);
    
    %Overland Reservoirs
    x(3,tstep+1) = x(3,tstep)*(1-a(5)) + M7(tstep);
    x(4,tstep+1) = x(4,tstep)*(1-a(5)) + a(5)*x(3,tstep);
    x(5,tstep+1) = x(5,tstep)*(1-a(5)) + a(5)*x(4,tstep);
%     if (isreal(x(:,tstep+1))==0)
%         error(['Here at tstep = ', num2str(tstep)]);
%     end
end

prediction.Q = M8;
prediction.SM = x(1,:);
prediction.ISU = x(2,:);
prediction.ISO = x(3:end,:);
prediction.M1 = M1;
prediction.M2 = M2;
prediction.M3 = M3;
prediction.M4 = M4;
prediction.M5 = M5;
prediction.M6 = M6;
prediction.M7 = M7;
prediction.NonLinear = NoNL;