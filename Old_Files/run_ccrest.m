clc; clear;

%Input Data: Ft.Cobb basin
data = dlmread('ts.07325800.crest.csv', ',', 1,1);
barea = 342;  %Basin Area in km^2

%Simulation Data and Settings
rundata.PET = data(:,4);
rundata.PET_period = datenum(2006,10,1,0,0,0):1/24:datenum(2007,10,01,0,0,0);
rundata.PP = data(:,3);
rundata.PP_period = datenum(2006,10,1,0,0,0):1/24:datenum(2007,10,01,0,0,0);
rundata.obsQ = data(:,2);
rundata.obsQ_period = datenum(2006,10,1,0,0,0):1/24:datenum(2007,10,01,0,0,0);
runset.barea = barea;
runset.tstep = 1; %Time Step in hours
runset.tperiod = datenum(2006,10,1,0,0,0):1/24:datenum(2007,10,01,0,0,0);
rundata.runset = runset;


%CREST Parameters
Parameters.PKE = 0.5099; %Multiplier to convert PET
Parameters.PIM = 14.2402; %The impervious area ratio. This ratio is stored as a percentage between 0 and 100.
Parameters.PWM = 108.1742; %The maximum soil water capacity (depth integrated pore space) of the soil layer. (mm)
Parameters.PFC = 9.7051; %The soil saturated hydraulic conductivity (Ksat ,mm/hr).
Parameters.LEAKO = 0.2977; %The overland reservoir discharge multiplier.
Parameters.LEAKI = 0.0060; %The interflow reservoir discharge multiplier.
Parameters.PB = 1; %The exponent of the variable infiltration curve.

%NASH Routing Parameters
Parameters.Nqo = 3;

%CCREST Approximation Parameters
alpha.alpha1 = 10^12;
beta.beta1 = 10^3;
alpha.alpha2 = 10^12;
beta.beta2 = 10^3;
alpha.alpha3 = 10^12;
beta.beta3 = 10^3;
alpha.alpha5 = 10^12;
beta.beta5 = 10^3;
alpha.alpha6 = 10^12;
beta.beta6 = 10^3;
alpha.alpha7 = 10^12; %sensitive when 10^0
beta.beta7 = 10^3;
alpha.alpha8 = 10^12; %sensitive
beta.beta8 = 10^3;
alpha.alpha9 = 10^12; %sensitive
beta.beta9 = 10^3;
alpha.alpha10 = 10^12;
beta.beta10 = 10^3;
alpha.alpha11 = 10^12;
beta.beta11 = 10^3;
alpha.alpha12 = 10^12;
beta.beta12 = 10^3;

Parameters.alpha = alpha;
Parameters.beta = beta;

%CREST States
States_0.IWU = 1;
States_0.ISO = [2 2 2];
States_0.ISU = 2;
rundata.States = States_0;

%Run CCREST
prediction = ccrest(rundata,Parameters);

%Run CREST
ref_prediction = crest(rundata,Parameters);

%Plot and Stats
st_idx = find(rundata.obsQ_period == runset.tperiod(1));
en_idx = find(rundata.obsQ_period == runset.tperiod(end));
rundata.obsQ = rundata.obsQ(st_idx:en_idx);
[bias(1),rmse(1),nsce(1)] = nanhydrostat(rundata.obsQ',ref_prediction.Q);
[bias(2),rmse(2),nsce(2)] = nanhydrostat(rundata.obsQ',prediction.Q);
[bias(3),rmse(3),nsce(3)] = nanhydrostat(ref_prediction.Q,prediction.Q);

plot(rundata.obsQ, 'DisplayName', 'Gauge', 'LineStyle', '--', 'Color', 'k'); hold all;
plot(ref_prediction.Q, 'DisplayName', 'CREST', 'Color', 'b');
plot(prediction.Q, 'DisplayName', 'CCREST', 'Color', 'r');
