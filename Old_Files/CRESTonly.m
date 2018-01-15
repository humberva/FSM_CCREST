%Rio Meta 
addpath('/hydros/humberva/Tools/HydroTools/hypro-and-water-libraries-for-matlab/hypro/');

[domain_grid,~] = geotiffread('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/dem_gt30w100n40_RioMeta.tif');

mapinfo = geotiffinfo('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/dem_gt30w100n40_RioMeta.tif');

%% PET data
%PETFolder = '/hydros/humberva/Colombia/EF5_Data/Rio_Meta/PET1km/FAO/';
%RioMeta1km_FAO.PET.10.tif

PETFolder = '/hydros/humberva/Colombia/EF5_Data/Rio_Meta/PET1km/FEWS/';
%RioMeta1km_FEWS.PET.07.tif

%% TRMM data
RainFolder = '/hydros/humberva/Colombia/3B42-V7/Gauge-Corrected/Monthly_Maps_All_1998-2016/1km_rioMeta/'; 
%Resampled_Rainfall_Totals_MAR_2016.tif'

%% CREST Parameters
[PIM,~] = geotiffread('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/imperv_RioMeta.tif');
PIM = PIM./100;
[PWM,~] = geotiffread('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/wm_50mm_param_1km_RioMeta_20170323.tif');
[PFC,~] = geotiffread('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/hydraulic_cond_1km_RioMeta_20170323.tif');
[PB,~] = geotiffread('/hydros/humberva/Colombia/EF5_Data/Rio_Meta/1km/new/b_param_1km_RioMeta_20170323.tif');
PKE = 1;

%% CREST Run
%pre-allocate temporary variables
[nrows,ncols] = size(domain_grid);

tstep = 24*30; %Monthly in hours

ListOfFiles_ET = dir(fullfile('CREST_Outputs/','*aET*'));
ListOfFiles_SM = dir(fullfile('CREST_Outputs/','*SM*'));
ListOfFiles_RO = dir(fullfile('CREST_Outputs/','*Runoff*'));

meses = ['ENE';'FEB';'MAR';'ABR';'MAY';'JUN';'JUL';'AGO';'SEP';'OCT';'NOV';'DIC'];

yyyy_i = 1998;
mm_i = 1;
cSM = zeros(nrows,ncols);

%Create simulation period
cyear = yyyy_i;
cmonth = mm_i;
period = [];
while (cyear < 1999)
    cmonth = cmonth + 1;
    period = [period; datenum(cyear,cmonth,1)];
    [cyear,cmonth,~,~,~,~] = datevec(period(end));
end

tic;
        
for period_i = 1:numel(period)
    [year,month,~,~,~,~] = datevec(period(period_i));
    %State variables
    cER = zeros(nrows,ncols);
    ctemX = zeros(nrows,ncols);
    excessET = zeros(nrows,ncols);
    W0 = zeros(nrows,ncols);

    %Read-in Rain and PET
    [Pin,~] = geotiffread([RainFolder, 'Resampled_Rainfall_Totals_', meses(month,:), '_', num2str(year), '.tif']);
   
    Pin(Pin<0)= NaN;
    %----caET = PET.*1;
    [PET,~] = geotiffread([PETFolder, 'RioMeta1km_FEWS.PET.', num2str(month, '%02.f'), '.tif']);
    %Water Balance
    %if (P(cstep) <= aET(cstep))
    %Convert precipitation from mm/month into mm
    cP = Pin.* 1; %tstep; 

    %Convert PET from mm/hr into mm, and compute actual ET
    caET = PKE.*(single(PET) .* tstep); 

    %This is the precip that makes it to the soil.
    cPSoil = (cP-caET) .* (1 - PIM);

    %Deal with grids that don't get Precip (i.e. PSoil <= 0)
    % cI(cPSoil <= 0) = 0;
    % cER(cPSoil <= 0) = 0;
    % cERO(cPSoil <= 0) = 0;
    % cERI(cPSoil <= 0) = 0;

    excessET(cPSoil <= 0) = (caET(cPSoil <= 0) - cP(cPSoil <= 0)) .* cSM(cPSoil <= 0) ./ PWM(cPSoil <= 0);
    % ctemX(cPSoil <= 0) = (caET(cPSoil <= 0) - cP(cPSoil <= 0)) .* cSM(cPSoil <= 0) ./ PWM(cPSoil <= 0);

    W0((cPSoil <= 0) & (excessET < cSM)) = cSM((cPSoil <= 0) & (excessET < cSM)) - excessET((cPSoil <= 0) & (excessET < cSM));
    W0((cPSoil <= 0) & (excessET >= cSM)) = 0;
    
    excessET((cPSoil <= 0) & (excessET >= cSM)) = cSM((cPSoil <= 0) & (excessET >= cSM));

    % W0((cPSoil <= 0) & (ctemX < cSM)) = cSM((cPSoil <= 0) & (ctemX < cSM)) - ctemX((cPSoil <= 0) & (ctemX < cSM));

    %Originally not commented. It should not be needed since before this point
    %W0 is either equal to cSM (which is less than or equal to WM) or less than
    %cSM
    % W0(W0 > PWM) = PWM(W0 > PWM);

    caET(cPSoil <= 0) = excessET(cPSoil <= 0) + cP(cPSoil <= 0); %cSM(cPSoil <= 0) - W0(cPSoil <= 0);
    cPSoil(cPSoil <= 0) = 0;

    %else

    %Deal with grids that do get Precip (i.e. PSoil > 0)
    %A. Deal with "unsaturated grids"
    vicPars.b = PB((cPSoil > 0) & (cSM < PWM));
    vicPars.Wm = PWM((cPSoil > 0) & (cSM < PWM));
    %cER - Excess Rainfall
    [cI((cPSoil > 0) & (cSM < PWM)),cER((cPSoil > 0) & (cSM < PWM)),W0((cPSoil > 0) & (cSM < PWM)),cperc((cPSoil > 0) & (cSM < PWM))] = dvic(cPSoil((cPSoil > 0) & (cSM < PWM)),vicPars,cSM((cPSoil > 0) & (cSM < PWM)));
    % [~,cER((cPSoil > 0) & (cSM < PWM)),W0((cPSoil > 0) & (cSM < PWM)),~] = dvic(cPSoil((cPSoil > 0) & (cSM < PWM)),vicPars,cSM((cPSoil > 0) & (cSM < PWM)));

    %Get Excess of Soil Mositure for Interflow
    interflowExcess = W0 - PWM;
    interflowExcess(interflowExcess < 0) = 0;

    %B. Deal with "saturated" grids
    cER((cPSoil > 0) & (cSM >= PWM)) = cPSoil((cPSoil > 0) & (cSM >= PWM));
    W0((cPSoil > 0) & (cSM >= PWM)) = PWM((cPSoil > 0) & (cSM >= PWM));
    % cI((cPSoil > 0) & (cSM >= PWM)) = 0;
    % cperc((cPSoil > 0) & (cSM >= PWM)) = 0;

    %ER Partitioning into Overland Flow and Interflow
    %EF5 Implementation
    ctemX(cPSoil > 0) = ((cSM(cPSoil > 0) + W0(cPSoil > 0)) ./ (2.*PWM(cPSoil > 0))) .* PFC(cPSoil > 0); %Calculate how much water can infiltrate
    
%         % ctemX(cPSoil > 0) = (1-((cSM(cPSoil > 0) + W0(cPSoil > 0)) ./ (2.*PWM(cPSoil > 0)))) .* PFC(cPSoil > 0); %Calculate how much water can infiltrate
%         cERI((cPSoil > 0) & (cER <= ctemX) & (rundata.facc < TH)) = cER((cPSoil > 0) & (cER <= ctemX) & (rundata.facc < TH));
%         cERI((cPSoil > 0) & (cER > ctemX) & (rundata.facc < TH)) = ctemX((cPSoil > 0) & (cER > ctemX) & (rundata.facc < TH));
%         cERI(cERI < 0) = 0; %check why getting negative values
%         cERO(cPSoil > 0) = cER(cPSoil > 0) - cERI(cPSoil > 0); % + (cP(cPSoil > 0) - caET(cPSoil > 0)) .* PIM(cPSoil > 0);
%         % cERO(cERO < 0) = 0; %check why getting negative values

    %Values for States after Water Balance
    cSM = W0;
    cSM(cSM >= PWM) = PWM(cSM >= PWM);
    

    cER(cER < 0 | isnan(cER) == 1) = -9999;
    geotiffwrite(['CREST_Outputs/CREST_Runoff_mm_month_', num2str(year), num2str(month, '%02.f'), '.tif'], cER, mapinfo.RefMatrix);

    caET(caET < 0 | isnan(caET) == 1) = -9999;
    geotiffwrite(['CREST_Outputs/CREST_ActualET_mm_month_', num2str(year), num2str(month, '%02.f'), '.tif'], caET, mapinfo.RefMatrix);

    cSM(cSM < 0 | isnan(cSM) == 1) = -9999;
    geotiffwrite(['CREST_Outputs/CREST_SoilMoisture_mm_month_', num2str(year), num2str(month, '%02.f'), '.tif'], cSM, mapinfo.RefMatrix);
end
