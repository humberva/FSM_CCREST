%CONUS domain
clc;
clear;

[domain_grid,~] = geotiffread('/Users/humberto.vergara/Documents/Trabajo/FLASH_Work/FLASH_parameters/CONUS_FlowGrids_06302016/flow_grids_092214/dem_new_mrms_domain.tif');
[rain,~] = geotiffread('/Users/humberto.vergara/Documents/Trabajo/Geomorphlogical_Study_CONUS/Database/Hydrologic_Spatial_Variability_Algorithms/HydroClimatology/PRISM_ppt_30yr_annual_1km_mrms_grid.tif');

mapinfo = geotiffinfo('/Users/humberto.vergara/Documents/Trabajo/FLASH_Work/FLASH_parameters/CONUS_FlowGrids_06302016/flow_grids_092214/dem_new_mrms_domain.tif');

%% PET data
PETFolder = '/Users/humberto.vergara/Documents/Trabajo/FLASH_Work/FLASH_parameters/';

%% CREST Parameters
ParamFolder = '/Users/humberto.vergara/Documents/Trabajo/FLASH_Work/FLASH_parameters/';
[PIM,~] = geotiffread([ParamFolder, 'im_usa.tif']);
PIM = double(PIM)./100;
[PWM,~] = geotiffread([ParamFolder, 'wm_usa.tif']);
PWM = double(PWM);
[PFC,~] = geotiffread([ParamFolder, 'ksat_usa.tif']);
PFC = double(PFC);
[PB,~] = geotiffread([ParamFolder, 'b_usa.tif']);
PB = double(PB);
parameters.PKE = 1;

compgrids = find(rain>0 & PIM >= 0 & PWM > 0 & PB > 0 & PFC > 0);

parameters.PIM = PIM(compgrids);
parameters.PWM = PWM(compgrids);
parameters.PFC = PFC(compgrids);
parameters.PB = PB(compgrids);

%% CREST Run
%pre-allocate temporary variables
[nrows,ncols] = size(domain_grid);

stepHours = 1; %Monthly in hours

meses = ['ENE';'FEB';'MAR';'ABR';'MAY';'JUN';'JUL';'AGO';'SEP';'OCT';'NOV';'DIC'];

yyyy_i = 1998;
mm_i = 7;

%Create simulation period
period = datenum(yyyy_i,mm_i,1):1/24:datenum(yyyy_i,mm_i,2);

[static_PET,~] = geotiffread([PETFolder, 'PET_', num2str(mm_i, '%02.f'), '_usa.tif']);
static_PET = double(static_PET)./24;

%Cold start
SM = zeros(nrows,ncols);
states.SM = SM(compgrids);

tic;

%Spatially Uniform Rain Rate
rain_rate = 50;
syntheticRain = zeros(nrows,ncols)+rain_rate;

time2saturation = zeros(nrows,ncols);

for period_i = 1:numel(period)
    [year,month,~,~,~,~] = datevec(period(period_i));

    %Read-in Rain and PET
    Pin = syntheticRain(compgrids);
    PET = static_PET(compgrids);

    %CREST
%     cSM1 = zeros(nrows,ncols);
%     aET1 = zeros(nrows,ncols);
%     cERI1 = zeros(nrows,ncols);
%     cERO1 = zeros(nrows,ncols);
%     infiltration1 = zeros(nrows,ncols);
%     for pix = 1:numel(compgrids)
%         parameters1.PIM = PIM(compgrids(pix));
%         parameters1.PWM = PWM(compgrids(pix));
%         parameters1.PFC = PFC(compgrids(pix));
%         parameters1.PB = PB(compgrids(pix));
%         parameters1.PKE = 1;
%         states1.SM = SM(compgrids(pix));
%         [cSM1(compgrids(pix)), aET1(compgrids(pix)), cERI1(compgrids(pix)), cERO1(compgrids(pix)), infiltration1(compgrids(pix))] = CRESTef5(stepHours, syntheticRain(compgrids(pix)), static_PET(compgrids(pix)), parameters1, states1);        
%     end 
    
    %[vec_cSM, vec_aET, vec_cERI, vec_cERO, vec_infiltration] = matCRESTef5(stepHours, Pin, PET, parameters, states);
    [vec_cSM, vec_aET, vec_cERI, vec_cERO, vec_infiltration] = ccrest(stepHours, Pin, PET, parameters, states);
    states.SM = vec_cSM;
    
    vec_saturated = vec_cSM./parameters.PWM;
    
    time2saturation(compgrids(vec_saturated < 1)) = time2saturation(compgrids(vec_saturated < 1)) + 1;
    
%     cSM = nan(nrows,ncols);
%     aET = nan(nrows,ncols);
%     cERI = nan(nrows,ncols);
%     cERO = nan(nrows,ncols);
%     infiltration = nan(nrows,ncols);
% 
%     infiltration(compgrids) = vec_infiltration;
%     cSM(compgrids) = vec_cSM;
%     aET(compgrids) = vec_aET;
%     cERI(compgrids) = vec_cERI;
%     cERO(compgrids) = vec_cERO;
%     
%     imagesc(cSM./PWM);
%     caxis([0 1]);
%     title(datestr(period(period_i)));
%     getframe();
end

% cSM1 = round(cSM1.*10000)./10000;
% aET1 = round(aET1.*10000)./10000;
% cERI1 = round(cERI1.*10000)./10000;
% cERO1 = round(cERO1.*10000)./10000;
% infiltration1 = round(infiltration1.*10000)./10000;

cSM = nan(nrows,ncols);
aET = nan(nrows,ncols);
cERI = nan(nrows,ncols);
cERO = nan(nrows,ncols);
infiltration = nan(nrows,ncols);

infiltration(compgrids) = vec_infiltration;
cSM(compgrids) = vec_cSM;
aET(compgrids) = vec_aET;
cERI(compgrids) = vec_cERI;
cERO(compgrids) = vec_cERO;

% infiltration(compgrids) = round(vec_infiltration.*10000)./10000;
% cSM(compgrids) = round(vec_cSM.*10000)./10000;
% aET(compgrids) = round(vec_aET.*10000)./10000;
% cERI(compgrids) = round(vec_cERI.*10000)./10000;
% cERO(compgrids) = round(vec_cERO.*10000)./10000;
toc;