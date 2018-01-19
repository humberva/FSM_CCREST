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

%% CREST Run
%pre-allocate temporary variables
[nrows,ncols] = size(domain_grid);

stepHours = 1; %Monthly in hours

meses = ['ENE';'FEB';'MAR';'ABR';'MAY';'JUN';'JUL';'AGO';'SEP';'OCT';'NOV';'DIC'];

yyyy_i = 1998;
mm_i = 7;

%Create simulation period
period = datenum(yyyy_i,mm_i,1):1/24:datenum(yyyy_i,mm_i,1);

n_steps = numel(period);

[static_PET,~] = geotiffread([PETFolder, 'PET_', num2str(mm_i, '%02.f'), '_usa.tif']);
static_PET = zeros(nrows,ncols)+0.9999; %double(static_PET)./24;

%Select pixels for computations
compgrids = find(static_PET > 0 & rain>0 & PIM >= 0 & PWM > 0 & PB > 0 & PFC > 0);
n_compgrids = numel(compgrids);

parameters.PIM = PIM(compgrids);
parameters.PWM = PWM(compgrids);
parameters.PFC = PFC(compgrids);
parameters.PB = PB(compgrids);

%Cold start
pctSM = zeros(nrows,ncols)+0.;
SM = pctSM.*PWM;

ccrest_out.cSM = zeros(n_steps,n_compgrids);
ccrest_out.aET = zeros(n_steps,n_compgrids);
ccrest_out.cERI = zeros(n_steps,n_compgrids);
ccrest_out.cERO = zeros(n_steps,n_compgrids);
ccrest_out.infiltration = zeros(n_steps,n_compgrids);

ccrest_states.SM = SM(compgrids);

crest_out.cSM = zeros(n_steps,n_compgrids);
crest_out.aET = zeros(n_steps,n_compgrids);
crest_out.cERI = zeros(n_steps,n_compgrids);
crest_out.cERO = zeros(n_steps,n_compgrids);
crest_out.infiltration = zeros(n_steps,n_compgrids);

crest_states.SM = SM(compgrids);
tic;

%Spatially Uniform Rain Rate
rain_rate = 1;
syntheticRain = zeros(nrows,ncols)+rain_rate;

time2saturation = zeros(nrows,ncols);

for period_i = 1:numel(period)
    [year,month,~,~,~,~] = datevec(period(period_i));

    %Read-in Rain and PET
    Pin = syntheticRain(compgrids);
    PET = static_PET(compgrids);

    %CREST
    [vec_cSM, vec_aET, vec_cERI, vec_cERO, vec_infiltration, misc] = matCRESTef5(stepHours, Pin, PET, parameters, crest_states);
    crest_out.cSM(period_i,:) = vec_cSM;
    crest_out.aET(period_i,:) = vec_aET;
    crest_out.cERI(period_i,:) = vec_cERI;
    crest_out.cERO(period_i,:) = vec_cERO;
    crest_out.infiltration(period_i,:) = vec_infiltration;
    
    crest_states.SM = vec_cSM;
    
    [vec_cSM, vec_aET, vec_cERI, vec_cERO, vec_infiltration, cc_misc] = ccrest(stepHours, Pin, PET, parameters, ccrest_states);
    ccrest_out.cSM(period_i,:) = vec_cSM;
    ccrest_out.aET(period_i,:) = vec_aET;
    ccrest_out.cERI(period_i,:) = vec_cERI;
    ccrest_out.cERO(period_i,:) = vec_cERO;
    ccrest_out.infiltration(period_i,:) = vec_infiltration;
    
    ccrest_states.SM = vec_cSM;
    
    vec_saturated = vec_cSM./parameters.PWM;
    
%     time2saturation(compgrids(vec_saturated < 1)) = time2saturation(compgrids(vec_saturated < 1)) + 1;
    
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

% cSM = nan(nrows,ncols);
% aET = nan(nrows,ncols);
% cERI = nan(nrows,ncols);
% cERO = nan(nrows,ncols);
% infiltration = nan(nrows,ncols);
% 
% infiltration(compgrids) = vec_infiltration;
% cSM(compgrids) = vec_cSM;
% aET(compgrids) = vec_aET;
% cERI(compgrids) = vec_cERI;
% cERO(compgrids) = vec_cERO;

% infiltration(compgrids) = round(vec_infiltration.*10000)./10000;
% cSM(compgrids) = round(vec_cSM.*10000)./10000;
% aET(compgrids) = round(vec_aET.*10000)./10000;
% cERI(compgrids) = round(vec_cERI.*10000)./10000;
% cERO(compgrids) = round(vec_cERO.*10000)./10000;
toc;
adjPET = parameters.PKE.*PET;
PSoil = ((Pin-adjPET).*(1-parameters.PIM))';
target_ix = find(Pin > adjPET & SM(compgrids) < parameters.PWM & PSoil'+misc.A < misc.Wmaxm);
numel(find(crest_out.infiltration(target_ix) > PSoil(target_ix)))