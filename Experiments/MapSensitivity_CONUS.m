clear;
clc;

par_folder = '/Users/humbertovergara/Documents/Trabajo/Geomorphlogical_Study_CONUS/';
morphdata = dlmread([par_folder, 'All_USGS_Morphometry_and_Stream_Network_Variables_latest_basinmorph_version_Oct2014.csv'], ',', 1,0);

surfdata = dlmread([par_folder, 'All_USGS_Surface_Parameters_Revised_020216.csv'], ',',1,0);
surfdata(surfdata<0) = NaN;

climdata = dlmread([par_folder, 'All_USGS_Climate_Precipitation_and_Temperature.csv'], ',', 1,0);
climdata(climdata(:,7)<0,7) = NaN;

qiconususgs = dlmread([par_folder, 'Dams_and_Lakes_QI/All_CONUS_USGS_PeakFlow_Qualification_Code.csv'], ',',1,0);
hsrh = dlmread('~/Documents/Trabajo/Geomorphlogical_Study_CONUS/All_USGS_HSRH.csv', ',',1,0);

climo_gage2 = dlmread('/Users/humbertovergara/Documents/Trabajo/Geomorphlogical_Study_CONUS/basinchar_and_report_sept_2011/ClimoData_GAGE2_CONUS_stations_FLASHdatabase.csv', ',',1,0);

Solratio = dlmread('CONUS_Sensitivity_Tests_Variance_Ratio_Summary_Table.csv', ',');
Soldiff = dlmread('CONUS_Sensitivity_Tests_Variance_Difference_Summary_Table.csv', ',');

basin_idx = nan(size(Solratio,1),1);
for i = 1:size(Solratio,1)
    basin_idx(i) = find(morphdata(:,1) == Solratio(i,1));
end

% subplot(1,2,1);
% usmapgeo(morphdata(basin_idx,5),morphdata(basin_idx,6), surfdata(basin_idx,6), 30, [],'quant', 'true');
% 
% subplot(1,2,2);
% usmapgeo(morphdata(basin_idx,5),morphdata(basin_idx,6), Solratio(:,2), 30, [],'quant', 'true');
% hbar = colorbar;
% set(hbar, 'FontSize', 14)
% ylabel(hbar, 'Percentile of MC to FS Variance Ratio', 'FontSize', 16);

qbias = dlmread('/Users/humbertovergara/Documents/Trabajo/FLASH_Work/ProFLASH/Study1_PostProcessor/daily_max_version/Corrected_Qratio_stats_all_stations_compiled.csv', ',', 1,0);
% raw_qbias = dlmread('/Users/humbertovergara/Documents/Trabajo/FLASH_Work/ProFLASH/Study1_PostProcessor/daily_max_version/Regression_Model_MD_LogMaxQ_Ratio.csv', ',', 1,0);
% qbias = nan(size(raw_qbias,1),2);
% qbias(:,1) = raw_qbias(:,3);
% qbias(:,2) = raw_qbias(:,2);

log_meanQbias = log(qbias(:,2));
[B,IX] = sort(log_meanQbias, 'ascend');

usmapgeo(morphdata(qbias(IX,1),5),morphdata(qbias(IX,1),6), log(qbias(IX,2)), 30, [],'quan', 'true');

sens_qbias = nan(1,1000);
for i = 1:1000
    idx = find(basin_idx(i) == (qbias(:,1)));
    if (isempty(idx) == 0)
        sens_qbias(i) = qbias(idx,4);
    end
end



%% Training dataset
% qc_idx = find(climo_gage2(qbias(:,1),14) >= 0 & climo_gage2(qbias(:,1),14) < 30);
qc_idx = find(qbias(:,1) > 0);
% X = [log(morphdata(qbias(qc_idx,1),2)), log(climdata(qbias(qc_idx,1),7)), climdata(qbias(qc_idx,1),8)];
X = [log(morphdata(qbias(qc_idx,1),2)), log(climdata(qbias(qc_idx,1),7)), climdata(qbias(qc_idx,1),8), sqrt(surfdata(qbias(qc_idx,1),6))];

% y = log(qbias(qc_idx,2));
y = qbias(qc_idx,2);

pred_y = simpleBiasCorrection_nnet(X');