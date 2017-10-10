%Process Multi-run 
clear;
clc;
load('SERVER_Results/Sensitivity_Characterization_TS_Alberto.mat')
P2 = load('SERVER_Results/Sensitivity_Characterization_TS_Alberto_Rain_times_2.mat');
P025 = load('SERVER_Results/Sensitivity_Characterization_TS_Alberto_Rain_times_0.25.mat');
load('SERVER_Results/Subset_1k_stations_parameters.mat');
load('SERVER_Results/Subset_1k_stations_aET_P_ratio_10years.mat');

bidx = find(uscrestpar(:,2) > 0); % & uscrestpar(:,2) < 5.5);
uscrestpar = uscrestpar(bidx,:);
usbasinlags = usbasinlags(bidx,:);
usbasinmorph = usbasinmorph(bidx,:);

event_start = datenum(2006,6,6,0,0,0);

%Filter out basins with Wmax > 200 - Some pixels had 2500 values (open
%water) that were averaged out and thus the estimates are negatively
%affected.
good_idx = find(uscrestpar(:,1)<=200);

%Initialize super matrix
X = [];

%U
U_CCmatrix = nan(5*7);
%Matrix indices of zero sensitivity
zero_sens = [21,22,26,28,29,30];

cont = 1;
for i = 1:5
    for j = 1:7
        subplot(5,7,cont);
        if (numel(find(abs(sub2ind([5 7], i,j)-zero_sens) == 0)) > 0)
            set(gca, 'Xlim', [-130 -65], 'Ylim', [25 52], 'Xtick', [], 'YTick', []);
            cont = cont + 1;
            continue;
        end
%         hist(reshape(CumSum_U(i,j,good_idx),1,numel(good_idx)));
%         hist(reshape(FreqDS_U(i,j,good_idx),1,numel(good_idx)));
%         hist(reshape(MeanU(i,j,good_idx),1,numel(good_idx)));
%         hist(reshape(PeakT_U(i,j,good_idx)-event_start,1,numel(good_idx)));
%         hist(reshape(Peak_U(i,j,good_idx),1,numel(good_idx)));
%         hist(reshape(VarCoeff_U(i,j,good_idx),1,numel(good_idx)));
%         hist(MeanSM(good_idx));

        % Maps
%         CVAR = reshape(Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4)));
%         CVAR = reshape(PeakT_U(i,j,good_idx)-event_start,1,numel(good_idx));
%         CVAR = reshape(MeanU(i,j,good_idx),1,numel(good_idx));
%         CVAR = reshape(CumSum_U(i,j,good_idx),1,numel(good_idx));
%         CVAR = reshape(FreqDS_U(i,j,good_idx),1,numel(good_idx));
%         CVAR = reshape(VarCoeff_U(i,j,good_idx),1,numel(good_idx));
        CVAR = reshape(P2.Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4)))./reshape(P025.Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4)));
%         usmapgeo(usbasinmorph(good_idx,3),usbasinmorph(good_idx,4),reshape(Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4))),20,[],'quant', 'false');
        [B,IX] = sort(CVAR, 'ascend');
        X(1:numel(B),1) = 1:numel(B); 
        Qtiles = (X./numel(B)).*100;
        
%         scatter(usbasinmorph(good_idx(IX),4),usbasinmorph(good_idx(IX),3),8,Qtiles, 'filled');
        scatter(usbasinmorph(good_idx(IX),4),usbasinmorph(good_idx(IX),3),8,CVAR, 'filled');
        colormap('jet');
%         caxis([0 100]);
        caxis([0 5]);
        set(gca, 'Xlim', [-130 -65], 'Ylim', [25 52], 'Xtick', [], 'YTick', []);

%         scatter(reshape(P2.Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4)))./reshape(Peak_U(i,j,good_idx),1,numel(usbasinmorph(good_idx,4))), P2.Peak_SM(good_idx)./Peak_SM(good_idx), 10, 'filled');

        %Create Super Matrix
%         X = [X, reshape(Peak_U(i,j,good_idx),numel(good_idx),1)];
        
        %Create covariance/correlation matrices
        U_CCmatrix(i,j) = corr(reshape(CumSum_U(i,j,good_idx),1,numel(good_idx))',reshape(MeanU(i,j,good_idx),1,numel(good_idx))', 'Type', 'Pearson');
        cont = cont + 1;
    end
end
set(gcf, 'Position', [1,74,1437,731]);

% %% Analytics
% close;
% y = Peak_Q(good_idx)'; %Peak_SM(good_idx)'; %Peak_Q(good_idx)';
% [b,se,pval,inmodel,stats,nextstep,history] = stepwisefit(X,y);
% yfit = X(stats.wasnan==0,:)*stats.B + stats.intercept;
% % scatter(y(stats.wasnan == 0),yfit, 10, 'k','filled');
% % ylabel('Predicted Response', 'FontSize', 14);
% % xlabel('Observed Response', 'FontSize', 14);
% corr(y(stats.wasnan==0),yfit, 'Type', 'Pearson')
% subplot(1,2,1);
% [~,~,~] = densityplot(y(stats.wasnan==0),yfit,50,50,30,0,[],'true','square','true');
% box on
% ylabel('Predicted Log Peak Flow', 'FontSize', 14);
% xlabel('Observed Log Peak Flow', 'FontSize', 14);
% set(gca, 'Xlim', [-11 10], 'Ylim', [-11 10], 'Xtick', -11:2:10,'Ytick', -11:2:10, 'FontSize',12)
% hold all
% plot(-12:1:12,-12:1:12, 'k')
% text(1,-8, ['Correlation = ', num2str(round(corr(y(stats.wasnan==0),log(crest_peaks(nooutliers(stats.wasnan==0))), 'Type', 'Pearson')*100)/100)], 'FontSize', 12);
% text(-10,9, 'a)', 'FontSize', 12);

% usmapgeo(usbasinmorph(good_idx,3),usbasinmorph(good_idx,4),uscrestpar(good_idx,2),20,[],[], 'false');
% title('Maximum Water Capacity', 'FontSize', 16)
% ylabel(gco, 'Water Capacity (mm)', 'FontSize', 14)
% set(gco, 'FontSize', 12)
% saveas(gcf, 'Soil_Max_Water_Cap_1k_basins_lt_200mm.fig')

% %% PDFs of Geophysical Factors
% parnames = [{'Wm'}, {'k_s_a_t'}, {'b'}, {'I_R'}, {'k_O'}, {'k_I'}, {'Mean P'}, {'Mean aET'}];
% parUnits = [{'mm'}, {'mm/hr'}, {'unitless'}, {'%'}, {'unitless'}, {'unitless'}, {'mm/hr'}, {'mm/hr'}];
% parScaler = [1, 1, 1, 100, 1, 1, 1, 1];
% 
ko = exp(-1.0582.*log(usbasinlags(good_idx,2).*10.8448)+0.8384);
ko(ko > 1) = 1;
ki = (uscrestpar(good_idx,2).*3/100)./(usbasinmorph(good_idx,6).*usbasinmorph(good_idx,5));
% 
% cont = 1;
% for i = 1:4
%     subplot(2,4,cont);
%     [counts, centers] = hist(uscrestpar(good_idx,i).*parScaler(cont),20);
%     bar(centers,(counts./sum(counts)).*100, 'BarWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
%     set(gca, 'FontSize', 12, 'Xlim', [min(uscrestpar(good_idx,i).*parScaler(cont)) max(uscrestpar(good_idx,i).*parScaler(cont))]);
%     ylabel('Relative Frequency (%)', 'FontSize', 14);
%     xlabel(parUnits{cont}, 'FontSize', 14);
%     title(parnames{cont}, 'FontSize', 14);
%     cont = cont + 1;
% end
% 
% subplot(2,4,cont);
% [counts, centers] = hist(ko,20);
% bar(centers,(counts./sum(counts)).*100, 'BarWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [min(ko) max(ko)]);
% ylabel('Relative Frequency (%)', 'FontSize', 14);
% xlabel(parUnits{cont}, 'FontSize', 14);
% title(parnames{cont}, 'FontSize', 14);
% cont = cont + 1;
% 
% subplot(2,4,cont);
% [counts, centers] = hist(ki,20);
% bar(centers,(counts./sum(counts)).*100, 'BarWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [min(ki) max(ki)]);
% ylabel('Relative Frequency (%)', 'FontSize', 14);
% xlabel(parUnits{cont}, 'FontSize', 14);
% title(parnames{cont}, 'FontSize', 14);
% cont = cont + 1;
% 
% subplot(2,4,cont);
% [counts, centers] = hist(mean_P,20);
% bar(centers,(counts./sum(counts)).*100, 'BarWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [min(mean_P) max(mean_P)]);
% ylabel('Relative Frequency (%)', 'FontSize', 14);
% xlabel(parUnits{cont}, 'FontSize', 14);
% title(parnames{cont}, 'FontSize', 14);
% cont = cont + 1;
% 
% subplot(2,4,cont);
% [counts, centers] = hist(mean_aET,20);
% bar(centers,(counts./sum(counts)).*100, 'BarWidth', 1, 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'k');
% set(gca, 'FontSize', 12, 'Xlim', [min(mean_aET) max(mean_aET)]);
% ylabel('Relative Frequency (%)', 'FontSize', 14);
% xlabel(parUnits{cont}, 'FontSize', 14);
% title(parnames{cont}, 'FontSize', 14);
% cont = cont + 1;