% vRF_plotModelPredictionPerformance.m
%
% Plots, across subj/ROIs, average model performance (R^2) for predicting
% held-out data using vRF parameters
%
% Calls vRF_testRFpredictions.m for each ROI, subj
%
% TCS 7/18/2020

function vRF_plotModelPredictionPerformance(subj,ret_sess,hex_sess,ROIs,VEthresh)

data_type = 'surf'; % unlikely we'll change this across runs...

if nargin < 1 || isempty(subj)
    subj = {'sub002','sub003','sub004'};
end

if nargin < 2 || isempty(ret_sess)
    ret_sess = {'barret01','barret01','barret01'};
end

if nargin < 3 || isempty(hex_sess)
    hex_sess = {'fspri_pilot01_map','fspri_pilot01_map','fspri_pilot01_map'};
end

if nargin < 4 || isempty(ROIs)
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2'};
end

if nargin < 5 || isempty(VEthresh)
    VEthresh = 0.3;
end

all_R2   = cell(length(subj),length(ROIs));
all_corr = cell(length(subj),length(ROIs));

for ss = 1:length(subj)
    for vv = 1:length(ROIs)
        
        [all_R2{ss,vv},all_corr{ss,vv}] = vRF_testRFpredictions(subj{ss},ret_sess{ss},hex_sess{ss},data_type,ROIs{vv},VEthresh);
        
    end
end

%% plot average R2 & corr across all voxels for each subj
% figure;
% plot(cellfun(@mean,all_R2).','o','MarkerSize',7,'LineWidth',1.5);
% set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'TickDir','out','Box','off');

figure; hold on;
fisherr2z = @(r) atanh(r);
all_corrz = cellfun(@(r) fisherr2z(r),all_corr,'UniformOutput',false);
plot(tanh(cellfun(@mean,all_corrz)).','o','MarkerSize',7,'LineWidth',1.5);
legend(subj);
set(gca,'XTick',1:length(ROIs),'XTickLabel',ROIs,'TickDir','out','Box','off','XTickLabelRotation',-45,'FontSize',12);
xlabel('Region of interest');
ylabel('Model prediction accuracy (\rho)');
xlim([0 length(ROIs)+1]);
ylim([-0.1 1]);

return