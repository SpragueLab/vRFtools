% Utility function to plot size vs ecc for best-fit RF models
%
% adapted from an NYU function (not used for publications)
% TCS 11/8/2021

% loop through subj, ROIs, plot size vs ecc


function vRF_plotSizeVsEcc(subj,sess,ROIs,proj_dir)
% TODO: extra input args to specify aspects of the sorting....


VE_thresh = 0.1; % greater than or equal to this

ecc_thresh = 12; % less than or equal to this


if nargin < 1 || isempty(subj)
    subj = {'sub001','sub004','sub005','sub009','sub010'};
end

if nargin < 2 || isempty(sess)
    sess = {'wmApple_pilot01_wmApple','wmApple_pilot01_wmApple','wmApple_pilot01_wmApple','wmApple_pilot01_wmApple','wmApple_pilot01_wmApple'}; % just one sess per subj
end

if nargin < 3 || isempty(ROIs)
    ROIs = {'V1','V2','V3'};    
end

if nargin < 4 || isempty(proj_dir)
    proj_dir = 'wmApple'; % where to look for data
end

n_bins = 10;

% need ROI x bin x subj

all_sig = nan(length(ROIs),n_bins,length(subj));

avg_color = lines(2);
avg_color = avg_color(2,:);

figure; 

for ss = 1:length(subj)
    for vv = 1:length(ROIs)
        
        subplot(length(subj),length(ROIs),(ss-1)*length(ROIs)+vv); hold on;
        
        this_rf = load(sprintf('%s/%s/%s_trialData/%s_%s_%s_surf_trialData.mat',vRF_loadRoot,proj_dir,proj_dir,subj{ss},sess{ss},ROIs{vv}),'rf');
        this_rf = this_rf.rf;
        
        which_vox = this_rf.ve >= VE_thresh & this_rf.ecc <= ecc_thresh;
        
        scatter(this_rf.ecc(which_vox),this_rf.sigma(which_vox),15,[0 0 0],'filled','MarkerFaceAlpha',0.15);
        
        tmp_sigma = this_rf.sigma(which_vox);
        %tmp_sigma = this_rf.sigma(which_vox)./sqrt(this_rf.exp(which_vox));
        
        [~,this_edges,this_bins] = histcounts(this_rf.ecc(which_vox),n_bins,'BinLimits',[0 ecc_thresh]);
        
        if ss == 1 && vv == 1
            bin_centers = mean([this_edges(1:end-1);this_edges(2:end)],1);
        end
        
        bu = unique(this_bins);
        for bb = 1:length(bu)
            all_sig(vv,bu(bb),ss) = mean(tmp_sigma(this_bins==bu(bb)));
        end
        
        plot(bin_centers,all_sig(vv,:,ss),'-','LineWidth',1.5,'Color',avg_color);
        
        if vv == 1
            ylabel(subj{ss});
        end
        
        if ss == 1
            title(ROIs{vv});
        end
        
        clear this_rf;
        
    end
end

%% now plot the avg over subj
figure;

for vv = 1:length(ROIs)
    subplot(1,length(ROIs),vv);
    hold on;
    
    plot(bin_centers,squeeze(nanmean(all_sig(vv,:,:),3)),'ko-','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',3);
    mye = squeeze(nanstd(all_sig(vv,:,:),[],3)) ./ squeeze(sqrt( sum( ~isnan(all_sig(vv,:,:)), 3 ) ));
    for ee = 1:length(mye)
        if ~isinf(mye(ee))
            plot(bin_centers(ee)*[1 1],nanmean(squeeze(all_sig(vv,ee,:)))+[-1 1]*mye(ee),'k-','LineWidth',1.5);
        end
    end
    
    title(ROIs{vv});
    xlim([0 ecc_thresh+1]);
    %set(gca,'XTick',[0 2 12 15]);
    
    if vv == 1
        ylabel('Size (\circ)');
        xlabel('Eccentricity (\circ)');
    else
        set(gca,'YTickLabel',[]);
    end
    
end

match_xlim(get(gcf,'Children')); match_ylim(get(gcf,'Children'));

set(get(gcf,'Children'),'TickDir','out');