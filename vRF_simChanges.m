% vRF_simChanges.m
%
% addresses issue raised by Sam Schwarzkopf w/ comparing RF parameters
% between conditions
%
% - loads RF params from example subj/ROI
% - generates 100 randomized versions of those parameters based on a few
%   noise models (a) just uniform noise added to x,y,sigma on each iter,
%   (b) noise scales with ecc (linear), (c) uniform noise increase on one
%   side of fixation, but not the other, (d) bias scales with noise
%   (nonlinearly?) - e.g., shift inwards for noisier voxels
% - computes pairwise comparisons between baseline (used to seed model) and
%   each iter; with different binning schemes(?)
% - -- scheme 1: use baseline (seed) to sort, compare different random
%      noise ("correct")
% - -- scheme 2: use one random noise condition to sort, look at changes
%      ("incorrect")
%
% TODO: rewrite so that uniform is a special case of linear slope...
% TODO: ecc > ecc_thresh, scale vectors to that length (so ecc == 15-eps)
% TODO: ecc plots (not just size)
% TODO: (related to above) separate out computing binned values and
%       plotting...
%
% TC Sprague 12/9/2020


subj = 'sub002';
sess = 'barret01';
which_ROI = 'V1';
which_RFs = 'surf';


fn = sprintf('%s/retinotopy/retinotopy_ROIdata/%s_%s_%s_%s.mat',vRF_loadRoot,subj,sess,which_ROI,which_RFs);
fprintf('Loading %s\n',fn);

load(fn,'rf');

% rf structure contains what we want

% limits for voxels we want to consider?
ve_thresh = 0.25; % variance explained >= this
ecc_thresh = 15;  % eccentricity <= this

goodvox = rf.ve >= ve_thresh & rf.ecc <= ecc_thresh;


% n_params (3) x n_vox params drawn from rf structure based on goodvox above
rf_baseline = [rf.x0(goodvox);rf.y0(goodvox);rf.sigma(goodvox)]; 

% number of simulation iterations
n_iter = 50;

noise_mode = 2; % 1 = uniform, 2 = scales w/ ecc

% MODE 1: uniform noise added to x, y, sigma of different magnitudes

if noise_mode == 1
    uni_noise_mag = 5*[0.5 0.5 0.25; 1.5 1.5 0.75]; % add Gaussian noise to each vox for x, y, sigma (in deg vis angle); each row is a different noise level value
    
    
    % cell array: for each noise level (rows above), n_params (3) x n_vox x n_iter - each parameter randomly wiggled
    rf_noisesim = cell(size(uni_noise_mag,1),1);
    for nn = 1:size(uni_noise_mag,1)
        rf_noisesim{nn} = nan(size(rf_baseline,1),size(rf_baseline,2),n_iter);
        for ii = 1:n_iter
            rf_noisesim{nn}(:,:,ii) = rf_baseline + randn(size(rf_baseline)) .* uni_noise_mag(nn,:).';
            rf_noisesim{nn}(3,rf_noisesim{nn}(3,:,ii)<0,ii) = 0.01;
        end
        
    end
    






% MODE 2: noise scales with eccentricity
elseif noise_mode == 2

    % additive gaussian noise w/ std dev = const + slope*ecc
    % Dumoulin & Wandell 2008 - .0625 size vs ecc slope in V1
    lin_noise_const = 2*[0.25 0.25 0.2;    1.5  1.5  0.5];  
    lin_noise_slope = 5*[0.1  0.1  0.0625; 0.25 0.25 0.15];

    % cell array: for each noise level (rows above), n_params (3) x n_vox x n_iter - each parameter randomly wiggled
    rf_noisesim = cell(size(lin_noise_const,1),1);

    tmp_ecc = sqrt(sum(rf_baseline([1 2],:).^2,1));
    
    for nn = 1:length(rf_noisesim)
        
        this_noise_scale = lin_noise_const(nn,:).' + lin_noise_slope(nn,:).' .* tmp_ecc;
        
        rf_noisesim{nn} = nan(size(rf_baseline,1),size(rf_baseline,2),n_iter);
        
        rf_noisesim{nn} = rf_baseline + randn([size(this_noise_scale) n_iter]) .* this_noise_scale;
        tmpsig = rf_noisesim{nn}(3,:,:);
        tmpsig(tmpsig<0) = 0.1;
        rf_noisesim{nn}(3,:,:) = tmpsig; clear tmpsig;
        
        % rescale ecc to be within ecc_thresh (like vRF param restrictions)
        tmpxy  = rf_noisesim{nn}([1 2],:,:);
        tmpecc = sqrt(sum(rf_noisesim{nn}([1 2],:,:).^2,1)); % 1 x nvox x niter
        % those that are >= ecc_thresh, compute a scale factor
        bigecc = tmpecc>ecc_thresh;
        tmpscale = tmpecc(bigecc)/ecc_thresh; % divide by tmpscale
        
        all_scale = ones(size(tmpecc));
        all_scale(bigecc) = 1./tmpscale;
        all_scale = repmat(all_scale,2,1,1);
        
        rf_noisesim{nn}([1 2],:,:) = tmpxy.*all_scale;
        
    end
    




% MODE 3: noise depends on hemisphere [[maybe this is a different scenario
% that deserves a different script...?]]
elseif noise_mode == 3 % TODO








end










% PLOT EVERYTHING

% first - scatterplot of size vs ecc with overlaid bins
% (one subplot per noise level; sorted within simulation) - just the first
% simulation iteration

noise_colors = lines(length(rf_noisesim));


figure;
for nn = 1:length(rf_noisesim)
    subplot(1,length(rf_noisesim),nn); hold on;
    
    this_ecc = sqrt(sum(rf_noisesim{nn}([1 2],:,1).^2,1));
    this_sig = rf_noisesim{nn}(3,:,1);
    
    plot(this_ecc,this_sig,'k.','MarkerSize',5);
    
    
    % bin voxels by this eccentricity
    [ecc_bin_idx,ecc_bin_edges] = discretize(this_ecc,0:1:ecc_thresh);
    bin_centers = mean([ecc_bin_edges(1:end-1);ecc_bin_edges(2:end)],1);
    
    binned_ecc = nan(length(bin_centers),1);
    binned_sig = nan(length(bin_centers),1);
    
    
    for bb = 1:length(bin_centers)
        thisidx = ecc_bin_idx==bb;
        binned_ecc(bb) = mean(this_ecc(thisidx));
        binned_sig(bb) = mean(this_sig(thisidx));
    end
    
    plot(binned_ecc,binned_sig,'o-','Markersize',7,'LineWidth',1.5,'MarkerFaceColor','w','Color',noise_colors(nn,:));
    
    
    ylabel('Size (\circ)'); xlabel('Eccentricity (\circ)');
    set(gca,'TickDir','out','LineWidth',1.5);

    
    hold off;
end


match_ylim(get(gcf,'Children'));

sgtitle('Eccentricity binned within dataset');



% second - comparison of size vs. ecc function differences when using
% baseline bin (left) and binning based on one condition (middle) and
% binning based on each condition individually (right)
% 
% first row - size/ecc plots for each condition (error bars across iter)
% second row - difference (error bars across iter)

rf_baseline_ecc = sqrt(sum(rf_baseline([1 2],:).^2,1));

[baseline_ecc_bin_idx,tmp_edges] = discretize(rf_baseline_ecc,0:1:ecc_thresh);
baseline_ecc_bin_centers = mean([tmp_edges(1:end-1);tmp_edges(2:end)],1);

ecc_bin_centers = baseline_ecc_bin_centers; % we'll just always use these...


figure;


% compare binned size vs ecc, binning using baseline
absax(1) = subplot(2,3,1); hold on;
binned_sig_bybaseline = nan(length(rf_noisesim),length(baseline_ecc_bin_centers),n_iter);
for nn = 1:length(rf_noisesim)
    for bb = 1:length(baseline_ecc_bin_centers)       
        thisidx = baseline_ecc_bin_idx==bb;
        binned_sig_bybaseline(nn,bb,:) = mean(rf_noisesim{nn}(3,thisidx,:),2);
    end
    thism = mean(binned_sig_bybaseline(nn,:,:),3);
    thisci = squeeze(prctile(binned_sig_bybaseline(nn,:,:),[2.5 97.5],3)); % n_bins x 2
    plot(baseline_ecc_bin_centers.*[1;1],thisci.','-','LineWidth',1.5,'Color',noise_colors(nn,:));
    plot(baseline_ecc_bin_centers,thism,'o-','MarkerSize',7,'LineWidth',1.5,'Color',noise_colors(nn,:),'MarkerFaceColor','w');
    clear thism thisci;
end

xlabel('Eccentricity (baseline; \circ)'); ylabel('Size (\circ)');
title('Binned based on baseline eccentricity');
set(gca,'TickDir','out');

% difference (end - 1 in case there are > 2 noise values used)
diffax(1) = subplot(2,3,4); hold on;

thism = mean(binned_sig_bybaseline(end,:,:),3) - mean(binned_sig_bybaseline(1,:,:),3);
thisci = squeeze(prctile(binned_sig_bybaseline(end,:,:)-binned_sig_bybaseline(1,:,:),[2.5 97.5],3)); % n_bins x 2

plot(baseline_ecc_bin_centers.*[1;1],thisci.','-','LineWidth',1.5,'Color','k');
plot(baseline_ecc_bin_centers,thism,'ko-','MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','w');
xlabel('Eccentricity (baseline; \circ)'); ylabel('\Delta size (\circ)');
title('Binned based on baseline eccentricity');
set(gca,'TickDir','out');
clear thisci thism;

% compare binned size vs ecc, binning using one condition (row 1)
absax(2) = subplot(2,3,2); hold on;
binned_sig_byonecond = nan(length(rf_noisesim),length(baseline_ecc_bin_centers),n_iter);

% do the binning per iteration
for ii = 1:n_iter
    
    tmp_ecc = sqrt(sum(rf_noisesim{1}([1 2],:,ii).^2,1)); % TODO - need to recompute on each iter...    
    [this_ecc_bin_idx] = discretize(tmp_ecc,0:1:ecc_thresh);

    for nn = 1:length(rf_noisesim)
        for bb = 1:length(ecc_bin_centers)
            thisidx = this_ecc_bin_idx==bb;
            binned_sig_byonecond(nn,bb,:) = mean(rf_noisesim{nn}(3,thisidx,:),2);
            clear thisidx;
        end
    end
    clear tmp_ecc this_ecc_bin_idx;
end

% plot
for nn = 1:length(rf_noisesim)
    thism  = mean(binned_sig_byonecond(nn,:,:),3);
    thisci = squeeze(prctile(binned_sig_byonecond(nn,:,:),[2.5 97.5],3));
    plot(ecc_bin_centers.*[1;1],thisci.','-','Color',noise_colors(nn,:),'LineWidth',1.5);
    plot(ecc_bin_centers,thism,'o-','MarkerSize',7,'LineWidth',1.5,'Color',noise_colors(nn,:),'MarkerFaceColor','w');
    clear thism thisci;
end
xlabel('Eccentricity (Cond 1; \circ)'); ylabel('Size (\circ)');
title('Binned based on single condition ecc');
set(gca,'TickDir','out');

% difference (end - 1 in case there are > 2 noise values used)
diffax(2) = subplot(2,3,5); hold on;

thism  = mean(binned_sig_byonecond(end,:,:),3) - mean(binned_sig_byonecond(1,:,:),3);
thisci = squeeze(prctile(binned_sig_byonecond(end,:,:) - binned_sig_byonecond(1,:,:),[2.5 97.5],3));

plot(ecc_bin_centers.*[1;1],thisci.','k-','LineWidth',1.5);
plot(ecc_bin_centers,thism,'ko-','MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','w');
clear thism thisci;

xlabel('Eccentricity (Cond 1; \circ)'); ylabel('\Delta size (\circ)');
title('Binned based on single condition ecc');
set(gca,'TickDir','out');

% compare binned size vs ecc, binning using each condition 
absax(3) = subplot(2,3,3); hold on;
binned_sig_byeachcond = nan(length(rf_noisesim),length(baseline_ecc_bin_centers),n_iter);

% do the binning per iteration
for ii = 1:n_iter
    

    for nn = 1:length(rf_noisesim)
        
        tmp_ecc = sqrt(sum(rf_noisesim{nn}([1 2],:,ii).^2,1)); % TODO - need to recompute on each iter...
        [this_ecc_bin_idx] = discretize(tmp_ecc,0:1:ecc_thresh);

        
        for bb = 1:length(ecc_bin_centers)
            thisidx = this_ecc_bin_idx==bb;
            binned_sig_byeachcond(nn,bb,:) = mean(rf_noisesim{nn}(3,thisidx,:),2);
            clear thisidx;
        end
        clear tmp_ecc this_ecc_bin_idx;
    end
    
end

% plot
for nn = 1:length(rf_noisesim)
    thism  = mean(binned_sig_byeachcond(nn,:,:),3);
    thisci = squeeze(prctile(binned_sig_byeachcond(nn,:,:),[2.5 97.5],3));
    plot(ecc_bin_centers.*[1;1],thisci.','-','LineWidth',1.5,'Color',noise_colors(nn,:));
    plot(ecc_bin_centers,thism,'o-','MarkerSize',7,'LineWidth',1.5,'Color',noise_colors(nn,:),'MarkerFaceColor','w');
    clear thism thisci;
end
xlabel('Eccentricity (each condition); \circ)'); ylabel('Size (\circ)');
title('Binned based on each condition ecc');
set(gca,'TickDir','out');

% difference (end - 1 in case there are > 2 noise values used)
diffax(3) = subplot(2,3,6); hold on;
thism = mean(binned_sig_byeachcond(end,:,:),3) - mean(binned_sig_byeachcond(1,:,:),3);
thisci = squeeze(prctile(binned_sig_byeachcond(end,:,:) - binned_sig_byeachcond(1,:,:),[2.5 97.5],3));
plot(ecc_bin_centers.*[1;1],thisci.','k-','LineWidth',1.5);
plot(ecc_bin_centers,thism,'ko-','MarkerSize',7,'LineWidth',1.5,'MarkerFaceColor','w');
xlabel('Eccentricity (Each condition 1; \circ)'); ylabel('\Delta size (\circ)');
title('Binned based on each condition ecc');
set(gca,'TickDir','out');


sgtitle('Size changes, binned by ecc');
match_ylim(absax);
match_ylim(diffax);
