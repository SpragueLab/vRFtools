% run_vRFfits.m
%
% example script for setting up vRF fits using vRF_fitCSS.m
%
% TCS 11/9/2021

show_plots = 0;

%% EXAMPLE 1: model fits on extracted ROI data
% first - we'll assume we've extracted ROI data, for simplicity (and
% testing)

% this is a file extracted using vRF_extractROIdata.m
ROI_fn = sprintf('%s/retinotopy/retinotopy_ROIdata/sub001_barret_contwidth01_V1_ss5.mat',vRF_loadRoot);
roi_data = load(ROI_fn);

% the above loads a file generated using vRF_extractROIdata.m, which gives
% us stacked 'raw' data from this ROI (each run stacked), and averaged data
% across all runs. also, there's a field that includes the stimulus images.
% These can be created from a vRF_stim session
% (https://github.com/SpragueLab/vRF_stim) using the vRF_make_stim_mask.m
% script in vRFtools

% we're interested in the following variables (stored as fields of roi_data):
% - stim_imgs - n_pixY x n_pixX x n_TRs - stimulus mask
% - stim_time - n_TRs (used for extracting sampling rate)
% - d_allz - n_TRs * n_runs x n_vox - z-scored data from each run
% - r_all  - n_TRs x 1 - run label for each timepoint

% average data across runs
% (NOTE: this is different from what we do in shell scripts - in those, we
% average the 'raw' data...)
% TODO: include a flag for whether to use the 'original' version, or
% whether to do this...

ru = unique(roi_data.r_all(:,1));
tmpd = nan(size(roi_data.d_allz,1)/length(ru),size(roi_data.d_allz,2),length(ru)); % 'pages' of each run
for rr = 1:length(ru)
    tmpd(:,:,rr) = roi_data.d_allz(roi_data.r_all(:,1)==ru(rr),:);
end

fitdata = mean(tmpd,3); % data to fit!

% take out the stimulus mask & coordinates
% version 1: use Y x X x tpts
stimmask = roi_data.stim_imgs;
maxdist = 8.75;
stimcoords = linspace(-maxdist,maxdist,size(stimmask,1)); % this i
% TODO: load the above....

samplingrate = roi_data.stim_time(2)-roi_data.stim_time(1);

[thisparams,thisnames,thispred] = vRF_fitCSS(fitdata,stimmask,stimcoords,samplingrate);

% plot a few example voxels
if show_plots == 1
    startvox = 151; nrows = 3; ncols = 5; voxidx = startvox:(startvox+nrows*ncols-1);
    figure;
    for ii = 1:length(voxidx)
        subplot(nrows,ncols,ii); hold on;
        plot(samplingrate*(1:size(fitdata,1)), fitdata (:,voxidx(ii)) );
        plot(samplingrate*(1:size(fitdata,1)), thispred(:,voxidx(ii)) );
        title(sprintf('Voxel %i',voxidx(ii)));
    end
end
% version 2: use (Y*X) x tpts
% TODO


% second - we'll generate some simulated data, add noise, and run the fits





% third - we'll run on the whole brain....(?)





% fourth - we'll show how to do resampling