% test_vRFfits.m
%
% example script demonstrating how to set up vRF fits using vRF_fitCSS.m
%
% TCS 11/9/2021

show_plots = 0;

%% EXAMPLE 1: model fits on extracted ROI data
% first - we'll assume we've extracted ROI data, for simplicity (and
% testing)

proj = 'retinotopy';
subj = 'sub001';
sess = 'barret_contwidth01';
ROI  = 'V1';
data_type = 'ss5';

% this is a file extracted using vRF_extractROIdata.m
ROI_fn = sprintf('%s/%s/%s_ROIdata/%s_%s_%s_%s.mat',vRF_loadRoot,proj,proj,subj,sess,ROI,data_type);
%ROI_fn = sprintf('%s/retinotopy/retinotopy_ROIdata/sub001_barret_contwidth01_V1_ss5.mat',vRF_loadRoot);
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

mask_mode = 2;

if mask_mode == 1
    % take out the stimulus mask & coordinates
    % version 1: use Y x X x tpts
    stimmask = flipud(roi_data.stim_imgs); % the mask saved here is in vista coords
    maxdist = 8.7462;
    stimcoords = linspace(-maxdist,maxdist,size(stimmask,1)); % this i
    samplingrate = roi_data.stim_time(2)-roi_data.stim_time(1);

elseif mask_mode == 2
    % use behavioral file for this subj to generate the stimulus mask &
    % coordinates
    behav_fs = sprintf('%s/%s/rawbehav/%s/%s/%s_r01_RF*.mat',vRF_loadRoot,proj,subj,sess,subj);
    tmpf = dir(behav_fs);
    behav_fn = sprintf('%s/%s',tmpf.folder,tmpf.name);
    thisbehav = load(behav_fn);
    thisbehav = thisbehav.p; % only p is saved, so let's just call thisbehav.p thisbehav

    [stimmask,stimcoords] = vRF_make_stim_mask(thisbehav,[],thisbehav.TR,151);

    samplingrate = thisbehav.TR;
end


[thisparams,thisnames,thispred,thisparams_grid,thispred_grid] = vRF_fitCSS(fitdata,stimmask,stimcoords,samplingrate);

fn2s = sprintf('%s/retinotopy/vRF_testfits/%s_%s_%s_%s_testfits_mask%i_res%i.mat',vRF_loadRoot,subj,sess,ROI,data_type,mask_mode,size(stimmask,1));
fprintf('saving to %s\n',fn2s);
save(fn2s,'thisparams','thisnames','thispred','thisparams_grid','thispred_grid','mask_mode','stimmask','fitdata','ROI_fn','stimcoords');

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