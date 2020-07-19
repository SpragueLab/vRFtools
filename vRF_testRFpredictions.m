% vRF_testRFpredictions.m
%
% use barret-based RF models to predict response to hexMap stimuli
%
% - in principle, we should be able to do this for any arbitrary stimulus -
% so in the future we'll just leverage the as-yet-not-existing
% vRF_makePredictions.m function, which takes in an RF model/s and
% 'stimulus' and generates predicted timecourse and/or trial sequence
%
% at present: use best-fit RF data to predict responses to IEM_hexMap
% stimuli, which were acquired as part of fspri
%
% IMPORTANT: for now, ALWAYS use RF params from ret_data - those are using
% gFit, others (hex) are using fFit...
% 
% TCS 7/16/2020


%% load data
ret_root = sprintf('%s/retinotopy/retinotopy_ROIdata',vRF_loadRoot);
hex_root = sprintf('%s/fspri_trialData',fspri_loadRoot);

data_type = 'surf';
ROI = 'V1';
subj = 'sub003';

ret_sess = 'barret01';
hex_sess = 'fspri_pilot01_map';

ret_data_fn = sprintf('%s/%s_%s_%s_%s_25mm.mat',ret_root,subj,ret_sess,ROI,data_type);
hex_data_fn = sprintf('%s/%s_%s_%s_%s_trialData.mat',hex_root,subj,hex_sess,ROI,data_type);

ret_data = load(ret_data_fn);
hex_data = load(hex_data_fn);

VEthresh = 0.5; 

which_vox  = ret_data.rf.ve>=VEthresh;
%which_vox2 = hex_data.rf.ve>=VEthresh;

%if sum(which_vox==which_vox2)~=length(which_vox)
%    error('uh oh...');
%end
clear which_vox2;


%% make RF profile for each voxel

res = 151; % 151 x 151 pix
fov = 24;  % horizontal/vertical field of view

scrpts = linspace(-fov/2,fov/2,res);
[scrx,scry] = meshgrid(scrpts,scrpts);
scrx = reshape(scrx,[numel(scrx),1]); scry = reshape(scry,[numel(scry),1]);

voxidx = find(which_vox);

% n_vox x n_pix (res^2)
% (note: can be vectorized...)
vox_mask = nan(length(voxidx),length(scrx));

for vv = 1:length(voxidx)
    thisr = hypot(scrx-ret_data.rf.x0(voxidx(vv)),scry-ret_data.rf.y0(voxidx(vv)));
    
    % this is where we make the RF image...
    % RF = exp( -.5 * ((Y ./ sigmaMajor).^2 + (X ./ sigmaMinor).^2));
    % (rfGaussian2d.m from vistasoft)
    vox_mask(vv,:) = exp(-0.5 * (thisr./ret_data.rf.sigma(voxidx(vv))).^2  );
end



%% make predictions for each trial of 'challenge' dataset

% TODO: verify pos/neg direction for both fspri and IEM_hexMap

% hex_data.p_all.size == radDeg

% first - make stim masks for each trial of IEM_hexMap
stim_mask = nan(size(hex_data.c_all,1),res.^2);

for tt = 1:size(hex_data.c_all,1)
    stim_mask(tt,:) = hypot(scrx-hex_data.stim_pos_all(tt,1),scry-hex_data.stim_pos_all(tt,2))<=hex_data.p_all.size;
end

% scale stim mask like mrVista does for RF estimation
stim_mask = stim_mask * (fov/res).^2;


% predicted responses: n_trials x n_vox
pred_resp_hex = (stim_mask*vox_mask.').^ret_data.rf.exp(voxidx);


%% process 'challenge' dataset to n_trials x n_vox responses

% TRs to use for evoked stimulus response
stim_range = [7 10]; 

hex_resp = mean(hex_data.dt_allz(:,:,hex_data.which_TRs>=stim_range(1) & hex_data.which_TRs<=stim_range(2)),3);

% limit hex_resp to only 'valid' voxels (those that have VE above our
% thresh above)
hex_resp = hex_resp(:,which_vox);

% limit hex_resp & pred_resp_hex to only trials with no targets
which_trials = hex_data.c_all(:,1)==0;

hex_resp = hex_resp(which_trials,:);
pred_resp_hex = pred_resp_hex(which_trials,:);

hex_pos  = hex_data.stim_pos_all(which_trials,:);


%% evaluate how strong prediction is (corr? regression?)

% correlation (rho)
pred_corr = corr(pred_resp_hex,hex_resp);
pred_corr = diag(pred_corr);

% regression R^2
pred_reg_R2 = nan(size(pred_resp_hex,2),1);

% I think I have to loop over voxels to do this...
for vv = 1:size(pred_resp_hex,2)
    [~,~,~,~,tmpstats] = regress(hex_resp(:,vv),[pred_resp_hex(:,vv) ones(size(pred_resp_hex,1),1)]);
    pred_reg_R2(vv) = tmpstats(1);
    clear tmptstats;
end

% TODO: plot visual field, color vox based on hex_map R2
figure;
subplot(1,2,1); hold on;
scatter(ret_data.rf.x0(voxidx),ret_data.rf.y0(voxidx),25,pred_reg_R2,'Filled','MarkerFaceAlpha',0.5);
colormap viridis;
axis square equal;
title('Predicted R^2 (hexMap)');

subplot(1,2,2); hold on;
scatter(ret_data.rf.x0(voxidx),ret_data.rf.y0(voxidx),25,ret_data.rf.ve(voxidx),'Filled','MarkerFaceAlpha',0.5);
colormap viridis;
axis square equal;
title('Model fit R^2 (ret)');



%% plot RF profile and measured response on each hexMap trial

[~,sorted_vox_idx] = sort(pred_reg_R2,'descend');

offset = sum(isnan(pred_reg_R2)); % because there's a nan because we're including one all-zeros vox...
vox_to_plot = sorted_vox_idx(offset+16:21); % 16, 17, 18, 19, 20, 21 - sub003, V1

circth = linspace(0,2*pi,1001);
circx  = cos(circth);
circy  = sin(circth);

text_pos = [0 10];

scatter_ax = nan(length(vox_to_plot),1);

% plot a bunch of these...
figure;
for vv = 1:length(vox_to_plot)
    
    % plot RF and predicted response for each trial
    subplot(length(vox_to_plot),3,1+(vv-1)*3); hold on;
    plot(ret_data.rf.x0(voxidx(vox_to_plot(vv)))+circx*ret_data.rf.sigma(voxidx(vox_to_plot(vv))),...
         ret_data.rf.y0(voxidx(vox_to_plot(vv)))+circy*ret_data.rf.sigma(voxidx(vox_to_plot(vv))),...
        'k-','LineWidth',1.5);
    
    % here - scatter is most efficient way to show the predicted response for
    % each stim position
    scatter(hex_pos(:,1),hex_pos(:,2),25,pred_resp_hex(:,vox_to_plot(vv)),'Filled','MarkerFaceAlpha',0.5); colormap viridis;
    text(text_pos(1),text_pos(2),'Predicted','HorizontalAlignment','center')
    axis square equal off;
    
    % plot RF and *measured* response on each trial
    subplot(length(vox_to_plot),3,2+(vv-1)*3); hold on;
    plot( ret_data.rf.x0(voxidx(vox_to_plot(vv)))+circx*ret_data.rf.sigma(voxidx(vox_to_plot(vv))),...
          ret_data.rf.y0(voxidx(vox_to_plot(vv)))+circy*ret_data.rf.sigma(voxidx(vox_to_plot(vv))),...
        'k-','LineWidth',1.5);
    
    % here - scatter is most efficient way to show the predicted response for
    % each stim position
    scatter(hex_pos(:,1),hex_pos(:,2),25,hex_resp(:,vox_to_plot(vv)),'Filled','MarkerFaceAlpha',0.5); colormap viridis;
    text(text_pos(1),text_pos(2),'Measured','HorizontalAlignment','center')
    axis square equal off;
    
    % plot scatterplot of measured vs predicted (scaled?)
    scatter_ax(vv) = subplot(length(vox_to_plot),3,vv*3); hold on;
    plot(pred_resp_hex(:,vox_to_plot(vv)),hex_resp(:,vox_to_plot(vv)),'.','MarkerSize',5,'Color',[0 0 0 0.5]);
    % add regression line?
    axis square equal;
    ylabel('Measured');
    xlabel('Predicted');
    
    
end
%match_xlim(scatter_ax); match_ylim(scatter_ax);
set(scatter_ax,'YLim',[-3 3],'XLim',[-0.2 3.5]);


% TODO: try predicted response for each point in time - regression coeff or
% R2 or something for each timepoint in hex_data