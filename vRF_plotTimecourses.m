% vRF_plotTimecourses.m
%
% Visualization script for vRF fits 
%
% Plots:
% - time course of selected voxels (see code)
% - best fits [TODO]
% - ...
%
% Things to consider:
% - maybe this is a helper function, and it's called by others which
%   actually load the ROIs?
%
%
%
% Tommy Sprague 7/13/2020 tsprague@ucsb.edu

vRF_file = 'Z:/projects/retinotopy/retinotopy_ROIdata/sub004_barret01_V1_ss5.mat';
thisdata = load(vRF_file);

% stimulus info
stim_file  = 'Z:/projects/retinotopy/sub003/barret01/sub003_barret01_vista/Stimuli/bar_stimulus_masks_1300ms_images.mat';
param_file = 'Z:/projects/retinotopy/sub003/barret01/sub003_barret01_vista/Stimuli/bar_stimulus_masks_1300ms_params.mat';

stim = load(stim_file);
tmp  = load(param_file);
params = tmp.stimulus;

nTRs  = size(stim.images,3);
myTR  = params.seqtiming(2)-params.seqtiming(1); % sec

frames_to_show = [1 9+(0:12)*24 size(stim.images,3)];

nvox = size(thisdata.d_all,2);

%% plot average timecourse for each of n_plot_vox centered left-to-right

ve_thresh = 0.85;

% save a version of RF that's sorted by X0

[~,x0idx] = sort(thisdata.rf.x0,'ascend');
rf_x0sorted = structfun(@(x) x(x0idx),thisdata.rf,'UniformOutput',false);

avg_x0sorted = thisdata.d_avgz(:,x0idx);

thisvox = rf_x0sorted.ve >= ve_thresh;

figure;imagesc(avg_x0sorted(:,thisvox));
figure;plot(rf_x0sorted.x0(thisvox),rf_x0sorted.y0(thisvox),'.')

voffset = .25;
figure; hold on;
thisvoxidx = find(thisvox);
for vv = 1:length(thisvoxidx)
    plot(1:304,avg_x0sorted(:,thisvoxidx(vv))+vv*voffset,'k-');
end
%axis ij;

%% plot average timecourse for each of n_plot_vox (around ring)

ve_thresh = 0.85;

% save a version of RF that's sorted by polar angle

[~,polidx] = sort(thisdata.rf.phase,'ascend');
rf_polsorted = structfun(@(x) x(polidx),thisdata.rf,'UniformOutput',false);

avg_polsorted = thisdata.d_avgz(:,polidx);

thisvox = rf_polsorted.ve >= ve_thresh;

%figure;imagesc(avg_polsorted(:,thisvox));
%figure;plot(rf_polsorted.x0(thisvox),rf_polsorted.y0(thisvox),'.')

circth = linspace(0,2*pi,1001);
circx  = cos(circth);
circy  = sin(circth);

subsample_by = 8;

thisvoxidx = find(thisvox);
thisvoxidx = thisvoxidx(1:subsample_by:end);
thishsv = hsv(1000);

% angle for each hsv value - adjust to match Jerde
thisang = mod(linspace(2*pi/1000,2*pi,1000)+pi,2*pi);

% figure out which thishsv is closest to each polar angle...
thiscolors = nan(length(thisvoxidx),3);
for vv = 1:length(thisvoxidx)
    % min absolute difference b/w
    % thisang-rf_polsorted.phase(thisvoxidx(vv))
    [~,thisidx] = min(abs( thisang - rf_polsorted.phase(thisvoxidx(vv)) ));
    thiscolors(vv,:) = thishsv(thisidx,:);
end

figure;

% vRF positions (sizes?)
subplot(1,3,1); hold on;
for vv = 1:length(thisvoxidx)
    plot((rf_polsorted.sigma(thisvoxidx(vv))*circx)+rf_polsorted.x0(thisvoxidx(vv)),...
         (rf_polsorted.sigma(thisvoxidx(vv))*circy)+rf_polsorted.y0(thisvoxidx(vv)),...
          '-','Color',thiscolors(vv,:));
end
axis equal square;
axis off;
plot(0,0,'k+','MarkerSize',8,'LineWidth',1.5);
plot(12*[-1 -1 1 1 -1],12*[-1 1 1 -1 -1],'k-');
plot([-10 -7],[-10.5 -10.5],'k-','LineWidth',2);


% timeseries
subplot(1,3,[2 3]); hold on;
voffset = .75; % how far apart lines are
for vv = 1:length(thisvoxidx)
    plot(params.seqtiming,avg_polsorted(:,thisvoxidx(vv))+vv*voffset,'-','Color',thiscolors(vv,:));
    if vv == 1
        ymin = min(avg_polsorted(:,thisvoxidx(vv))+vv*voffset);
    elseif vv == length(thisvoxidx)
        ymax = max(avg_polsorted(:,thisvoxidx(vv))+vv*voffset);
    end
end

% show vertical lines for onsets of each bar sweep
yrange = get(gca,'YLim');
plot(params.seqtiming(frames_to_show).*[1;1],[ymin ymax],'-','Color',[0.5 0.5 0.5]);


% mini-ax
mini_ax_loc = [-5 -2]; % origin in x,y (x = sec, y = Z-score)
mini_ax_val = [30 3];  % how long axis is (x = sec, y = Z-score)

plot([0 0 mini_ax_val(1)]+mini_ax_loc(1),[mini_ax_val(2) 0 0]+mini_ax_loc(2),'k-','LineWidth',2.5);
axis off;

set(gcf,'Position',[360         794        1387         544]);


%% plot stim apertures 

figure;

for ii = 1:length(frames_to_show)
    subplot(1,length(frames_to_show),ii); hold on;
    imagesc(flipud(stim.images(:,:,frames_to_show(ii)))); colormap gray; axis square equal tight;
    set(gca,'CLim',[0 1],'XTick',[],'YTick',[],'Box','on');
    %axis off;
    
end