% vRF_fitCSS.m
%
% function for fitting CSS vRF model (Kay et al, 2013) to measured
% responses (typically fMRI voxel timecourse), given stimulus sequence
%
% Refactored/reworked version of mrVista's fitting pipeline, attempting to
% remove as much overhead as possible. Operates explicitly on
% matrix/matrices of measured signals, with no reference to brains.
% Optimized to be quick and use GPU, so, in principle, we can do
% resampling, etc. 
%
%
%
% INPUT:
% - fitdata:  data to be fit (n_tpts x n_vox)
% - stimmask: stimulus sequence/mask
%      either: n_pixY x n_pixX x n_tpts OR n_pixY*n_pixX x n_tpts
%      (NOTE: each 'page' should be viewable with imagesc - so first dim,
%      rows, is Y...)
% - stimcoords: coordinates for each pixel of the stimulus array
%      either: n_pixX/Y x 1, if n_pixX==n_pixY, or {n_pixX x 1, n_pixY x 1}
% - samplingrate: time between measurements in fitdata and stimmask, in s
%
% OUTPUT: 
% - bestfit_params: n_params x n_vox matrix of best-fit parameters
% - param_names:    n_params x 1 string array of names of each parameter
% - bestfit_pred:   n_tpts x n_vox matrix of the best-fit model prediction
%                   (timeseries), of identical shape to input data
%
%
%
% TODO:
% - support for non-square stimulus images
%
% TCS 11/9/2021



function [bestfit_params, param_names, bestfit_pred] = vRF_fitCSS(fitdata,stimmask,stimcoords,samplingrate)

% params that govern memory usage, etc
models_per_iter = 5000; % for generating model predictions, how many at a time?


% some 'default' params (TODO)
gridparams_dxy  = 0.5; % 0.25; % for gridfit, how closely to sample X,Y
gridparams_dsig = 0.7; % 0.35; % for gridfit, how closely to sample sigma
gridparams_dexp = 0.15;  % for gridfit, how closely to sample exponent


% check params...
% TODO

% make sure stimcoords is cell array of X, Y stimulus coordinates &
% reformat if necessary
if ~iscell(stimcoords)
    if ismember(2,size(stimcoords))
        % TODO - this is a guess....need to verify...
        stimcoords = {stimcoords(:,1),stimcoords(:,2)};
    elseif ismember(1,size(stimcoords))
        stimcoords = {stimcoords,stimcoords};
    else
        % TODO: make this good...
        error('oops, bad stimcoords')
    end
end



%% generate parameter grid for initial 'coarse' search

% x, y - default to square grid, constrained to a circle with radius equal
% to maximum hypotenuse of stimulus grid (e.g., corners)
maxecc = max(hypot(stimcoords{1},stimcoords{2}));

% (trick here is to make sure we include 0,0...
% tmpgridlim = gridparams_dxy * ceil(maxecc/gridparams_dxy)
n_steps_xy = ceil(maxecc/gridparams_dxy); % how many steps out from 0 (NOT including 0)
gridlims_xy = [-1 1] .* n_steps_xy*gridparams_dxy;

xpts = linspace(gridlims_xy(1),gridlims_xy(2),n_steps_xy*2+1);
ypts = linspace(gridlims_xy(1),gridlims_xy(2),n_steps_xy*2+1);

% sigma - small (~0.1 dva) to very big (~75% screen size)
gridlims_sig = [0.1 (0.75*maxecc)];
sigpts = gridlims_sig(1):gridparams_dsig:gridlims_sig(2);

% exponent - small (0.1) to linear (1)
gridlims_exp = [0.1 1];
exppts = gridlims_exp(1):gridparams_dexp:gridlims_exp(2);

% grid of everythign except x,y (which need to stay together)
[gridx,gridy,gridsig,gridexp] = ndgrid(xpts,ypts,sigpts,exppts);
%gridx = gridx(:);gridy = gridy(:);gridsig = gridsix(:); gridexp = gridexp(:);

gridparams = [gridx(:) gridy(:) gridsig(:) gridexp(:)];

tmpecc = hypot(gridparams(:,1),gridparams(:,2));

% parameter grid - length(gridx) * length(sigpts) * length(exppts) x [x y sig exp]
% (where gridx would be a censored version of x,y points based on our threshold)
gridparams = gridparams(tmpecc<=maxecc,:); clear tmpecc;


%% use parameter grid to generate RF profiles

[scrptsX,scrptsY] = meshgrid(stimcoords{1},stimcoords{2});
scrptsX = scrptsX(:); scrptsY = scrptsY(:); % resahpe to column vectors
scrpts = [scrptsX.'; scrptsY.'];

% NOTE: the below need to be broken up a bit - maybe in groups of 25k
% models?

% reshape stimulus mask into n_pixels x n_tpts, normalize by pixel area
% (resp/deg^2)
stimmask_mat = reshape(stimmask,length(scrptsX),size(stimmask,3));
stimmask_mat = stimmask_mat * (stimcoords{1}(2) - stimcoords{1}(1)) * (stimcoords{2}(2) - stimcoords{2}(1));


% compute HRF convolution kernel, normalize
hrf_tpts = 0:samplingrate:35;
hrf_kernel = rmHrfTwogammas(hrf_tpts); % use default params...
hrf_kernel = hrf_kernel./(sum(hrf_kernel).*samplingrate); % normalize kernel according to volume

%% generate predictions
niter = ceil(size(gridparams,1)/models_per_iter);
startidx = 1;

allpredbold = nan(size(gridparams,1),size(stimmask_mat,2));

for iter = 1:niter

    % pick only a subset of models to work with now
    thisidx = startidx:(startidx+models_per_iter-1);
    thisidx = thisidx(thisidx<=size(gridparams,1)); % cull any dangling rows

    % use gridfit functions to create RF profiles for each grid point
    thisparams = gridparams(thisidx,:);
    thisfilt = vRF_2dGaussian_nonlinear(scrpts,thisparams);

    % generate predicted timeseries for each gridfit point
    thispred = (thisfilt * stimmask_mat) .^ thisparams(:,4);

    % convolve with HRF
    % for now, use rmHrfTwogammas.m, which is what we've used from vista
    % before
    %thispredbold = conv(thispred,hrf_kernel,'same');
    %thispredbold = filter(hrf_kernel(:),1,thispred,[],2); % based on vista
    allpredbold(thisidx,:) = filter(hrf_kernel(:),1,thispred,[],2); % based on vista

    clear thisfilt thisparams thispred;
    startidx = thisidx(end)+1;
end
clear startidx;

%% run gridfit routines








bestfit_params = nan;
param_names = {'x0','y0','sigma','exp','amp','baseline','ve'}; 
bestfit_pred = nan;








return