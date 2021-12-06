% function vRF_make_stim_mask.m
%
% Generates stimulus mask based on saved parameters from vRF_stim PTB
% functions for use with voxel/population receptive field analysis code.
% Either can save a file containing images and parameters (in mrVista
% format), if a filename is provided, or can return the stimulus mask
% directly (if an output argument is used). 
%
% NOTE: mrVista assumes i,j coordinates (+y is down), while vRF_stim and 
% vRFtools assume x,y coordinates (+y is up). Thus, the saved file will be
% upside-down, while the output stimulus mask will be right-side up based
% on vRFtools assumptions. 
% 
% INPUT:
% - p:   structure, saved from vRF_stim stimulus presentation scripts
% - fn_out: string, containing the prefix for the filename to save (all but
%           the .mat) - will save out vista-style images/params files
% - TR:  duration , in s, of the sampling rate (TR). If given in ms, will
%        convert to s
% - res: pixel resolution of stimulus mask to create. defaults to fairly
%        high-resolution 270x270. Use a smaller value for faster
%        computation, bigger value for more precision
%
% OUTPUT:
% - stim_mask:  matrix (n_ypix, n_xpix, n_tpts) of stimulus mask on each TR
%               NOTE: +y is UP ('top' rows of the matrix)
%               NOTE: n_ypix == n_xpix right now
% - stim_coords: vector of coordinates corresponding to rows AND columns of
%                stim_mask (at present, square image required)
%
% NOTE: p.start_wait, end_wait should be integer # of TRs... (as should
% step_dur)
%
% TCS, 6/12/2017
% 
% Updates:
% - TCS 12/5/2021: 
%   * added optional output argument for using this as a function
%   * only save files if a filename is defined
%   * updated docs above
%   * output argument is (x,y), saved _images file is (i,j)


function [stim_mask, stim_coords] = vRF_make_stim_mask(p,fn_out,TR,res)


if nargin < 4 || isempty(res)
    res = 270;
end
  
% if TR in ms...
if TR > 100
    TR = TR/1000;
end

dir_all = nan(size(p.bar_pos,1),1);
size_all = nan(size(p.bar_pos,1),1);
idx = 1;

for bb = 1:size(p.seq,1)
    size_all(idx:(idx+p.n_steps-1)) = p.seq(bb,1)*ones(p.n_steps,1);
    dir_all(idx:(idx+p.n_steps-1)) = p.seq(bb,2)*ones(p.n_steps,1);
    idx = idx+p.n_steps;
end

gg = linspace(min(min(p.bar_pos - size_all/2,[],2),[],1),max(max(p.bar_pos + size_all/2,[],2),[],1),res);

TRs_per_step = round(p.step_dur/TR);

n_TRs = round((p.expt_end-p.expt_start)/TR);
fprintf('In total, %i TRs\n',n_TRs);


% ------- PARAMS FILE ------
% need a _params file, which has stimulus.seq and stimulus.seqtiming
% stimulus.seq is 1:n_TRs

stimulus.seq = 1:n_TRs;

% stimulus.seqtiming is TR*(stimulus.seq-1)
stimulus.seqtiming = TR*(stimulus.seq-1);

if ~isempty(fn_out)
    fn_seq = sprintf('%s_params.mat',fn_out);
    save(fn_seq,'stimulus');
end


% ------ IMAGES FILE ------

images = zeros(res,res,n_TRs);


start_TR = 1+round(p.start_wait/TR); % figure out when to start (1+TR/p.start_wait)


% keep track of where we are in iamges
this_TR = start_TR;

for tt = 1:size(p.bar_pos,1)
    
    proto_img = gg>=(sum(p.bar_pos(tt,:)) - size_all(tt)/2) & gg<=(sum(p.bar_pos(tt,:)) + size_all(tt)/2);
    
    % only do verbose output if running to save a file (typically,
    % manually)
    if nargout == 0
        fprintf('Step %i, %i pixels stimulated\n',tt,sum(proto_img));
    end
    
    if dir_all(tt) < 2.5 % up/down - horizontal
        images(:,:,this_TR:(this_TR+TRs_per_step-1)) = repmat(proto_img.',1,res,TRs_per_step);
    else
        images(:,:,this_TR:(this_TR+TRs_per_step-1)) = repmat(proto_img,res,1, TRs_per_step);
    end
    
    this_TR = this_TR+TRs_per_step;
    
end

stim_mask = images;      % stim_mask is in x,y....
stim_coords = gg;
images = flipud(images); % because vista wants i,j coords

if ~isempty(fn_out)
    fn_imgs = sprintf('%s_images.mat',fn_out);

    save(fn_imgs,'images');
end



return
