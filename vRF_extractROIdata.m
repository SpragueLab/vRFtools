% vRF_extractROIdata.m
% (adapted from fspri_extractROIdata.m)
%
% saves data into ROI files for each subj, session
%
% GOAL: save out files containing retinotopy timecouress for further
% testing alongside a 'neutral' set of RFs (NOT NECESSARILY those fit to
% the data contained within the file!!!!)
%
% IF: sess_rf == sess, then the RF params are those fit to the extracted
% data
% OTHERWISE: saved RF params are based on sess_rf; saved data is based on
% sess
%
%
% stores:
% - best-fit RFs (rf. struct) - from sess_rf (VISTA)
% - data that can be used to fit a model (bar_width_1.nii.gz files...) [z-scored and psc]
%   (from VISTA) --> if sess_rf matches sess, then this is the data that
%   was used to fit the model
% - 'raw' data from each imaging run [z-scored and psc]
% - stimulus sequence [yes, 7/19/2020]
%
% for vRFattn:
% - want 'neutral' RFs, so session files will save data for each condition,
%   and RFs from attending continuous bar
%
%
%
% voxel selection IDENTICAL to fspri script, so those 'mapping' data can be
% loaded to compare
%
% NOTE: we'll probably want to move some of this to a 'project' directory,
% rather than general tools. but for now, we'll just assume that sometimes
% you want to extract all the RF timecourses...will prove useful eventually
%
% NOTE: sFit/fFit computed at UCSB is a bit wonky at present - we'll just
% use gFit for visualization for now
%
% Tommy Sprague, 7/13/2020 tsprague@ucsb.edu

function vRF_extractROIdata(subj,sess,sess_rf,ROIs)

task_dir = 'retinotopy';
root = vRF_loadRoot;

if nargin < 1 || isempty(subj)
    % NOTE: for sub004, need to use _25mm files...

    subj = {'sub003','sub004','sub005'}; % vRFattn subj

end

% this is the sessions to extract data from...
if nargin < 2 || isempty(sess)
    
    % vRFattn: (VSS)
    sess = {{'barret_attnCont_sess01','barret_attnSeq1_sess01','barret_attnSeq2_sess01'},{'barret_attnCont_sess01','barret_attnSeq1_sess01','barret_attnSeq2_sess01','barret_attnCont_sess02','barret_attnSeq1_sess02','barret_attnSeq2_sess02'},{'barret_attnCont_sess01','barret_attnSeq1_sess01','barret_attnSeq2_sess01'}};

end

% heuristic for defining RF session - can make more explicit, but this
% should work ok
if nargin < 3 || isempty(sess_rf)
    
    % by default, use the 'first' session in the list for the 'neutral' RF
    % params (from vista)
    for ss = 1:length(subj)
        sess_rf{ss} = sess{ss}{1};
    end
end

% should work...
root_rf = sprintf('%s/retinotopy/',root);


% in case ROIs live somewhere else...
root_ROI = root_rf;
ROI_dir = 'rois_25mm'; % which ROIs are we using? (vox size...)


if nargin < 4 || isempty(ROIs)
    %ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','TO1','TO2','IPS0','IPS1','IPS2','IPS3'};
    ROIs = {'V1','V2','V3','V3AB','hV4','VO1','VO2','LO1','LO2','IPS0','IPS1','IPS2'}; % wmApple
end

roi_str = cell(size(ROIs));

%task_TRs = 304; % for sub004's 1.3 s TR data, 304 TRs
%task_TRs = 448; % for sub003's earlier version, 448 TRs (& sub005)
%task_TRs = 460; % all other subj (sub006, sub007, sub008, sub010, sub011, sub012) - 460 TRs

% which functional files do we want? func or surf
%func_type = 'surf_25mm'; % 'surf' or 'func' or 'ss5' or 'surf_25mm' (sub004)
func_type = 'surf'; % all but sub004

%func_file = sprintf('%s_volreg_normPctDet*.nii.gz',func_type); % search for this in SUBJ/SESS
func_file = sprintf('r*_%s.nii.gz',func_type); % search for SUBJ_SESS.<this> in SUBJ/SESS

%stim_file = sprintf('%s/%s/',root,task_dir,subj);
%stim_file   = 'bar_stimulus_masks_1300ms_images.mat';
%params_file = 'bar_stimulus_masks_1300ms_params.mat';
stim_file   = 'bar_stimulus_masks_750ms_images.mat';
params_file = 'bar_stimulus_masks_750ms_params.mat';

% which RF file to load/store?           
rf_nii = sprintf('RF_%s-gFit.nii.gz',func_type); % unsmoothed data; smoothed ersion (C2F enabled) used to make ROIs

% where to look in the RF file for params (is there a way to pull this out
% dynamically from nii? doesn't seem stored anywhere...)
RF_paramIdx.phase = 1;
RF_paramIdx.ve    = 2;
RF_paramIdx.ecc   = 3;
RF_paramIdx.sigma = 4;
RF_paramIdx.exp   = 5;
RF_paramIdx.x0    = 6;
RF_paramIdx.y0    = 7;
RF_paramIdx.b     = 8;

rf_fields = fieldnames(RF_paramIdx);


for ss = 1:length(subj)
    
    
    % load all ROIs, best-fit RF data

    rf_file = sprintf('%s%s/%s/%s_%s_vista/%s',root_rf,subj{ss},sess_rf{ss},subj{ss},sess_rf{ss},rf_nii);
    %rf_file = sprintf('%s%s/%s/%s_%s_vista/%s',root_rf,subj{ss},sess_rf{ss},subj{ss},sess_rf{ss},rf_nii);
    fprintf('Loading RF params from %s\n',rf_file);
    rfnii = niftiRead(rf_file);
    
    roi_nii = cell(length(ROIs),1);
    %roi_rf_params = cell(length(ROIs),1); 
    
    roi_rf_struct = cell(length(ROIs),1); 
    for rr = 1:length(ROIs)
        
        
        % is this ROI a cell? if not, make it a cell
        if ~iscell(ROIs{rr})
            this_ROIs = {ROIs{rr}};
        else
            this_ROIs = ROIs{rr};
        end
        
        % for saving files
        roi_str{rr} = horzcat(this_ROIs{:});
        roi_nii{rr} = cell(length(this_ROIs),1);

        % initialize roi_rf_struct{rr} fields
        roi_rf_struct{rr} = struct;
        for ff = 1:length(rf_fields)
            roi_rf_struct{rr}.(rf_fields{ff}) = [];
        end
        roi_rf_struct{rr}.coordIdx = [];
        roi_rf_struct{rr}.ROIidx = []; % which ROI each voxel comes from (from thisROI)
        
        for tt = 1:length(this_ROIs)
            
            
            % updated TCS 9/14/2017 - ROIs will live w/in expt dir, not RF dir,
            % because of different resolution...
            %roifn = sprintf('%s%s/%s/%s_%s_vista/roi/bilat.%s.nii.gz',root_rf,subj{ss},sess_rf{ss},subj{ss},sess_rf{ss},ROIs{rr});
            roifn = sprintf('%s%s/%s/bilat.%s.nii.gz',root_ROI,subj{ss},ROI_dir,this_ROIs{tt});
            fprintf('Loading ROI data from %s\n',roifn);
            roi_nii{rr}{tt} = niftiRead(roifn);
            
            % extract RF params from that ROI
            roi_rf_params = niftiExtract(rfnii,roi_nii{rr}{tt});
            
            for ff = 1:length(rf_fields)
                roi_rf_struct{rr}.(rf_fields{ff}) = horzcat(roi_rf_struct{rr}.(rf_fields{ff}),roi_rf_params(RF_paramIdx.(rf_fields{ff}),:));
            end
            
            
            % for putting back into RAI nii file
            roi_rf_struct{rr}.coordIdx = horzcat(roi_rf_struct{rr}.coordIdx, find(roi_nii{rr}{tt}.data(:)~=0).');

            roi_rf_struct{rr}.ROIidx = horzcat(roi_rf_struct{rr}.ROIidx,ones(1,size(roi_rf_params,2))*tt);

            
            clear roifn;
        end
    end
    
    
    % hmmm....we're goign to want to put everything together
    % eventually...but I guess for now I'll save out different mat files
    % per session; can fix later on
    
    
    for sess_idx = 1:length(sess{ss})
        
        
        % load stimuli & params
        stim_fn   = sprintf('%s/%s/%s/%s/%s_%s_vista/Stimuli/%s',root,task_dir,subj{ss},sess{ss}{sess_idx},subj{ss},sess{ss}{sess_idx},stim_file  );
        params_fn = sprintf('%s/%s/%s/%s/%s_%s_vista/Stimuli/%s',root,task_dir,subj{ss},sess{ss}{sess_idx},subj{ss},sess{ss}{sess_idx},params_file);
        
        stimdata   = load(stim_fn);
        paramsdata = load(params_fn);
        stim_imgs = stimdata.images;
        stim_time = paramsdata.stimulus.seqtiming;
        clear stimdata paramsdata;
        
        % NOTE: not presently tested w/ 1300 ms data...
        % set task_TRs based on subj, session
        if strcmpi(subj{ss},'sub005') || (strcmpi(subj{ss},'sub003') && strcmpi(sess{ss}{sess_idx},'barret_contwidth_pilot01'))
            task_TRs = 448;
        else
            task_TRs = 460; % most subj
        end
        
        
        
        % load 'raw' data
               

        
        % loop through all funcPrefix files within subj, sess and determine
        % which is map and which is task based on num TRs
        this_func_files = dir(sprintf('%s%s/%s/%s/%s_%s.%s',root,task_dir,subj{ss},sess{ss}{sess_idx},subj{ss},sess{ss}{sess_idx},func_file));
        
        sess_nii = cell(length(this_func_files),1);
        for ff = 1:length(this_func_files)
            thisfn = sprintf('%s%s/%s/%s/%s',root,task_dir,subj{ss},sess{ss}{sess_idx},this_func_files(ff).name);
            fprintf('loading %s\n',thisfn);
            sess_nii{ff} = niftiRead(thisfn);
            clear thisfn;
        end
        
        % get # of TRs from each run this session
        nTRs = cellfun(@(x) x.dim(4),sess_nii);
        
        this_task_runs = find(nTRs==task_TRs);
        

        r_all = nan(task_TRs*numel(this_task_runs),2); % run, subrun
        sess_all = sess_idx * ones(task_TRs*numel(this_task_runs),1);
        
        % load average timeseries (the one input to vistasoft)
        % [NOTE: this won't necessarily be identical to averaging
        % individual runs...]
        % d_avg, d_avgz
        
        fn_avg = sprintf('%s%s/%s/%s/%s_%s_vista/bar_seq_1_%s.nii.gz',root,task_dir,subj{ss},sess{ss}{sess_idx},subj{ss},sess{ss}{sess_idx},func_type);
        fprintf('Loading average timecourses from %s\n',fn_avg);
        avg_nii = niftiRead(fn_avg);
        
        
        
        % extract, for each ROI, the time series from each bar mapping run &
        % concatenate
        
        fn_all = cell(length(this_task_runs),1); % empty for now
        
        for rr = 1:length(ROIs)
            
            % since we're indexing into original ROIs, need to turn
            % single-ROIs into cells again (like above)
            if ~iscell(ROIs{rr})
                this_ROIs = {ROIs{rr}};
            else
                this_ROIs = ROIs{rr};
            end
            
            % all runs concatenated
            d_all  = nan(task_TRs*numel(this_task_runs),size(roi_rf_struct{rr}.coordIdx,2));
            d_allz = nan(task_TRs*numel(this_task_runs),size(roi_rf_struct{rr}.coordIdx,2));

            % all runs averaged (loaded from vista version)
            d_avg  = nan(task_TRs,size(roi_rf_struct{rr}.coordIdx,2));
            d_avgz = nan(task_TRs,size(roi_rf_struct{rr}.coordIdx,2));
            
            startidx_TR = 1;
            
            for tt = 1:length(this_task_runs)
                
                fprintf('TASK run %i: %s\n',tt,sess_nii{this_task_runs(tt)}.fname);

                % which TRs we're extracting
                thisidx = startidx_TR:(startidx_TR+task_TRs-1);

                %
                if rr == 1
                    fn_all{tt} = sess_nii{this_task_runs(tt)}.fname; % save the files we loaded from
                end
                
                
                for subroi_idx = 1:length(this_ROIs)
                    
                    mdata = niftiExtract(sess_nii{this_task_runs(tt)},roi_nii{rr}{subroi_idx});
                    
                    
                    d_all( thisidx,roi_rf_struct{rr}.ROIidx==subroi_idx) = mdata; % no manipulation
                    d_allz(thisidx,roi_rf_struct{rr}.ROIidx==subroi_idx) = zscore(mdata,0,1); % Z-SCORE!!!
                    r_all(thisidx) = tt;
                    
                    % only on first pass....
                    if tt == 1
                        mdata = niftiExtract(avg_nii,roi_nii{rr}{subroi_idx});
                        d_avg( :,roi_rf_struct{rr}.ROIidx==subroi_idx) = mdata;
                        d_avgz(:,roi_rf_struct{rr}.ROIidx==subroi_idx) = zscore(mdata,0,1);
                    end
                    
                    clear mdata;
                    
                end
                
                startidx_TR = thisidx(end)+1;
                clear thisidx;
                
            end
            
            % save
            rf = roi_rf_struct{rr};

            this_sess = sess{ss}{sess_idx};
            this_sess_rf = sess_rf{ss};
            
            
            
            fn2s = sprintf('%s%s/%s_ROIdata/%s_%s_%s_%s.mat',root,task_dir,task_dir,subj{ss},sess{ss}{sess_idx},roi_str{rr},func_type);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'d_all','d_allz','d_avg','d_avgz','r_all','sess_all','r_all','fn_all','rf_file','rf','func_file','fn_avg','stim_imgs','stim_time','stim_fn','params_fn','this_sess','this_sess_rf');
            clear rf d_all d_allz d_avg d_avgz this_sess this_sess_rf;
            
        end
        
        % TODO: load best-fit data...(model predictions)
        
    end
    clear sess_nii avg_nii;
end



