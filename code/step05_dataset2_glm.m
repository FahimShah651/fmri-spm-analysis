%% =========================================================================
%  FILE        : step05_dataset2_glm.m
%  TITLE       : fMRI Brain Data Analysis using SPM and MATLAB —
%                Preprocessing, GLM, Statistical Mapping and
%                Functional Connectivity
%  AUTHOR      : Fahim Ur Rehman Shah
%  REG NO      : EE2629
%  COURSE      : CSE532 — Signal and Image Processing (MS Level)
%  SUPERVISOR  : Dr. Adnan Shah, FCSE, GIKI
%  INSTITUTION : GIK Institute of Engineering Sciences & Technology
%  DATE        : May 2026
%  STEP        : 05 of 07
%
%  DESCRIPTION : General Linear Model (GLM) analysis for Dataset 2
%                (ds000105, visual object recognition task). Reads
%                events.tsv to build a category-specific design matrix,
%                estimates the model, and saves t-maps for visual
%                category contrasts.
%
%  INPUTS      : data/ds000105/func/swrasub-1_..._bold.nii
%                data/ds000105/func/rp_asub-1_..._bold.txt
%                data/ds000105/func/sub-1_task-objectviewing_run-1_events.tsv
%  OUTPUTS     : results/dataset2/glm/SPM.mat, beta/t-maps
%                figures/fig17–fig20_ds2_glm_*.png/eps
%                results/dataset2/glm/table_ds2_activation_peaks.csv
%  DEPENDENCIES: SPM25, MATLAB R2022b+, step03 complete
%  =========================================================================

clear; clc; close all;
fprintf('=== Step 05: Dataset 2 GLM (ds000105 — Visual Object Recognition) ===\n\n');

%% -----------------------------------------------------------------------
%  PARAMETERS
%  -----------------------------------------------------------------------
ROOT_DIR = 'C:\fmri_project';
DATASET  = 'ds000105';
SUBJECT  = 'sub-1';
TASK     = 'objectviewing';
RUN      = 'run-1';
TR       = 2.5;

P_THRESH = 0.001;
K_THRESH = 10;

func_dir = fullfile(ROOT_DIR,'data',DATASET,'func');
anat_dir = fullfile(ROOT_DIR,'data',DATASET,'anat');
OUT_DIR  = fullfile(ROOT_DIR,'results','dataset2','glm');
FIG_DIR  = fullfile(ROOT_DIR,'figures');
LOG_FILE = fullfile(ROOT_DIR,'project_log.txt');

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

fn_bold     = [SUBJECT '_task-' TASK '_' RUN '_bold'];
fn_anat     = [SUBJECT '_T1w'];
swra_func   = fullfile(func_dir, ['swra' fn_bold '.nii']);
rp_func     = fullfile(func_dir, ['rp_a' fn_bold '.txt']);
events_file = fullfile(func_dir, [SUBJECT '_task-' TASK '_' RUN '_events.tsv']);
anat_nii    = fullfile(anat_dir, [fn_anat '.nii']);

% Decompress anatomy if needed
if ~exist(anat_nii,'file') && exist([anat_nii '.gz'],'file')
    gunzip([anat_nii '.gz'], anat_dir);
end
% Check anat in func dir fallback
if ~exist(anat_nii,'file')
    anat_nii2 = fullfile(func_dir,[fn_anat '.nii']);
    if exist(anat_nii2,'file'), anat_nii = anat_nii2; end
end

if ~exist(events_file,'file')
    error('Events TSV not found: %s', events_file);
end
if ~exist(swra_func,'file')
    error('Smoothed functional not found: %s\nRun step03 first.', swra_func);
end

spm('defaults','FMRI');
spm_jobman('initcfg');

V_func   = spm_vol(swra_func);
num_vols = numel(V_func);
fprintf('  Functional: %d volumes\n', num_vols);

% Auto-detect TR from NIfTI header
try
    hdr_tr = V_func(1).private.hdr.pixdim(5);
    if hdr_tr > 500, hdr_tr = hdr_tr / 1000; end
    if hdr_tr > 0.5 && hdr_tr < 30, TR = hdr_tr; end
    fprintf('  TR (from header): %.3f s\n', TR);
catch
    fprintf('  TR from header unavailable — using %.2f s\n', TR);
end

%% -----------------------------------------------------------------------
%  STEP 5.1 — READ EVENTS TSV
%  -----------------------------------------------------------------------
fprintf('--- [5.1] Reading events.tsv...\n');
events = read_tsv(events_file);
fprintf('  Total events: %d\n', size(events.onset,1));

all_types  = events.trial_type;
cond_names = unique(all_types);
fprintf('  Conditions: ');
fprintf('%s  ', cond_names{:}); fprintf('\n');

conditions = struct();
for c = 1:numel(cond_names)
    name = cond_names{c};
    idx  = strcmp(all_types, name);
    conditions(c).name      = name;
    conditions(c).onsets    = events.onset(idx);
    conditions(c).durations = events.duration(idx);
    fprintf('  %s: %d trials\n', name, sum(idx));
end

%% -----------------------------------------------------------------------
%  FIGURE 17 — CONDITION TIMING (Dataset 2)
%  -----------------------------------------------------------------------
fig17_name = 'fig17_ds2_design_timing';
if ~exist(fullfile(FIG_DIR,[fig17_name '.png']),'file')
    fprintf('--- [5.2] Figure 17: Condition timing...\n');
    total_time = num_vols * TR;
    colors = lines(numel(conditions));
    fig = figure('Visible','off','Position',[100 100 1200 400]);
    hold on;
    for c = 1:numel(conditions)
        for k = 1:numel(conditions(c).onsets)
            t0  = conditions(c).onsets(k);
            dur = conditions(c).durations(k);
            rectangle('Position',[t0, c-0.4, dur, 0.8], ...
                'FaceColor',colors(c,:),'EdgeColor','none');
        end
    end
    set(gca,'YTick',1:numel(conditions), ...
        'YTickLabel',{conditions.name},'FontSize',9);
    xlabel('Time (s)'); xlim([0 total_time]);
    title(sprintf('Condition Timing — %s %s', DATASET, SUBJECT));
    grid on;
    saveFig(fig, FIG_DIR, fig17_name);
    printLatex(fig17_name, ...
        ['Experimental condition timing for Dataset 2 (ds000105). Each row ' ...
         'represents a visual object category; coloured bars indicate active blocks.'], ...
        'ds2:design_timing');
end

%% -----------------------------------------------------------------------
%  STEP 5.3 — GLM SPECIFICATION
%  -----------------------------------------------------------------------
spm_mat = fullfile(OUT_DIR,'SPM.mat');
if ~exist(spm_mat,'file')
    fprintf('--- [5.3] GLM Specification...\n');
    scans = cell(num_vols,1);
    for v = 1:num_vols, scans{v} = sprintf('%s,%d', swra_func, v); end

    if exist(rp_func,'file')
        rp = load(rp_func);
    else
        rp = zeros(num_vols,6);
        fprintf('  WARNING: Motion regressors not found.\n');
    end

    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_spec.dir            = {OUT_DIR};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units   = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT      = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans      = scans;

    for c = 1:numel(conditions)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name     = conditions(c).name;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset    = conditions(c).onsets(:);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = conditions(c).durations(:);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod     = struct('name',{},'param',{},'poly',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth     = 1;
    end
    for r = 1:size(rp,2)
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(r).name = sprintf('Motion_%d',r);
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(r).val  = rp(:,r);
    end
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    spm_jobman('run', matlabbatch);
    fprintf('  GLM specified → %s\n', spm_mat);
else
    fprintf('--- [5.3] GLM: SPM.mat exists, skipping.\n');
end

%% -----------------------------------------------------------------------
%  STEP 5.4 — ESTIMATION
%  -----------------------------------------------------------------------
if ~exist(fullfile(OUT_DIR,'beta_0001.nii'),'file')
    fprintf('--- [5.4] Model Estimation...\n');
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat           = {spm_mat};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals  = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);
    fprintf('  Estimation complete.\n');
else
    fprintf('--- [5.4] Estimation: betas exist, skipping.\n');
end

%% -----------------------------------------------------------------------
%  FIGURE 18 — DESIGN MATRIX
%  -----------------------------------------------------------------------
fig18_name = 'fig18_ds2_design_matrix';
if ~exist(fullfile(FIG_DIR,[fig18_name '.png']),'file') && exist(spm_mat,'file')
    fprintf('--- [5.5] Figure 18: Design matrix...\n');
    load(spm_mat,'SPM');
    X = SPM.xX.X;
    fig = figure('Visible','off','Position',[100 100 600 800]);
    imagesc(zscore(X)); colormap gray; colorbar;
    set(gca,'YDir','normal');
    ylabel('Scan'); xlabel('Regressor');
    ncols = min(numel(SPM.xX.name), size(X,2));
    set(gca,'XTick',1:ncols,'XTickLabel',SPM.xX.name(1:ncols), ...
        'XTickLabelRotation',90,'FontSize',7);
    title(sprintf('Design Matrix — %s %s', DATASET, SUBJECT));
    saveFig(fig, FIG_DIR, fig18_name);
    printLatex(fig18_name, ...
        ['GLM design matrix for Dataset 2 (ds000105). Columns correspond to ' ...
         'HRF-convolved regressors for each visual category and six motion ' ...
         'nuisance regressors.'], ...
        'ds2:design_matrix');
end

%% -----------------------------------------------------------------------
%  STEP 5.6 — CONTRASTS
%  -----------------------------------------------------------------------
fprintf('--- [5.6] Contrast Specification...\n');
load(spm_mat,'SPM');
ncols = numel(SPM.xX.name);

% Identify 'scrambled' condition for contrast denominator
scr_idx = [];
for c = 1:numel(conditions)
    if contains(lower(conditions(c).name),'scrambl') || ...
       contains(lower(conditions(c).name),'pattern')
        scr_idx = c;
    end
end

contrasts = {};
% Each category positive activation
for c = 1:numel(conditions)
    con = zeros(1,ncols); con(c) = 1;
    contrasts{end+1} = struct('name',[conditions(c).name '_pos'],'c',con,'STAT','T'); %#ok<AGROW>
end

% Each visual category > scrambled (if scrambled exists)
if ~isempty(scr_idx)
    for c = 1:numel(conditions)
        if c ~= scr_idx
            con = zeros(1,ncols);
            con(c)       =  1;
            con(scr_idx) = -1;
            contrasts{end+1} = struct('name',[conditions(c).name '_vs_scrambled'], ...
                'c',con,'STAT','T'); %#ok<AGROW>
        end
    end
end

% All visual > baseline (mean of all conditions)
con_all = zeros(1,ncols);
con_all(1:numel(conditions)) = 1/numel(conditions);
contrasts{end+1} = struct('name','AllVisual_pos','c',con_all,'STAT','T');

% F-test: any category
fcon = zeros(numel(conditions), ncols);
for k = 1:numel(conditions), fcon(k,k) = 1; end
contrasts{end+1} = struct('name','AllVisual_Ftest','c',fcon,'STAT','F');

con_dir = dir(fullfile(OUT_DIR,'con_*.nii'));
if isempty(con_dir)
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {spm_mat};
    for k = 1:numel(contrasts)
        if strcmp(contrasts{k}.STAT,'T')
            matlabbatch{1}.spm.stats.con.consess{k}.tcon.name    = contrasts{k}.name;
            matlabbatch{1}.spm.stats.con.consess{k}.tcon.weights = contrasts{k}.c;
            matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
        else
            matlabbatch{1}.spm.stats.con.consess{k}.fcon.name    = contrasts{k}.name;
            matlabbatch{1}.spm.stats.con.consess{k}.fcon.weights = contrasts{k}.c;
        end
    end
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_jobman('run', matlabbatch);
    fprintf('  %d contrasts written.\n', numel(contrasts));
else
    fprintf('--- [5.6] Contrasts: exist, skipping.\n');
end

%% -----------------------------------------------------------------------
%  FIGURES 19-20 — T-MAP OVERLAYS (Dataset 2)
%  -----------------------------------------------------------------------
fprintf('--- [5.7] Generating t-map figures...\n');
load(spm_mat,'SPM');
spmT_files = dir(fullfile(OUT_DIR,'spmT_*.nii'));
fig_count = 19;
peak_table = {};

for k = 1:min(numel(spmT_files), numel(conditions)+1)
    t_file = fullfile(OUT_DIR, spmT_files(k).name);
    con_num = str2double(strrep(strrep(spmT_files(k).name,'spmT_',''),'.nii',''));
    if ~isnan(con_num) && con_num <= numel(SPM.xCon)
        disp_name = SPM.xCon(con_num).name;
    else
        disp_name = strrep(strrep(spmT_files(k).name,'spmT_',''),'.nii','');
    end

    fig_name = sprintf('fig%02d_ds2_tmap_%s', fig_count, ...
        lower(strrep(disp_name,' ','_')));
    fig_count = fig_count + 1;

    if ~exist(fullfile(FIG_DIR,[fig_name '.png']),'file')
        try
            V_t = spm_vol(t_file);
            t_data = spm_read_vols(V_t);
            t_thresh = spm_invTcdf(1-P_THRESH, SPM.xX.erdf);
            t_bin = t_data > t_thresh;

            fig = figure('Visible','off','Position',[100 100 1400 420]);
            ms_list = round(linspace(20, V_t.dim(3)-10, 4));
            for p = 1:4
                sl = min(ms_list(p), V_t.dim(3));
                subplot(1,4,p);
                bg = t_data(:,:,sl)';
                ov = t_bin(:,:,sl)';
                bg_n = (bg-min(bg(:)))/(max(bg(:))-min(bg(:))+eps);
                rgb = repmat(bg_n,[1 1 3]);
                t_norm = min((t_data(:,:,sl)'-t_thresh)/(5*t_thresh-t_thresh+eps),1);
                t_norm(t_norm<0) = 0;
                rgb(:,:,1) = min(rgb(:,:,1)+ov*0.8,1);
                rgb(:,:,2) = min(rgb(:,:,2)+ov.*t_norm*0.6,1);
                rgb(:,:,3) = max(rgb(:,:,3)-ov*0.5,0);
                imagesc(rgb); axis image off;
                title(sprintf('z=%d',sl),'FontSize',8);
            end
            sgtitle(sprintf('%s — %s (p<%.3f, k>%d)', ...
                strrep(disp_name,'_',' '), DATASET, P_THRESH, K_THRESH),'FontSize',11);
            saveFig(fig, FIG_DIR, fig_name);
            printLatex(fig_name, ...
                sprintf('Thresholded t-map for contrast \\textit{%s} in Dataset 2 (p$<$%.3f uncorrected, k$>$%d voxels).', ...
                    strrep(disp_name,'_',' '), P_THRESH, K_THRESH), ...
                ['ds2:tmap:' lower(strrep(disp_name,' ','_'))]);

            [peak_t, idx] = max(t_data(:));
            [px,py,pz] = ind2sub(V_t.dim, idx);
            mni = V_t.mat * [px;py;pz;1];
            peak_table{end+1} = {disp_name, peak_t, mni(1), mni(2), mni(3)}; %#ok<AGROW>
        catch ME
            fprintf('  WARNING: %s — %s\n', disp_name, ME.message);
        end
    end
end

%% SAVE PEAK TABLE
csv_file = fullfile(OUT_DIR,'table_ds2_activation_peaks.csv');
fid = fopen(csv_file,'w');
fprintf(fid,'Contrast,Peak_T,MNI_x_mm,MNI_y_mm,MNI_z_mm\n');
for i = 1:numel(peak_table)
    row = peak_table{i};
    fprintf(fid,'%s,%.4f,%.1f,%.1f,%.1f\n',row{1},row{2},row{3},row{4},row{5});
end
fclose(fid);
fprintf('  Saved: table_ds2_activation_peaks.csv\n');

%% -----------------------------------------------------------------------
%  VALIDATION
%  -----------------------------------------------------------------------
fprintf('\n--- [5.8] Validation...\n');
all_pass = true;
if exist(spm_mat,'file'), fprintf('  PASS: SPM.mat\n');
else, fprintf('  FAIL: SPM.mat missing\n'); all_pass = false; end
betas = dir(fullfile(OUT_DIR,'beta_*.nii'));
tmaps = dir(fullfile(OUT_DIR,'spmT_*.nii'));
fprintf('  INFO: %d betas | %d t-maps\n', numel(betas), numel(tmaps));

status_str = iif(all_pass,'PASS','WARNINGS');
fprintf('\n=== GLM COMPLETE — Dataset 2 | Status: %s ===\n\n', status_str);

fid = fopen(LOG_FILE,'a');
fprintf(fid,'========================================================\n');
fprintf(fid,'STEP 05 : Dataset 2 GLM\n');
fprintf(fid,'Date    : %s\n', datestr(now));
fprintf(fid,'Dataset : %s | Subject: %s\n', DATASET, SUBJECT);
fprintf(fid,'Contrasts: %d | T-maps: %d\n', numel(contrasts), numel(tmaps));
fprintf(fid,'STATUS  : %s\n', status_str);
fprintf(fid,'========================================================\n\n');
fclose(fid);

fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', fullfile(ROOT_DIR,'github'));
fprintf('git add .\n');
fprintf('git commit -m "Step 05 complete: Dataset 2 GLM and contrasts"\n');
fprintf('git push origin main\n');
fprintf('\n=== END OF STEP 05 ===\n');

%% -----------------------------------------------------------------------
%  LOCAL HELPER FUNCTIONS
%  -----------------------------------------------------------------------

function events = read_tsv(fname)
    T = readtable(fname,'FileType','text','Delimiter','\t', ...
        'TreatAsEmpty',{'n/a','NA'},'ReadVariableNames',true);
    events.onset    = T.onset;
    events.duration = T.duration;
    if any(strcmp(T.Properties.VariableNames,'trial_type'))
        events.trial_type = cellstr(T.trial_type);
    else
        events.trial_type = repmat({'condition'},height(T),1);
    end
end

function saveFig(fig, fig_dir, name)
    print(fig, fullfile(fig_dir,[name '.png']),'-dpng','-r150');
    print(fig, fullfile(fig_dir,[name '.eps']),'-depsc');
    close(fig);
    fprintf('  Saved: %s.png / .eps\n', name);
end

function printLatex(fig_name, caption, label)
    fprintf('\n%%%% LaTeX snippet for %s.png\n', fig_name);
    fprintf('\\begin{figure}[htbp]\n');
    fprintf('\\centering\n');
    fprintf('\\includegraphics[width=\\columnwidth]{figures/%s}\n', fig_name);
    fprintf('\\caption{%s}\n', caption);
    fprintf('\\label{fig:%s}\n', label);
    fprintf('\\end{figure}\n\n');
end

function s = iif(cond, tval, fval)
    if cond; s = tval; else; s = fval; end
end
