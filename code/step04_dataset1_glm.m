%% =========================================================================
%  FILE        : step04_dataset1_glm.m
%  TITLE       : fMRI Brain Data Analysis using SPM and MATLAB —
%                Preprocessing, GLM, Statistical Mapping and
%                Functional Connectivity
%  AUTHOR      : Fahim Ur Rehman Shah
%  REG NO      : EE2629
%  COURSE      : CSE532 — Signal and Image Processing (MS Level)
%  SUPERVISOR  : Dr. Adnan Shah, FCSE, GIKI
%  INSTITUTION : GIK Institute of Engineering Sciences & Technology
%  DATE        : May 2026
%  STEP        : 04 of 07
%
%  DESCRIPTION : General Linear Model (GLM) analysis for Dataset 1
%                (ds000114, motor task). Reads events.tsv to build the
%                design matrix, estimates the model using SPM HRF, defines
%                motor contrasts, and saves thresholded t-maps.
%
%  INPUTS      : data/ds000114/func/swrasub-01_..._bold.nii
%                data/ds000114/func/rp_asub-01_..._bold.txt
%                data/ds000114/func/sub-01_..._events.tsv
%  OUTPUTS     : results/dataset1/glm/SPM.mat
%                figures/fig11–fig15_ds1_glm_*.png/eps
%                results/dataset1/glm/table_ds1_activation_peaks.csv
%  DEPENDENCIES: SPM25, MATLAB R2022b+, step02 complete
%  =========================================================================

clear; clc; close all;
fprintf('=== Step 04: Dataset 1 GLM (ds000114 — Motor Task) ===\n\n');

%% -----------------------------------------------------------------------
%  PARAMETERS
%  -----------------------------------------------------------------------
ROOT_DIR = 'C:\fmri_project';
DATASET  = 'ds000114';
SUBJECT  = 'sub-01';
SESSION  = 'ses-test';
TASK     = 'fingerfootlips';
RUN      = 'run-1';
TR       = 2.5;            % seconds

% Thresholding
P_THRESH  = 0.001;         % uncorrected p-value
K_THRESH  = 10;            % minimum cluster extent (voxels)

func_dir = fullfile(ROOT_DIR,'data',DATASET,'func');
OUT_DIR  = fullfile(ROOT_DIR,'results','dataset1','glm');
FIG_DIR  = fullfile(ROOT_DIR,'figures');
LOG_FILE = fullfile(ROOT_DIR,'project_log.txt');

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Key filenames
fn_bold     = [SUBJECT '_' SESSION '_task-' TASK '_' RUN '_bold'];
swra_func   = fullfile(func_dir, ['swra' fn_bold '.nii']);
rp_func     = fullfile(func_dir, ['rp_a' fn_bold '.txt']);
events_file = fullfile(func_dir, [SUBJECT '_' SESSION '_task-' TASK '_' RUN '_events.tsv']);

% Check events file (try alternate locations)
if ~exist(events_file,'file')
    events_file = fullfile(func_dir, [SUBJECT '_task-' TASK '_' RUN '_events.tsv']);
end
if ~exist(events_file,'file')
    error('Events TSV not found. Expected: %s', events_file);
end
if ~exist(swra_func,'file')
    error('Smoothed functional not found: %s\nRun step02 first.', swra_func);
end

spm('defaults','FMRI');
spm_jobman('initcfg');

V_func   = spm_vol(swra_func);
num_vols = numel(V_func);
fprintf('  Functional: %d volumes\n', num_vols);

% Auto-detect TR from NIfTI header (overrides hardcoded value)
try
    hdr_tr = V_func(1).private.hdr.pixdim(5);
    if hdr_tr > 500, hdr_tr = hdr_tr / 1000; end
    if hdr_tr > 0.5 && hdr_tr < 30, TR = hdr_tr; end
    fprintf('  TR (from header): %.3f s\n', TR);
catch
    fprintf('  TR from header unavailable — using %.2f s\n', TR);
end

%% -----------------------------------------------------------------------
%  STEP 4.1 — READ EVENTS TSV
%  -----------------------------------------------------------------------
fprintf('--- [4.1] Reading events.tsv...\n');
events = read_tsv(events_file);
fprintf('  Events file: %s\n', events_file);
fprintf('  Total events: %d\n', size(events,1));

% Extract unique conditions (excluding 'instructions' / 'rest' if present)
all_types   = events.trial_type;
cond_names  = unique(all_types);
fprintf('  Conditions found: ');
fprintf('%s  ', cond_names{:}); fprintf('\n');

% Build per-condition onset/duration vectors (in seconds)
conditions = struct();
for c = 1:numel(cond_names)
    name = cond_names{c};
    idx  = strcmp(all_types, name);
    conditions(c).name     = name;
    conditions(c).onsets   = events.onset(idx);
    conditions(c).durations = events.duration(idx);
    fprintf('  %s: %d trials, mean duration=%.1fs\n', ...
        name, sum(idx), mean(events.duration(idx)));
end

%% -----------------------------------------------------------------------
%  STEP 4.2 — FIGURE 11: DESIGN MATRIX PREVIEW
%  (Shows condition timing before SPM estimation)
%  -----------------------------------------------------------------------
fig11_name = 'fig11_ds1_design_timing';
if ~exist(fullfile(FIG_DIR,[fig11_name '.png']),'file')
    fprintf('--- [4.2] Figure 11: Condition timing...\n');
    scan_time = (0:num_vols-1)*TR;
    total_time = num_vols * TR;
    colors = lines(numel(conditions));
    fig = figure('Visible','off','Position',[100 100 1200 350]);
    hold on;
    for c = 1:numel(conditions)
        for k = 1:numel(conditions(c).onsets)
            t0 = conditions(c).onsets(k);
            dur = conditions(c).durations(k);
            rectangle('Position',[t0, c-0.4, dur, 0.8], ...
                'FaceColor',colors(c,:),'EdgeColor','none');
        end
    end
    set(gca,'YTick',1:numel(conditions), ...
        'YTickLabel',{conditions.name},'FontSize',10);
    xlabel('Time (s)'); xlim([0 total_time]);
    title(sprintf('Condition Timing — %s %s', DATASET, SUBJECT));
    grid on;
    saveFig(fig, FIG_DIR, fig11_name);
    printLatex(fig11_name, ...
        ['Experimental condition timing for Dataset 1 (ds000114). Each row ' ...
         'shows one condition; coloured bars indicate active blocks.'], ...
        'ds1:design_timing');
end

%% -----------------------------------------------------------------------
%  STEP 4.3 — SPM GLM SPECIFICATION
%  -----------------------------------------------------------------------
spm_mat = fullfile(OUT_DIR, 'SPM.mat');
if ~exist(spm_mat,'file')
    fprintf('--- [4.3] GLM Specification...\n');

    % Build scan list
    scans = cell(num_vols,1);
    for v = 1:num_vols
        scans{v} = sprintf('%s,%d', swra_func, v);
    end

    % Load realignment parameters (motion regressors)
    if exist(rp_func,'file')
        rp = load(rp_func);
        fprintf('  Adding %d motion regressors\n', size(rp,2));
    else
        rp = zeros(num_vols,6);
        fprintf('  WARNING: Motion regressors not found.\n');
    end

    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_spec.dir          = {OUT_DIR};
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT    = TR;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t  = 16;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

    matlabbatch{1}.spm.stats.fmri_spec.sess.scans  = scans;

    % Conditions
    for c = 1:numel(conditions)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).name     = conditions(c).name;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).onset    = conditions(c).onsets(:);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).duration = conditions(c).durations(:);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).tmod     = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).pmod     = struct('name',{},'param',{},'poly',{});
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(c).orth     = 1;
    end

    % Motion regressors as nuisance covariates
    for r = 1:size(rp,2)
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(r).name = ...
            sprintf('Motion_%d', r);
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress(r).val  = rp(:,r);
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf        = 128;
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';

    spm_jobman('run', matlabbatch);
    fprintf('  GLM specified → %s\n', spm_mat);
else
    fprintf('--- [4.3] GLM Specification: SPM.mat exists, skipping.\n');
end

%% -----------------------------------------------------------------------
%  STEP 4.4 — MODEL ESTIMATION
%  -----------------------------------------------------------------------
spm_mat_est = fullfile(OUT_DIR, 'beta_0001.nii');
if ~exist(spm_mat_est,'file')
    fprintf('--- [4.4] Model Estimation...\n');
    clear matlabbatch;
    matlabbatch{1}.spm.stats.fmri_est.spmmat          = {spm_mat};
    matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
    spm_jobman('run', matlabbatch);
    fprintf('  Estimation complete.\n');
else
    fprintf('--- [4.4] Model Estimation: betas exist, skipping.\n');
end

%% -----------------------------------------------------------------------
%  FIGURE 12 — DESIGN MATRIX (SPM-generated)
%  -----------------------------------------------------------------------
fig12_name = 'fig12_ds1_design_matrix';
if ~exist(fullfile(FIG_DIR,[fig12_name '.png']),'file') && exist(spm_mat,'file')
    fprintf('--- [4.5] Figure 12: Design matrix...\n');
    load(spm_mat, 'SPM');
    X = SPM.xX.X;
    fig = figure('Visible','off','Position',[100 100 600 800]);
    imagesc(zscore(X));
    colormap gray; colorbar;
    set(gca,'YDir','normal');
    ylabel('Scan'); xlabel('Regressor');
    % Label columns
    ncols = min(numel(SPM.xX.name), size(X,2));
    set(gca,'XTick',1:ncols,'XTickLabel',SPM.xX.name(1:ncols), ...
        'XTickLabelRotation',90,'FontSize',7);
    title(sprintf('Design Matrix — %s %s', DATASET, SUBJECT));
    saveFig(fig, FIG_DIR, fig12_name);
    printLatex(fig12_name, ...
        ['GLM design matrix for Dataset 1 (ds000114). Columns show the HRF-convolved ' ...
         'regressors for each motor condition plus six motion-parameter nuisance regressors. ' ...
         'Each row corresponds to one acquired volume.'], ...
        'ds1:design_matrix');
end

%% -----------------------------------------------------------------------
%  STEP 4.6 — CONTRAST SPECIFICATION AND T-MAPS
%  -----------------------------------------------------------------------
fprintf('--- [4.6] Contrast Specification...\n');
load(spm_mat, 'SPM');
cond_labels = {SPM.xX.name{1:numel(conditions)}};  % first N columns = conditions

% Identify condition indices by name
get_idx = @(nm) find(strcmpi(cond_labels, nm));

% Build contrast vectors (length = total regressors)
ncols = numel(SPM.xX.name);

% Motor contrasts
contrasts = {};
motor_names = {'foot','finger','lips','lips_tongue','fingerfootlips','all_motor'};
present_motor = {};
for c = 1:numel(conditions)
    nm = lower(conditions(c).name);
    if ~contains(nm,'instruct') && ~contains(nm,'rest') && ~contains(nm,'fix')
        present_motor{end+1} = conditions(c).name;  %#ok<AGROW>
    end
end

% Positive activation for each motor condition
for c = 1:numel(conditions)
    nm = conditions(c).name;
    idx_c = c;  % conditions are in order in design matrix
    con = zeros(1, ncols);
    con(idx_c) = 1;
    contrasts{end+1} = struct('name', [nm '_pos'], 'c', con, 'STAT','T'); %#ok<AGROW>
end

% All motor conditions combined (exclude 'instructions')
motor_idx = [];
for c = 1:numel(conditions)
    nm = lower(conditions(c).name);
    if ~contains(nm,'instruct') && ~contains(nm,'rest') && ~contains(nm,'fix')
        motor_idx(end+1) = c; %#ok<AGROW>
    end
end
if ~isempty(motor_idx)
    con = zeros(1,ncols);
    con(motor_idx) = 1/numel(motor_idx);
    contrasts{end+1} = struct('name','AllMotor_pos','c',con,'STAT','T');
end

% F-contrast: any motor > baseline
if ~isempty(motor_idx)
    fcon = zeros(numel(motor_idx), ncols);
    for k = 1:numel(motor_idx), fcon(k, motor_idx(k)) = 1; end
    contrasts{end+1} = struct('name','AllMotor_Ftest','c',fcon,'STAT','F');
end

% Write contrasts
con_dir = dir(fullfile(OUT_DIR,'con_*.nii'));
if isempty(con_dir)
    clear matlabbatch;
    matlabbatch{1}.spm.stats.con.spmmat = {spm_mat};
    for k = 1:numel(contrasts)
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.name    = contrasts{k}.name;
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.weights = contrasts{k}.c;
        matlabbatch{1}.spm.stats.con.consess{k}.tcon.sessrep = 'none';
        if strcmp(contrasts{k}.STAT,'F')
            matlabbatch{1}.spm.stats.con.consess{k} = rmfield(...
                matlabbatch{1}.spm.stats.con.consess{k},'tcon');
            matlabbatch{1}.spm.stats.con.consess{k}.fcon.name    = contrasts{k}.name;
            matlabbatch{1}.spm.stats.con.consess{k}.fcon.weights = contrasts{k}.c;
        end
    end
    matlabbatch{1}.spm.stats.con.delete = 0;
    spm_jobman('run', matlabbatch);
    fprintf('  %d contrasts written.\n', numel(contrasts));
else
    fprintf('--- [4.6] Contrasts: already exist, skipping.\n');
end

%% -----------------------------------------------------------------------
%  FIGURES 13-15 — T-MAP OVERLAYS ON ANATOMY
%  -----------------------------------------------------------------------
fprintf('--- [4.7] Generating activation map figures...\n');
load(spm_mat,'SPM');

% Get T-map files
spmT_files = dir(fullfile(OUT_DIR,'spmT_*.nii'));

% Load anatomy for overlay
fn_anat = [SUBJECT '_' SESSION '_T1w'];
anat_nii = fullfile(ROOT_DIR,'data',DATASET,'anat',[fn_anat '.nii']);
if ~exist(anat_nii,'file')
    anat_gz = [anat_nii '.gz'];
    if exist(anat_gz,'file'), gunzip(anat_gz, fileparts(anat_nii)); end
end

fig_count = 13;
peak_table = {};

for k = 1:min(numel(spmT_files), numel(conditions)+1)
    t_file = fullfile(OUT_DIR, spmT_files(k).name);
    con_name = strrep(spmT_files(k).name, 'spmT_','');
    con_name = strrep(con_name, '.nii','');

    % Get contrast name from SPM
    con_num = str2double(con_name);
    if ~isnan(con_num) && con_num <= numel(SPM.xCon)
        disp_name = SPM.xCon(con_num).name;
    else
        disp_name = con_name;
    end

    fig_name = sprintf('fig%02d_ds1_tmap_%s', fig_count, ...
        lower(strrep(disp_name,' ','_')));
    fig_count = fig_count + 1;

    if ~exist(fullfile(FIG_DIR,[fig_name '.png']),'file')
        fprintf('  Generating: %s\n', fig_name);
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
                % Normalise background
                bg_n = (bg - min(bg(:)))/(max(bg(:))-min(bg(:))+eps);
                rgb = repmat(bg_n,[1 1 3]);
                % Overlay activations in red-yellow
                t_norm = min((t_data(:,:,sl)' - t_thresh) / (5*t_thresh - t_thresh + eps), 1);
                t_norm(t_norm<0) = 0;
                rgb(:,:,1) = min(rgb(:,:,1) + ov * 0.8, 1);
                rgb(:,:,2) = min(rgb(:,:,2) + ov .* t_norm * 0.6, 1);
                rgb(:,:,3) = max(rgb(:,:,3) - ov * 0.5, 0);
                imagesc(rgb); axis image off;
                title(sprintf('z=%d', sl),'FontSize',8);
            end
            sgtitle(sprintf('%s — %s (p<%.3f, k>%d)', ...
                disp_name, DATASET, P_THRESH, K_THRESH),'FontSize',11);
            saveFig(fig, FIG_DIR, fig_name);
            printLatex(fig_name, ...
                sprintf('Thresholded t-map for contrast \\textit{%s} in Dataset 1 (p$<$%.3f uncorrected, k$>$%d voxels). Activations shown in red-yellow overlay on the mean EPI background.', ...
                    strrep(disp_name,'_',' '), P_THRESH, K_THRESH), ...
                ['ds1:tmap:' lower(strrep(disp_name,' ','_'))]);

            % Collect peak voxel info
            [peak_t, idx] = max(t_data(:));
            [px,py,pz] = ind2sub(V_t.dim, idx);
            mni = V_t.mat * [px;py;pz;1];
            peak_table{end+1} = {disp_name, peak_t, mni(1), mni(2), mni(3)}; %#ok<AGROW>
        catch ME
            fprintf('  WARNING: Could not plot %s: %s\n', disp_name, ME.message);
        end
    end
end

%% -----------------------------------------------------------------------
%  SAVE ACTIVATION PEAK TABLE (CSV)
%  -----------------------------------------------------------------------
fprintf('--- [4.8] Saving peak activation table...\n');
csv_file = fullfile(OUT_DIR,'table_ds1_activation_peaks.csv');
fid = fopen(csv_file,'w');
fprintf(fid,'Contrast,Peak_T,MNI_x_mm,MNI_y_mm,MNI_z_mm\n');
for i = 1:numel(peak_table)
    row = peak_table{i};
    fprintf(fid,'%s,%.4f,%.1f,%.1f,%.1f\n', ...
        row{1}, row{2}, row{3}, row{4}, row{5});
end
fclose(fid);
fprintf('  Saved: table_ds1_activation_peaks.csv\n');

%% -----------------------------------------------------------------------
%  VALIDATION
%  -----------------------------------------------------------------------
fprintf('\n--- [4.9] Validation...\n');
all_pass = true;
if exist(spm_mat,'file')
    fprintf('  PASS: SPM.mat exists\n');
else
    fprintf('  FAIL: SPM.mat missing\n'); all_pass = false;
end
betas = dir(fullfile(OUT_DIR,'beta_*.nii'));
fprintf('  INFO: %d beta images found\n', numel(betas));
tmaps = dir(fullfile(OUT_DIR,'spmT_*.nii'));
fprintf('  INFO: %d t-maps found\n', numel(tmaps));
cons  = dir(fullfile(OUT_DIR,'con_*.nii'));
fprintf('  INFO: %d contrast images found\n', numel(cons));

status_str = iif(all_pass,'PASS','WARNINGS');
fprintf('\n=== GLM COMPLETE — Dataset 1 | Status: %s ===\n\n', status_str);

%% PROJECT LOG
fid = fopen(LOG_FILE,'a');
fprintf(fid,'========================================================\n');
fprintf(fid,'STEP 04 : Dataset 1 GLM\n');
fprintf(fid,'Date    : %s\n', datestr(now));
fprintf(fid,'Dataset : %s | Subject: %s\n', DATASET, SUBJECT);
fprintf(fid,'Contrasts: %d | T-maps: %d\n', numel(contrasts), numel(tmaps));
fprintf(fid,'Threshold: p<%.4f uncorr, k>%d vox\n', P_THRESH, K_THRESH);
fprintf(fid,'STATUS  : %s\n', status_str);
fprintf(fid,'========================================================\n\n');
fclose(fid);

fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', fullfile(ROOT_DIR,'github'));
fprintf('git add .\n');
fprintf('git commit -m "Step 04 complete: Dataset 1 GLM and contrasts"\n');
fprintf('git push origin main\n');
fprintf('\n=== END OF STEP 04 ===\n');

%% -----------------------------------------------------------------------
%  LOCAL HELPER FUNCTIONS
%  -----------------------------------------------------------------------

function events = read_tsv(fname)
    % Read BIDS events.tsv into a struct with fields onset, duration, trial_type
    fid = fopen(fname,'r');
    header_line = fgetl(fid);
    fclose(fid);
    cols = strsplit(strtrim(header_line), '\t');

    % Use readtable for robust TSV parsing
    T = readtable(fname, 'FileType','text', 'Delimiter','\t', ...
        'TreatAsEmpty',{'n/a','NA'}, 'ReadVariableNames',true);

    events.onset      = T.onset;
    events.duration   = T.duration;
    if any(strcmp(T.Properties.VariableNames,'trial_type'))
        events.trial_type = cellstr(T.trial_type);
    else
        events.trial_type = repmat({'condition'}, height(T), 1);
    end
end

function saveFig(fig, fig_dir, name)
    print(fig, fullfile(fig_dir,[name '.png']), '-dpng','-r150');
    print(fig, fullfile(fig_dir,[name '.eps']), '-depsc');
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
