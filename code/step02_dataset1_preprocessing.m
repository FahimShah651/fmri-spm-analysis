%% =========================================================================
%  FILE        : step02_dataset1_preprocessing.m
%  TITLE       : fMRI Brain Data Analysis using SPM and MATLAB —
%                Preprocessing, GLM, Statistical Mapping and
%                Functional Connectivity
%  AUTHOR      : Fahim Ur Rehman Shah
%  REG NO      : EE2629
%  COURSE      : CSE532 — Signal and Image Processing (MS Level)
%  SUPERVISOR  : Dr. Adnan Shah, FCSE, GIKI
%  INSTITUTION : GIK Institute of Engineering Sciences & Technology
%  DATE        : May 2026
%  STEP        : 02 of 07
%
%  DESCRIPTION : Full SPM preprocessing pipeline for Dataset 1 (ds000114,
%                Motor/finger-foot-lips task). Steps executed in order:
%                Slice Timing → Realign (Estimate+Reslice) → Coregister
%                → Segment → Normalise to MNI → Smooth.
%                Each step is skipped if its output already exists.
%
%  INPUTS      : data/ds000114/func/sub-01_ses-test_task-fingerfootlips_run-1_bold.nii[.gz]
%                data/ds000114/anat/sub-01_ses-test_T1w.nii[.gz]
%  OUTPUTS     : swrasub-01_..._bold.nii  (smoothed, normalised functional)
%                rp_asub-01_..._bold.txt  (6 motion parameters)
%                results/dataset1/preprocessing/*.mat, *.csv
%                figures/fig01–fig05_ds1_*.png/eps
%  DEPENDENCIES: SPM25, MATLAB R2022b+
%  =========================================================================

clear; clc; close all;
fprintf('=== Step 02: Dataset 1 Preprocessing (ds000114 — Motor Task) ===\n\n');

%% -----------------------------------------------------------------------
%  PARAMETERS
%  -----------------------------------------------------------------------
ROOT_DIR    = 'C:\fmri_project';
DATASET     = 'ds000114';
SUBJECT     = 'sub-01';
SESSION     = 'ses-test';
TASK        = 'fingerfootlips';
RUN         = 'run-1';

% Scan parameters (from OpenNeuro ds000114 dataset descriptor)
TR          = 2.5;       % Repetition time (seconds)
NUM_SLICES  = 40;        % Slices per volume
SMOOTH_FWHM = [6 6 6];   % Gaussian smoothing kernel (mm)

% Directory paths
func_dir = fullfile(ROOT_DIR,'data',DATASET,'func');
anat_dir = fullfile(ROOT_DIR,'data',DATASET,'anat');
OUT_DIR  = fullfile(ROOT_DIR,'results','dataset1','preprocessing');
FIG_DIR  = fullfile(ROOT_DIR,'figures');
LOG_FILE = fullfile(ROOT_DIR,'project_log.txt');

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Base filenames (no extension, no prefix)
fn_bold = [SUBJECT '_' SESSION '_task-' TASK '_' RUN '_bold'];
fn_anat = [SUBJECT '_' SESSION '_T1w'];

% Source files
func_gz  = fullfile(func_dir, [fn_bold '.nii.gz']);
func_nii = fullfile(func_dir, [fn_bold '.nii']);
anat_gz  = fullfile(anat_dir, [fn_anat '.nii.gz']);
anat_nii = fullfile(anat_dir, [fn_anat '.nii']);

% Derived filenames — SPM prefix stacking: a(ST) → ra(Realign) → wra(Norm) → swra(Smooth)
st_func   = fullfile(func_dir, ['a'      fn_bold '.nii']);   % slice-timed
ra_func   = fullfile(func_dir, ['ra'     fn_bold '.nii']);   % realigned
wra_func  = fullfile(func_dir, ['wra'    fn_bold '.nii']);   % normalised
swra_func = fullfile(func_dir, ['swra'   fn_bold '.nii']);   % smoothed  ← final output
mean_func = fullfile(func_dir, ['meana'  fn_bold '.nii']);   % mean image from realign
rp_func   = fullfile(func_dir, ['rp_a'   fn_bold '.txt']);   % motion parameters
y_field   = fullfile(anat_dir, ['y_'     fn_anat '.nii']);   % deformation field

% Initialise SPM
spm('defaults','FMRI');
spm_jobman('initcfg');

%% -----------------------------------------------------------------------
%  STEP 2.1 — DECOMPRESS SOURCE FILES
%  -----------------------------------------------------------------------
fprintf('--- [2.1] Decompressing source NIfTI files...\n');
if ~exist(func_nii,'file') && exist(func_gz,'file')
    fprintf('  Decompressing %s\n', func_gz);
    gunzip(func_gz, func_dir);
end
if ~exist(anat_nii,'file') && exist(anat_gz,'file')
    fprintf('  Decompressing %s\n', anat_gz);
    gunzip(anat_gz, anat_dir);
end
if ~exist(func_nii,'file'), error('Functional NIfTI not found: %s', func_nii); end
if ~exist(anat_nii,'file'), error('Anatomical NIfTI not found: %s', anat_nii); end

V_func = spm_vol(func_nii);
V_anat = spm_vol(anat_nii);
num_vols = numel(V_func);
fprintf('  Functional : %d volumes | voxels %dx%dx%d\n', ...
    num_vols, V_func(1).dim(1), V_func(1).dim(2), V_func(1).dim(3));
fprintf('  Anatomical : voxels %dx%dx%d\n', ...
    V_anat.dim(1), V_anat.dim(2), V_anat.dim(3));

% Auto-detect scan parameters from actual NIfTI data (overrides hardcoded values)
NUM_SLICES = V_func(1).dim(3);   % z-dimension = number of acquired slices
fprintf('  NUM_SLICES (from data): %d\n', NUM_SLICES);

% Try to read TR from NIfTI pixdim[4]
try
    hdr_tr = V_func(1).private.hdr.pixdim(5);
    if hdr_tr > 500, hdr_tr = hdr_tr / 1000; end   % ms → s
    if hdr_tr > 0.5 && hdr_tr < 30
        TR = hdr_tr;
        fprintf('  TR (from header): %.3f s\n', TR);
    else
        fprintf('  TR from header out of range — using %.2f s\n', TR);
    end
catch
    fprintf('  TR from header unavailable — using %.2f s\n', TR);
end
fprintf('  Final parameters: TR=%.3f s, NUM_SLICES=%d\n\n', TR, NUM_SLICES);

%% -----------------------------------------------------------------------
%  FIGURE 01 — RAW EPI SAMPLE (before any preprocessing)
%  -----------------------------------------------------------------------
fig01_name = 'fig01_ds1_raw_epi_sample';
if ~exist(fullfile(FIG_DIR, [fig01_name '.png']), 'file')
    fprintf('--- [2.2] Figure 01: Raw EPI sample slices...\n');
    fig = figure('Visible','off','Position',[100 100 1200 420]);
    cuts = [0.25 0.50 0.75];
    for k = 1:3
        sl = max(1, round(V_func(1).dim(3) * cuts(k)));
        subplot(1,3,k);
        img = spm_slice_vol(V_func(1), spm_matrix([0 0 sl]), V_func(1).dim(1:2), 0);
        imagesc(img'); colormap(gca,'gray'); axis image off; colorbar;
        title(sprintf('Axial slice %d', sl), 'FontSize', 10);
    end
    sgtitle(sprintf('Raw EPI — %s %s (Volume 1)', DATASET, SUBJECT), 'FontSize', 12);
    saveFig(fig, FIG_DIR, fig01_name);
    printLatex(fig01_name, ...
        ['Raw EPI sample (Volume 1) for Dataset 1 (ds000114), subject ' SUBJECT ...
         '. Three axial slices at 25\%, 50\%, and 75\% of the brain volume, ' ...
         'before any preprocessing.'], ...
        'ds1:raw_epi');
end

%% -----------------------------------------------------------------------
%  STEP 2.3 — SLICE TIMING CORRECTION
%  -----------------------------------------------------------------------
if ~exist(st_func, 'file')
    fprintf('--- [2.3] Slice Timing Correction (TR=%.2fs, %d slices)...\n', TR, NUM_SLICES);
    scans = cell(num_vols, 1);
    for v = 1:num_vols
        scans{v} = sprintf('%s,%d', func_nii, v);
    end
    clear matlabbatch;
    matlabbatch{1}.spm.temporal.st.scans{1}  = scans;
    matlabbatch{1}.spm.temporal.st.nslices   = NUM_SLICES;
    matlabbatch{1}.spm.temporal.st.tr        = TR;
    matlabbatch{1}.spm.temporal.st.ta        = TR - TR/NUM_SLICES;
    matlabbatch{1}.spm.temporal.st.so        = 1:NUM_SLICES;   % sequential ascending
    matlabbatch{1}.spm.temporal.st.refslice  = round(NUM_SLICES/2);
    matlabbatch{1}.spm.temporal.st.prefix    = 'a';
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', st_func);
else
    fprintf('--- [2.3] Slice Timing: output exists, skipping.\n');
end
V_st = spm_vol(st_func);
fprintf('  Volumes after slice timing: %d\n', numel(V_st));

%% -----------------------------------------------------------------------
%  STEP 2.4 — REALIGNMENT (Estimate & Reslice)
%  -----------------------------------------------------------------------
if ~exist(ra_func, 'file')
    fprintf('--- [2.4] Realignment (Estimate & Reslice)...\n');
    scans_ra = cell(numel(V_st), 1);
    for v = 1:numel(V_st)
        scans_ra{v} = sprintf('%s,%d', st_func, v);
    end
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.realign.estwrite.data{1}           = scans_ra;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality  = 0.9;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep      = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm     = 5;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm      = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp   = 2;
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.weight   = '';
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which    = [2 1];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp   = 4;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap     = [0 0 0];
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.mask     = 1;
    matlabbatch{1}.spm.spatial.realign.estwrite.roptions.prefix   = 'r';
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', ra_func);
else
    fprintf('--- [2.4] Realignment: output exists, skipping.\n');
end

% Load motion parameters (SPM writes rp_ file with same stem as input)
if ~exist(rp_func, 'file')
    rp_search = dir(fullfile(func_dir, 'rp_a*.txt'));
    if ~isempty(rp_search)
        rp_func = fullfile(func_dir, rp_search(1).name);
    end
end
if exist(rp_func, 'file')
    rp = load(rp_func);
    fprintf('  Motion params loaded: %d timepoints x 6 params\n', size(rp,1));
else
    rp = zeros(numel(V_st), 6);
    fprintf('  WARNING: Motion params file not found — using zeros.\n');
end
V_ra = spm_vol(ra_func);
fprintf('  Realigned volumes: %d\n', numel(V_ra));

%% -----------------------------------------------------------------------
%  FIGURE 02 — MOTION PARAMETERS
%  -----------------------------------------------------------------------
fig02_name = 'fig02_ds1_motion_parameters';
if ~exist(fullfile(FIG_DIR, [fig02_name '.png']), 'file')
    fprintf('--- [2.5] Figure 02: Motion parameters...\n');
    nv = size(rp, 1);
    fig = figure('Visible','off','Position',[100 100 1200 520]);
    subplot(2,1,1);
    plot(1:nv, rp(:,1), 'r-', 1:nv, rp(:,2), 'g-', 1:nv, rp(:,3), 'b-', 'LineWidth',1.5);
    xlabel('Volume'); ylabel('Translation (mm)');
    legend({'x','y','z'}, 'Location','northeast');
    title('Translations'); grid on; xlim([1 nv]);

    subplot(2,1,2);
    plot(1:nv, rp(:,4)*180/pi, 'r-', 1:nv, rp(:,5)*180/pi, 'g-', ...
         1:nv, rp(:,6)*180/pi, 'b-', 'LineWidth',1.5);
    xlabel('Volume'); ylabel('Rotation (degrees)');
    legend({'pitch','roll','yaw'}, 'Location','northeast');
    title('Rotations'); grid on; xlim([1 nv]);

    sgtitle(sprintf('Head Motion Parameters — %s %s', DATASET, SUBJECT), 'FontSize',12);
    saveFig(fig, FIG_DIR, fig02_name);
    printLatex(fig02_name, ...
        ['Six rigid-body head-motion parameters for Dataset 1 (ds000114), subject ' SUBJECT ...
         '. Top panel: three translations in mm ($x$, $y$, $z$). ' ...
         'Bottom panel: three rotations in degrees (pitch, roll, yaw).'], ...
        'ds1:motion_params');
end

%% -----------------------------------------------------------------------
%  STEP 2.6 — COREGISTRATION (mean EPI → T1 anatomical)
%  -----------------------------------------------------------------------
% Locate mean functional image produced by realignment
if ~exist(mean_func, 'file')
    mf_search = dir(fullfile(func_dir, 'mean*.nii'));
    if ~isempty(mf_search)
        mean_func = fullfile(func_dir, mf_search(1).name);
        fprintf('  Mean image found: %s\n', mf_search(1).name);
    else
        % Fall back to first realigned volume as reference
        mean_func = ra_func;
        fprintf('  WARNING: mean image not found, using first volume as reference.\n');
    end
end

fprintf('--- [2.6] Coregistration (mean EPI → T1)...\n');
clear matlabbatch;
matlabbatch{1}.spm.spatial.coreg.estimate.ref              = {[anat_nii ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source           = {[mean_func ',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other            = {ra_func};
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep     = [4 2];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol     = ...
    [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm    = [7 7];
spm_jobman('run', matlabbatch);
fprintf('  Coregistration complete (headers updated).\n');

%% -----------------------------------------------------------------------
%  FIGURE 03 — COREGISTRATION CHECK (T1 + mean EPI overlay)
%  -----------------------------------------------------------------------
fig03_name = 'fig03_ds1_coregistration';
if ~exist(fullfile(FIG_DIR, [fig03_name '.png']), 'file')
    fprintf('--- [2.7] Figure 03: Coregistration overlay check...\n');
    V_anat_v = spm_vol(anat_nii);
    V_mean_v = spm_vol(mean_func);
    fig = figure('Visible','off','Position',[100 100 1200 420]);
    fracs = [0.30 0.50 0.70];
    for k = 1:3
        sl = max(1, round(V_anat_v.dim(3) * fracs(k)));
        subplot(1,3,k);
        as = spm_slice_vol(V_anat_v, spm_matrix([0 0 sl]), V_anat_v.dim(1:2), 1);
        ms = spm_slice_vol(V_mean_v, spm_matrix([0 0 sl]), V_anat_v.dim(1:2), 1);
        as = norm01(as);
        ms = norm01(ms);
        overlay = cat(3, as*0.7 + ms*0.3, as*0.7, as*0.7);
        imagesc(overlay); axis image off;
        title(sprintf('Axial z=%d', sl), 'FontSize',10);
    end
    sgtitle('Coregistration Check: T1 (gray) + mean EPI (red)', 'FontSize',12);
    saveFig(fig, FIG_DIR, fig03_name);
    printLatex(fig03_name, ...
        ['Coregistration check for Dataset 1. Three axial slices showing T1 anatomical ' ...
         '(greyscale) with mean EPI overlaid in red. Correct alignment is indicated by ' ...
         'brain-boundary correspondence between the two modalities.'], ...
        'ds1:coregistration');
end

%% -----------------------------------------------------------------------
%  STEP 2.8 — SEGMENTATION (New Segment — generates deformation field)
%  -----------------------------------------------------------------------
if ~exist(y_field, 'file')
    fprintf('--- [2.8] Segmentation (New Segment)...\n');
    tpm_path = fullfile(spm('Dir'), 'tpm', 'TPM.nii');
    ngaus_vals = [1 1 2 3 4 2];
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[anat_nii ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write    = [1 1];
    for t = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm    = ...
            {[tpm_path ',' num2str(t)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus  = ngaus_vals(t);
        matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [double(t<6) 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
    end
    matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg     = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg  = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm    = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp    = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write   = [1 1];
    spm_jobman('run', matlabbatch);
    fprintf('  Segmentation complete → %s\n', y_field);
else
    fprintf('--- [2.8] Segmentation: deformation field exists, skipping.\n');
end

if ~exist(y_field, 'file')
    % Search for any y_*.nii in anat dir
    df_s = dir(fullfile(anat_dir, 'y_*.nii'));
    if ~isempty(df_s)
        y_field = fullfile(anat_dir, df_s(1).name);
        fprintf('  Found deformation field: %s\n', df_s(1).name);
    else
        error('Deformation field y_*.nii not found in %s', anat_dir);
    end
end

%% -----------------------------------------------------------------------
%  STEP 2.9 — NORMALISE TO MNI SPACE (2mm isotropic)
%  -----------------------------------------------------------------------
if ~exist(wra_func, 'file')
    fprintf('--- [2.9] Normalisation to MNI space...\n');
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {y_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {ra_func};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    fprintf('  Normalisation complete → %s\n', wra_func);
else
    fprintf('--- [2.9] Normalisation: output exists, skipping.\n');
end
V_wra = spm_vol(wra_func);
fprintf('  Normalised: %d volumes | voxels %dx%dx%d\n', ...
    numel(V_wra), V_wra(1).dim(1), V_wra(1).dim(2), V_wra(1).dim(3));

%% -----------------------------------------------------------------------
%  FIGURE 04 — tSNR BEFORE vs AFTER NORMALISATION
%  -----------------------------------------------------------------------
fig04_name = 'fig04_ds1_tsnr_comparison';
if ~exist(fullfile(FIG_DIR, [fig04_name '.png']), 'file')
    fprintf('--- [2.10] Figure 04: tSNR before/after normalisation...\n');
    fprintf('  Computing tSNR maps (may take a moment)...\n');
    tsnr_raw  = compute_tsnr(V_ra);
    tsnr_norm = compute_tsnr(V_wra);

    ms_ra  = min(round(size(tsnr_raw,3)/2),  size(tsnr_raw,3));
    ms_wra = min(round(size(tsnr_norm,3)/2), size(tsnr_norm,3));
    clim = [0 80];

    fig = figure('Visible','off','Position',[100 100 1200 500]);
    subplot(1,2,1);
    imagesc(tsnr_raw(:,:,ms_ra)'); caxis(clim);
    colormap parula; axis image off; colorbar;
    title(sprintf('tSNR Before Normalise (slice %d)', ms_ra), 'FontSize',10);

    subplot(1,2,2);
    imagesc(tsnr_norm(:,:,ms_wra)'); caxis(clim);
    colormap parula; axis image off; colorbar;
    title(sprintf('tSNR After Normalise (slice %d)', ms_wra), 'FontSize',10);

    sgtitle(sprintf('Temporal SNR Maps — %s %s', DATASET, SUBJECT), 'FontSize',12);
    saveFig(fig, FIG_DIR, fig04_name);
    printLatex(fig04_name, ...
        ['Temporal signal-to-noise ratio (tSNR = mean/std across time) for Dataset 1 ' ...
         'before (left) and after (right) MNI normalisation. ' ...
         'Higher values (yellow) indicate more stable signal.'], ...
        'ds1:tsnr_norm');
end

%% -----------------------------------------------------------------------
%  STEP 2.11 — SMOOTHING (Gaussian kernel FWHM 6mm)
%  -----------------------------------------------------------------------
if ~exist(swra_func, 'file')
    fprintf('--- [2.11] Smoothing (FWHM = [%s] mm)...\n', num2str(SMOOTH_FWHM));
    nv = numel(V_wra);
    files_smooth = cell(nv, 1);
    for v = 1:nv
        files_smooth{v} = sprintf('%s,%d', wra_func, v);
    end
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.smooth.data   = files_smooth;
    matlabbatch{1}.spm.spatial.smooth.fwhm   = SMOOTH_FWHM;
    matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{1}.spm.spatial.smooth.im     = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
    fprintf('  Smoothing complete → %s\n', swra_func);
else
    fprintf('--- [2.11] Smoothing: output exists, skipping.\n');
end
V_smooth = spm_vol(swra_func);
fprintf('  Smoothed: %d volumes | voxels %dx%dx%d\n', ...
    numel(V_smooth), V_smooth(1).dim(1), V_smooth(1).dim(2), V_smooth(1).dim(3));

%% -----------------------------------------------------------------------
%  FIGURE 05 — SMOOTHING EFFECT ON tSNR
%  -----------------------------------------------------------------------
fig05_name = 'fig05_ds1_smoothing_effect';
if ~exist(fullfile(FIG_DIR, [fig05_name '.png']), 'file')
    fprintf('--- [2.12] Figure 05: Smoothing effect on tSNR...\n');
    fprintf('  Computing tSNR for smoothed data...\n');
    tsnr_smooth = compute_tsnr(V_smooth);

    % Re-use tsnr_norm if already in workspace, else recompute
    if ~exist('tsnr_norm','var')
        tsnr_norm = compute_tsnr(V_wra);
    end

    ms2 = min(round(size(tsnr_smooth,3)/2), size(tsnr_smooth,3));
    sl_norm   = tsnr_norm(:,:,min(ms2,size(tsnr_norm,3)))';
    sl_smooth = tsnr_smooth(:,:,ms2)';
    sl_diff   = sl_smooth - sl_norm;
    clim1 = [0 max(quantile(sl_norm(:),0.99), quantile(sl_smooth(:),0.99))];

    fig = figure('Visible','off','Position',[100 100 1400 420]);
    subplot(1,3,1);
    imagesc(sl_norm);  caxis(clim1); colormap parula; axis image off; colorbar;
    title('tSNR Before Smoothing', 'FontSize',10);
    subplot(1,3,2);
    imagesc(sl_smooth); caxis(clim1); colormap parula; axis image off; colorbar;
    title(sprintf('tSNR After Smoothing (FWHM=%dmm)', SMOOTH_FWHM(1)), 'FontSize',10);
    subplot(1,3,3);
    imagesc(sl_diff);  colormap jet; axis image off; colorbar;
    title('Difference (Smooth - Unsmooth)', 'FontSize',10);

    sgtitle(sprintf('Smoothing Effect — %s %s', DATASET, SUBJECT), 'FontSize',12);
    saveFig(fig, FIG_DIR, fig05_name);
    printLatex(fig05_name, ...
        ['Effect of Gaussian spatial smoothing (FWHM = 6\,mm) on temporal SNR for Dataset 1. ' ...
         'Left: tSNR before smoothing. Centre: tSNR after smoothing. ' ...
         'Right: difference map (positive values indicate tSNR gain from smoothing).'], ...
        'ds1:smoothing');
end

%% -----------------------------------------------------------------------
%  STEP 2.13 — VALIDATION
%  -----------------------------------------------------------------------
fprintf('\n--- [2.13] Validation...\n');
all_pass = true;

% Check key output files
checks = {
    swra_func, sprintf('Smoothed output: %s', ['swra' fn_bold '.nii'])
    rp_func,   'Motion parameters file'
    y_field,   'Deformation field'
};
for i = 1:size(checks,1)
    if exist(checks{i,1},'file')
        info = dir(checks{i,1});
        fprintf('  PASS: %s (%.1f MB)\n', checks{i,2}, info.bytes/1e6);
    else
        fprintf('  FAIL: %s — NOT FOUND\n', checks{i,2});
        all_pass = false;
    end
end

% Volume count
V_final = spm_vol(swra_func);
if numel(V_final) == num_vols
    fprintf('  PASS: Volume count = %d\n', numel(V_final));
else
    fprintf('  WARN: Expected %d volumes, got %d\n', num_vols, numel(V_final));
end

% Voxel size (should be ~2mm isotropic after normalisation)
vs = sqrt(sum(V_final(1).mat(1:3,1:3).^2));
fprintf('  INFO: Voxel size = [%.2f %.2f %.2f] mm\n', vs(1), vs(2), vs(3));

% Motion quality check
max_t_mm  = max(abs(rp(:,1:3)), [], 'all');
max_r_deg = max(abs(rp(:,4:6)) * 180/pi, [], 'all');
fd = sum(abs(diff(rp(:,1:3))),2) + sum(abs(diff(rp(:,4:6)))*50,2);  % approximate FD
fprintf('  INFO: Max translation = %.2f mm | Max rotation = %.2f deg\n', max_t_mm, max_r_deg);
fprintf('  INFO: Mean framewise displacement = %.2f mm\n', mean(fd));
if max_t_mm > 3.0
    fprintf('  WARN: Large head motion (>3mm)\n');
end
if max_r_deg > 3.0
    fprintf('  WARN: Large rotation (>3 deg)\n');
end

status_str = iif(all_pass, 'PASS', 'WARNINGS');
fprintf('\n=== PREPROCESSING COMPLETE — Dataset 1 | Status: %s ===\n\n', status_str);

%% -----------------------------------------------------------------------
%  STEP 2.14 — SAVE RESULTS TO RESULTS FOLDER
%  -----------------------------------------------------------------------
fprintf('--- [2.14] Saving results...\n');

% Copy key outputs
dest_func = fullfile(OUT_DIR, 'swra_bold_ds1.nii');
if ~exist(dest_func,'file') || dir(dest_func).bytes < dir(swra_func).bytes
    copyfile(swra_func, dest_func);
    fprintf('  Copied smoothed func to results/\n');
end
if exist(rp_func,'file')
    copyfile(rp_func, fullfile(OUT_DIR,'rp_bold_ds1.txt'));
end

% Summary .mat
preproc = struct();
preproc.dataset       = DATASET;
preproc.subject       = SUBJECT;
preproc.TR            = TR;
preproc.num_slices    = NUM_SLICES;
preproc.num_volumes   = num_vols;
preproc.smooth_fwhm   = SMOOTH_FWHM;
preproc.smoothed_file = swra_func;
preproc.max_trans_mm  = max_t_mm;
preproc.max_rot_deg   = max_r_deg;
preproc.mean_fd_mm    = mean(fd);
preproc.status        = status_str;
preproc.timestamp     = datestr(now);
save(fullfile(OUT_DIR, 'result_ds1_preprocessing_summary.mat'), '-struct', 'preproc');
fprintf('  Saved: result_ds1_preprocessing_summary.mat\n');

% Motion CSV
csv_file = fullfile(OUT_DIR, 'table_ds1_motion_parameters.csv');
fid = fopen(csv_file, 'w');
fprintf(fid, 'Volume,TransX_mm,TransY_mm,TransZ_mm,RotPitch_deg,RotRoll_deg,RotYaw_deg,FD_mm\n');
fd_all = [0; fd];  % prepend 0 for first volume
for v = 1:size(rp,1)
    fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
        v, rp(v,1), rp(v,2), rp(v,3), ...
        rp(v,4)*180/pi, rp(v,5)*180/pi, rp(v,6)*180/pi, ...
        fd_all(v));
end
fclose(fid);
fprintf('  Saved: table_ds1_motion_parameters.csv\n');

%% -----------------------------------------------------------------------
%  PROJECT LOG
%  -----------------------------------------------------------------------
fid = fopen(LOG_FILE, 'a');
fprintf(fid, '========================================================\n');
fprintf(fid, 'STEP 02 : Dataset 1 Preprocessing\n');
fprintf(fid, 'Date    : %s\n', datestr(now));
fprintf(fid, 'Dataset : %s | Subject: %s\n', DATASET, SUBJECT);
fprintf(fid, 'TR      : %.2f s | Slices: %d | Volumes: %d\n', TR, NUM_SLICES, num_vols);
fprintf(fid, 'Output  : %s\n', swra_func);
fprintf(fid, 'Motion  : max=%.2f mm, max=%.2f deg, meanFD=%.2f\n', ...
    max_t_mm, max_r_deg, mean(fd));
fprintf(fid, 'STATUS  : %s\n', status_str);
fprintf(fid, '========================================================\n\n');
fclose(fid);
fprintf('  Appended to project_log.txt\n');

%% -----------------------------------------------------------------------
%  GIT COMMIT COMMANDS
%  -----------------------------------------------------------------------
fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', fullfile(ROOT_DIR,'github'));
fprintf('git add .\n');
fprintf('git commit -m "Step 02 complete: Dataset 1 preprocessing (%s)"\n', status_str);
fprintf('git push origin main\n');
fprintf('\n=== END OF STEP 02 ===\n');

%% -----------------------------------------------------------------------
%  LOCAL HELPER FUNCTIONS
%  -----------------------------------------------------------------------

function tsnr = compute_tsnr(V)
    nv = numel(V);
    d0 = spm_read_vols(V(1));
    [nx,ny,nz] = size(d0);
    all_vols = zeros(nx, ny, nz, nv, 'single');
    for k = 1:nv
        all_vols(:,:,:,k) = single(spm_read_vols(V(k)));
    end
    mu = mean(all_vols, 4);
    sg = std(all_vols, 0, 4);
    sg(sg < single(eps)) = single(eps);
    tsnr = double(mu ./ sg);
end

function saveFig(fig, fig_dir, name)
    print(fig, fullfile(fig_dir, [name '.png']), '-dpng', '-r150');
    print(fig, fullfile(fig_dir, [name '.eps']), '-depsc');
    close(fig);
    fprintf('  Saved: %s.png / .eps\n', name);
end

function out = norm01(img)
    mn = min(img(:));  mx = max(img(:));
    out = (img - mn) / (mx - mn + eps);
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
