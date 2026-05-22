%% =========================================================================
%  FILE        : step03_dataset2_preprocessing.m
%  TITLE       : fMRI Brain Data Analysis using SPM and MATLAB —
%                Preprocessing, GLM, Statistical Mapping and
%                Functional Connectivity
%  AUTHOR      : Fahim Ur Rehman Shah
%  REG NO      : EE2629
%  COURSE      : CSE532 — Signal and Image Processing (MS Level)
%  SUPERVISOR  : Dr. Adnan Shah, FCSE, GIKI
%  INSTITUTION : GIK Institute of Engineering Sciences & Technology
%  DATE        : May 2026
%  STEP        : 03 of 07
%
%  DESCRIPTION : Full SPM preprocessing pipeline for Dataset 2 (ds000105,
%                Visual object recognition task). Steps executed in order:
%                Slice Timing → Realign → Coregister → Segment
%                → Normalise → Smooth.
%                Same pipeline as step02 but adapted for ds000105 filenames.
%
%  INPUTS      : data/ds000105/func/sub-1_task-objectviewing_run-1_bold.nii[.gz]
%                data/ds000105/anat/sub-1_T1w.nii[.gz]
%  OUTPUTS     : swrasub-1_..._bold.nii (smoothed, normalised functional)
%                rp_asub-1_..._bold.txt (6 motion parameters)
%                results/dataset2/preprocessing/*.mat, *.csv
%                figures/fig06–fig10_ds2_*.png/eps
%  DEPENDENCIES: SPM25, MATLAB R2022b+
%  =========================================================================

clear; clc; close all;
fprintf('=== Step 03: Dataset 2 Preprocessing (ds000105 — Visual Task) ===\n\n');

%% -----------------------------------------------------------------------
%  PARAMETERS
%  -----------------------------------------------------------------------
ROOT_DIR    = 'C:\fmri_project';
DATASET     = 'ds000105';
SUBJECT     = 'sub-1';       % NOTE: no leading zero in ds000105
TASK        = 'objectviewing';
RUN         = 'run-1';

% Scan parameters (from OpenNeuro ds000105 dataset descriptor)
TR          = 2.5;       % Repetition time (seconds)
NUM_SLICES  = 40;        % Slices per volume
SMOOTH_FWHM = [6 6 6];   % Gaussian smoothing kernel (mm)

% Directory paths
func_dir = fullfile(ROOT_DIR,'data',DATASET,'func');
anat_dir = fullfile(ROOT_DIR,'data',DATASET,'anat');
OUT_DIR  = fullfile(ROOT_DIR,'results','dataset2','preprocessing');
FIG_DIR  = fullfile(ROOT_DIR,'figures');
LOG_FILE = fullfile(ROOT_DIR,'project_log.txt');

if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

% Base filenames (no extension, no prefix) — ds000105 has no session label
fn_bold = [SUBJECT '_task-' TASK '_' RUN '_bold'];
fn_anat = [SUBJECT '_T1w'];

% Source files
func_gz  = fullfile(func_dir, [fn_bold '.nii.gz']);
func_nii = fullfile(func_dir, [fn_bold '.nii']);
anat_gz  = fullfile(anat_dir, [fn_anat '.nii.gz']);
anat_nii = fullfile(anat_dir, [fn_anat '.nii']);

% Check if anat is also in func dir (some datasets store it there)
if ~exist(anat_gz,'file') && ~exist(anat_nii,'file')
    anat_gz2  = fullfile(func_dir, [fn_anat '.nii.gz']);
    anat_nii2 = fullfile(func_dir, [fn_anat '.nii']);
    if exist(anat_gz2,'file')
        anat_gz  = anat_gz2;
        anat_nii = anat_nii2;
        anat_dir = func_dir;  % update anat_dir
        fprintf('  NOTE: T1w found in func/ directory.\n');
    end
end

% Derived filenames — SPM prefix convention
st_func   = fullfile(func_dir, ['a'      fn_bold '.nii']);
ra_func   = fullfile(func_dir, ['ra'     fn_bold '.nii']);
wra_func  = fullfile(func_dir, ['wra'    fn_bold '.nii']);
swra_func = fullfile(func_dir, ['swra'   fn_bold '.nii']);
mean_func = fullfile(func_dir, ['meana'  fn_bold '.nii']);
rp_func   = fullfile(func_dir, ['rp_a'   fn_bold '.txt']);
y_field   = fullfile(anat_dir, ['y_'     fn_anat '.nii']);

% Initialise SPM
spm('defaults','FMRI');
spm_jobman('initcfg');

%% -----------------------------------------------------------------------
%  STEP 3.1 — DECOMPRESS SOURCE FILES
%  -----------------------------------------------------------------------
fprintf('--- [3.1] Decompressing source NIfTI files...\n');
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

% Try to read TR from NIfTI pixdim[4] (SPM stores it 1-indexed as pixdim(5))
try
    hdr_tr = V_func(1).private.hdr.pixdim(5);
    if hdr_tr > 500            % header stores ms — convert
        hdr_tr = hdr_tr / 1000;
        fprintf('  TR from header (converted from ms): %.3f s\n', hdr_tr);
    end
    if hdr_tr > 0.5 && hdr_tr < 30   % sanity: must be 0.5–30 s
        TR = hdr_tr;
        fprintf('  TR (from header): %.3f s\n', TR);
    else
        fprintf('  TR from header out of range (%.3f) — using %.2f s\n', hdr_tr, TR);
    end
catch
    fprintf('  TR from header unavailable — using %.2f s\n', TR);
end
fprintf('  Final parameters: TR=%.3f s, NUM_SLICES=%d\n\n', TR, NUM_SLICES);

%% -----------------------------------------------------------------------
%  FIGURE 06 — RAW EPI SAMPLE (Dataset 2)
%  -----------------------------------------------------------------------
fig06_name = 'fig06_ds2_raw_epi_sample';
if ~exist(fullfile(FIG_DIR, [fig06_name '.png']), 'file')
    fprintf('--- [3.2] Figure 06: Raw EPI sample slices...\n');
    fig = figure('Visible','off','Position',[100 100 1200 420]);
    cuts = [0.25 0.50 0.75];
    for k = 1:3
        sl = max(1, round(V_func(1).dim(3) * cuts(k)));
        subplot(1,3,k);
        img = spm_slice_vol(V_func(1), spm_matrix([0 0 sl]), V_func(1).dim(1:2), 0);
        imagesc(img'); colormap(gca,'gray'); axis image off; colorbar;
        title(sprintf('Axial slice %d', sl), 'FontSize',10);
    end
    sgtitle(sprintf('Raw EPI — %s %s (Volume 1)', DATASET, SUBJECT), 'FontSize',12);
    saveFig(fig, FIG_DIR, fig06_name);
    printLatex(fig06_name, ...
        ['Raw EPI sample (Volume 1) for Dataset 2 (ds000105), subject ' SUBJECT ...
         '. Three axial slices at 25\%, 50\%, and 75\% of the brain volume, ' ...
         'before any preprocessing.'], ...
        'ds2:raw_epi');
end

%% -----------------------------------------------------------------------
%  STEP 3.3 — SLICE TIMING CORRECTION
%  -----------------------------------------------------------------------
if ~exist(st_func, 'file')
    fprintf('--- [3.3] Slice Timing Correction (TR=%.2fs, %d slices)...\n', TR, NUM_SLICES);
    scans = cell(num_vols, 1);
    for v = 1:num_vols
        scans{v} = sprintf('%s,%d', func_nii, v);
    end
    clear matlabbatch;
    matlabbatch{1}.spm.temporal.st.scans{1}  = scans;
    matlabbatch{1}.spm.temporal.st.nslices   = NUM_SLICES;
    matlabbatch{1}.spm.temporal.st.tr        = TR;
    matlabbatch{1}.spm.temporal.st.ta        = TR - TR/NUM_SLICES;
    matlabbatch{1}.spm.temporal.st.so        = 1:NUM_SLICES;
    matlabbatch{1}.spm.temporal.st.refslice  = round(NUM_SLICES/2);
    matlabbatch{1}.spm.temporal.st.prefix    = 'a';
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', st_func);
else
    fprintf('--- [3.3] Slice Timing: output exists, skipping.\n');
end
V_st = spm_vol(st_func);
fprintf('  Volumes after slice timing: %d\n', numel(V_st));

%% -----------------------------------------------------------------------
%  STEP 3.4 — REALIGNMENT (Estimate & Reslice)
%  -----------------------------------------------------------------------
if ~exist(ra_func, 'file')
    fprintf('--- [3.4] Realignment (Estimate & Reslice)...\n');
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
    fprintf('--- [3.4] Realignment: output exists, skipping.\n');
end

% Load motion parameters
if ~exist(rp_func,'file')
    rp_search = dir(fullfile(func_dir,'rp_a*.txt'));
    if ~isempty(rp_search)
        rp_func = fullfile(func_dir, rp_search(1).name);
    end
end
if exist(rp_func,'file')
    rp = load(rp_func);
    fprintf('  Motion params: %d x 6\n', size(rp,1));
else
    rp = zeros(numel(V_st), 6);
    fprintf('  WARNING: Motion params not found — using zeros.\n');
end
V_ra = spm_vol(ra_func);

%% -----------------------------------------------------------------------
%  FIGURE 07 — MOTION PARAMETERS (Dataset 2)
%  -----------------------------------------------------------------------
fig07_name = 'fig07_ds2_motion_parameters';
if ~exist(fullfile(FIG_DIR, [fig07_name '.png']), 'file')
    fprintf('--- [3.5] Figure 07: Motion parameters...\n');
    nv = size(rp,1);
    fig = figure('Visible','off','Position',[100 100 1200 520]);
    subplot(2,1,1);
    plot(1:nv, rp(:,1),'r-', 1:nv, rp(:,2),'g-', 1:nv, rp(:,3),'b-','LineWidth',1.5);
    xlabel('Volume'); ylabel('Translation (mm)');
    legend({'x','y','z'},'Location','northeast');
    title('Translations'); grid on; xlim([1 nv]);
    subplot(2,1,2);
    plot(1:nv, rp(:,4)*180/pi,'r-', 1:nv, rp(:,5)*180/pi,'g-', ...
         1:nv, rp(:,6)*180/pi,'b-','LineWidth',1.5);
    xlabel('Volume'); ylabel('Rotation (degrees)');
    legend({'pitch','roll','yaw'},'Location','northeast');
    title('Rotations'); grid on; xlim([1 nv]);
    sgtitle(sprintf('Head Motion Parameters — %s %s', DATASET, SUBJECT),'FontSize',12);
    saveFig(fig, FIG_DIR, fig07_name);
    printLatex(fig07_name, ...
        ['Six rigid-body head-motion parameters for Dataset 2 (ds000105), subject ' SUBJECT ...
         '. Top: three translations in mm. Bottom: three rotations in degrees.'], ...
        'ds2:motion_params');
end

%% -----------------------------------------------------------------------
%  STEP 3.6 — COREGISTRATION
%  -----------------------------------------------------------------------
if ~exist(mean_func,'file')
    mf_s = dir(fullfile(func_dir,'mean*.nii'));
    if ~isempty(mf_s), mean_func = fullfile(func_dir,mf_s(1).name); end
end
if ~exist(mean_func,'file'), mean_func = ra_func; end

fprintf('--- [3.6] Coregistration (mean EPI → T1)...\n');
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
fprintf('  Coregistration complete.\n');

%% -----------------------------------------------------------------------
%  FIGURE 08 — COREGISTRATION CHECK (Dataset 2)
%  -----------------------------------------------------------------------
fig08_name = 'fig08_ds2_coregistration';
if ~exist(fullfile(FIG_DIR, [fig08_name '.png']), 'file')
    fprintf('--- [3.7] Figure 08: Coreg overlay check...\n');
    V_anat_v = spm_vol(anat_nii);
    V_mean_v = spm_vol(mean_func);
    fig = figure('Visible','off','Position',[100 100 1200 420]);
    fracs = [0.30 0.50 0.70];
    for k = 1:3
        sl = max(1, round(V_anat_v.dim(3)*fracs(k)));
        subplot(1,3,k);
        as = spm_slice_vol(V_anat_v, spm_matrix([0 0 sl]), V_anat_v.dim(1:2), 1);
        ms = spm_slice_vol(V_mean_v, spm_matrix([0 0 sl]), V_anat_v.dim(1:2), 1);
        as = norm01(as); ms = norm01(ms);
        overlay = cat(3, as*0.7+ms*0.3, as*0.7, as*0.7);
        imagesc(overlay); axis image off;
        title(sprintf('Axial z=%d',sl),'FontSize',10);
    end
    sgtitle('Coregistration Check: T1 (gray) + mean EPI (red) — DS2','FontSize',12);
    saveFig(fig, FIG_DIR, fig08_name);
    printLatex(fig08_name, ...
        ['Coregistration check for Dataset 2 (ds000105). T1 anatomical (greyscale) ' ...
         'with mean EPI overlaid in red after mutual-information registration.'], ...
        'ds2:coregistration');
end

%% -----------------------------------------------------------------------
%  STEP 3.8 — SEGMENTATION
%  -----------------------------------------------------------------------
if ~exist(y_field,'file')
    fprintf('--- [3.8] Segmentation (New Segment)...\n');
    tpm_path = fullfile(spm('Dir'),'tpm','TPM.nii');
    ngaus_vals = [1 1 2 3 4 2];
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.preproc.channel.vols     = {[anat_nii ',1']};
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.preproc.channel.write    = [1 1];
    for t = 1:6
        matlabbatch{1}.spm.spatial.preproc.tissue(t).tpm   = {[tpm_path ',' num2str(t)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(t).ngaus = ngaus_vals(t);
        matlabbatch{1}.spm.spatial.preproc.tissue(t).native = [double(t<6) 0];
        matlabbatch{1}.spm.spatial.preproc.tissue(t).warped = [0 0];
    end
    matlabbatch{1}.spm.spatial.preproc.warp.mrf    = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.reg    = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm   = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp   = 3;
    matlabbatch{1}.spm.spatial.preproc.warp.write  = [1 1];
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', y_field);
else
    fprintf('--- [3.8] Segmentation: deformation field exists, skipping.\n');
end
if ~exist(y_field,'file')
    df_s = dir(fullfile(anat_dir,'y_*.nii'));
    if ~isempty(df_s), y_field = fullfile(anat_dir,df_s(1).name); end
end

%% -----------------------------------------------------------------------
%  STEP 3.9 — NORMALISE TO MNI
%  -----------------------------------------------------------------------
if ~exist(wra_func,'file')
    fprintf('--- [3.9] Normalisation to MNI...\n');
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {y_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = {ra_func};
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = [-78 -112 -70; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', wra_func);
else
    fprintf('--- [3.9] Normalisation: output exists, skipping.\n');
end
V_wra = spm_vol(wra_func);
fprintf('  Normalised: %d vols | %dx%dx%d voxels\n', ...
    numel(V_wra), V_wra(1).dim(1), V_wra(1).dim(2), V_wra(1).dim(3));

%% -----------------------------------------------------------------------
%  FIGURE 09 — tSNR COMPARISON (Dataset 2)
%  -----------------------------------------------------------------------
fig09_name = 'fig09_ds2_tsnr_comparison';
if ~exist(fullfile(FIG_DIR, [fig09_name '.png']), 'file')
    fprintf('--- [3.10] Figure 09: tSNR before/after normalisation...\n');
    tsnr_raw  = compute_tsnr(V_ra);
    tsnr_norm = compute_tsnr(V_wra);
    ms1 = round(size(tsnr_raw,3)/2);
    ms2 = round(size(tsnr_norm,3)/2);
    clim = [0 80];
    fig = figure('Visible','off','Position',[100 100 1200 500]);
    subplot(1,2,1);
    imagesc(tsnr_raw(:,:,min(ms1,size(tsnr_raw,3)))'); caxis(clim);
    colormap parula; axis image off; colorbar;
    title(sprintf('tSNR Before Normalise (slice %d)',ms1),'FontSize',10);
    subplot(1,2,2);
    imagesc(tsnr_norm(:,:,ms2)'); caxis(clim);
    colormap parula; axis image off; colorbar;
    title(sprintf('tSNR After Normalise (slice %d)',ms2),'FontSize',10);
    sgtitle(sprintf('Temporal SNR — %s %s', DATASET, SUBJECT),'FontSize',12);
    saveFig(fig, FIG_DIR, fig09_name);
    printLatex(fig09_name, ...
        ['Temporal SNR for Dataset 2 (ds000105) before (left) and after (right) ' ...
         'MNI normalisation.'], 'ds2:tsnr_norm');
end

%% -----------------------------------------------------------------------
%  STEP 3.11 — SMOOTHING
%  -----------------------------------------------------------------------
if ~exist(swra_func,'file')
    fprintf('--- [3.11] Smoothing (FWHM=[%s] mm)...\n', num2str(SMOOTH_FWHM));
    nv = numel(V_wra);
    files_smooth = cell(nv,1);
    for v = 1:nv, files_smooth{v} = sprintf('%s,%d', wra_func, v); end
    clear matlabbatch;
    matlabbatch{1}.spm.spatial.smooth.data   = files_smooth;
    matlabbatch{1}.spm.spatial.smooth.fwhm   = SMOOTH_FWHM;
    matlabbatch{1}.spm.spatial.smooth.dtype  = 0;
    matlabbatch{1}.spm.spatial.smooth.im     = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
    fprintf('  Done → %s\n', swra_func);
else
    fprintf('--- [3.11] Smoothing: output exists, skipping.\n');
end
V_smooth = spm_vol(swra_func);
fprintf('  Smoothed: %d vols | %dx%dx%d voxels\n', ...
    numel(V_smooth), V_smooth(1).dim(1), V_smooth(1).dim(2), V_smooth(1).dim(3));

%% -----------------------------------------------------------------------
%  FIGURE 10 — SMOOTHING EFFECT (Dataset 2)
%  -----------------------------------------------------------------------
fig10_name = 'fig10_ds2_smoothing_effect';
if ~exist(fullfile(FIG_DIR, [fig10_name '.png']), 'file')
    fprintf('--- [3.12] Figure 10: Smoothing effect on tSNR...\n');
    tsnr_smooth = compute_tsnr(V_smooth);
    if ~exist('tsnr_norm','var'), tsnr_norm = compute_tsnr(V_wra); end
    ms3 = round(size(tsnr_smooth,3)/2);
    sl_n = tsnr_norm(:,:,min(ms3,size(tsnr_norm,3)))';
    sl_s = tsnr_smooth(:,:,ms3)';
    sl_d = sl_s - sl_n;
    clim1 = [0 max(quantile(sl_n(:),0.99), quantile(sl_s(:),0.99))];
    fig = figure('Visible','off','Position',[100 100 1400 420]);
    subplot(1,3,1); imagesc(sl_n);  caxis(clim1); colormap parula; axis image off; colorbar;
    title('tSNR Before Smoothing','FontSize',10);
    subplot(1,3,2); imagesc(sl_s); caxis(clim1); colormap parula; axis image off; colorbar;
    title(sprintf('tSNR After Smoothing (FWHM=%dmm)',SMOOTH_FWHM(1)),'FontSize',10);
    subplot(1,3,3); imagesc(sl_d); colormap jet; axis image off; colorbar;
    title('Difference','FontSize',10);
    sgtitle(sprintf('Smoothing Effect — %s %s', DATASET, SUBJECT),'FontSize',12);
    saveFig(fig, FIG_DIR, fig10_name);
    printLatex(fig10_name, ...
        ['Effect of Gaussian smoothing (FWHM = 6\,mm) on tSNR for Dataset 2 ' ...
         '(ds000105). Layout identical to Fig.~\ref{fig:ds1:smoothing}.'], ...
        'ds2:smoothing');
end

%% -----------------------------------------------------------------------
%  STEP 3.13 — VALIDATION
%  -----------------------------------------------------------------------
fprintf('\n--- [3.13] Validation...\n');
all_pass = true;
checks = {
    swra_func, ['Smoothed output: swra' fn_bold '.nii']
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
V_final = spm_vol(swra_func);
if numel(V_final)==num_vols
    fprintf('  PASS: Volume count = %d\n', numel(V_final));
else
    fprintf('  WARN: Expected %d vols, got %d\n', num_vols, numel(V_final));
end
vs = sqrt(sum(V_final(1).mat(1:3,1:3).^2));
fprintf('  INFO: Voxel size = [%.2f %.2f %.2f] mm\n', vs(1),vs(2),vs(3));
max_t_mm  = max(abs(rp(:,1:3)), [],'all');
max_r_deg = max(abs(rp(:,4:6))*180/pi, [],'all');
fd = sum(abs(diff(rp(:,1:3))),2) + sum(abs(diff(rp(:,4:6)))*50,2);
fprintf('  INFO: Max translation=%.2f mm | Max rotation=%.2f deg\n', max_t_mm, max_r_deg);
fprintf('  INFO: Mean FD = %.2f mm\n', mean(fd));

status_str = iif(all_pass,'PASS','WARNINGS');
fprintf('\n=== PREPROCESSING COMPLETE — Dataset 2 | Status: %s ===\n\n', status_str);

%% -----------------------------------------------------------------------
%  STEP 3.14 — SAVE RESULTS
%  -----------------------------------------------------------------------
fprintf('--- [3.14] Saving results...\n');
dest_func = fullfile(OUT_DIR,'swra_bold_ds2.nii');
if ~exist(dest_func,'file') || dir(dest_func).bytes < dir(swra_func).bytes
    copyfile(swra_func, dest_func);
end
if exist(rp_func,'file'), copyfile(rp_func, fullfile(OUT_DIR,'rp_bold_ds2.txt')); end

preproc = struct();
preproc.dataset     = DATASET;
preproc.subject     = SUBJECT;
preproc.TR          = TR;
preproc.num_slices  = NUM_SLICES;
preproc.num_volumes = num_vols;
preproc.smooth_fwhm = SMOOTH_FWHM;
preproc.smoothed_file = swra_func;
preproc.max_trans_mm  = max_t_mm;
preproc.max_rot_deg   = max_r_deg;
preproc.mean_fd_mm    = mean(fd);
preproc.status        = status_str;
preproc.timestamp     = datestr(now);
save(fullfile(OUT_DIR,'result_ds2_preprocessing_summary.mat'),'-struct','preproc');
fprintf('  Saved: result_ds2_preprocessing_summary.mat\n');

csv_file = fullfile(OUT_DIR,'table_ds2_motion_parameters.csv');
fid = fopen(csv_file,'w');
fprintf(fid,'Volume,TransX_mm,TransY_mm,TransZ_mm,RotPitch_deg,RotRoll_deg,RotYaw_deg,FD_mm\n');
fd_all = [0; fd];
for v = 1:size(rp,1)
    fprintf(fid,'%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n',...
        v, rp(v,1),rp(v,2),rp(v,3), ...
        rp(v,4)*180/pi, rp(v,5)*180/pi, rp(v,6)*180/pi, fd_all(v));
end
fclose(fid);
fprintf('  Saved: table_ds2_motion_parameters.csv\n');

%% PROJECT LOG
fid = fopen(LOG_FILE,'a');
fprintf(fid,'========================================================\n');
fprintf(fid,'STEP 03 : Dataset 2 Preprocessing\n');
fprintf(fid,'Date    : %s\n', datestr(now));
fprintf(fid,'Dataset : %s | Subject: %s\n', DATASET, SUBJECT);
fprintf(fid,'TR      : %.2f s | Slices: %d | Volumes: %d\n', TR, NUM_SLICES, num_vols);
fprintf(fid,'Output  : %s\n', swra_func);
fprintf(fid,'Motion  : max=%.2f mm, max=%.2f deg, meanFD=%.2f\n', ...
    max_t_mm, max_r_deg, mean(fd));
fprintf(fid,'STATUS  : %s\n', status_str);
fprintf(fid,'========================================================\n\n');
fclose(fid);

fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', fullfile(ROOT_DIR,'github'));
fprintf('git add .\n');
fprintf('git commit -m "Step 03 complete: Dataset 2 preprocessing (%s)"\n', status_str);
fprintf('git push origin main\n');
fprintf('\n=== END OF STEP 03 ===\n');

%% -----------------------------------------------------------------------
%  LOCAL HELPER FUNCTIONS
%  -----------------------------------------------------------------------

function tsnr = compute_tsnr(V)
    nv = numel(V); d0 = spm_read_vols(V(1));
    [nx,ny,nz] = size(d0);
    all_vols = zeros(nx,ny,nz,nv,'single');
    for k = 1:nv, all_vols(:,:,:,k) = single(spm_read_vols(V(k))); end
    mu = mean(all_vols,4); sg = std(all_vols,0,4);
    sg(sg<single(eps)) = single(eps);
    tsnr = double(mu./sg);
end

function saveFig(fig, fig_dir, name)
    print(fig, fullfile(fig_dir,[name '.png']), '-dpng','-r150');
    print(fig, fullfile(fig_dir,[name '.eps']), '-depsc');
    close(fig);
    fprintf('  Saved: %s.png / .eps\n', name);
end

function out = norm01(img)
    mn = min(img(:)); mx = max(img(:));
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
