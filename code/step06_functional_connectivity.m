%% =========================================================================
%  FILE        : step06_functional_connectivity.m
%  TITLE       : fMRI Brain Data Analysis using SPM and MATLAB —
%                Preprocessing, GLM, Statistical Mapping and
%                Functional Connectivity
%  AUTHOR      : Fahim Ur Rehman Shah
%  REG NO      : EE2629
%  COURSE      : CSE532 — Signal and Image Processing (MS Level)
%  SUPERVISOR  : Dr. Adnan Shah, FCSE, GIKI
%  INSTITUTION : GIK Institute of Engineering Sciences & Technology
%  DATE        : May 2026
%  STEP        : 06 of 07
%
%  DESCRIPTION : Seed-based functional connectivity analysis for both
%                datasets. For Dataset 1 (motor task) the seed is placed
%                in primary motor cortex (M1). For Dataset 2 (visual task)
%                the seed is in primary visual cortex (V1).
%                Method: extract mean seed time series → correlate with
%                every voxel → Fisher z-transform → save connectivity maps.
%
%  INPUTS      : results/dataset1/preprocessing/swra_bold_ds1.nii
%                results/dataset2/preprocessing/swra_bold_ds2.nii
%                (or the preprocessed files in the data/ directories)
%  OUTPUTS     : results/dataset1/connectivity/conn_map_ds1.nii
%                results/dataset2/connectivity/conn_map_ds2.nii
%                figures/fig21–fig23_connectivity_*.png/eps
%                results/*/connectivity/table_*_connectivity_peaks.csv
%  DEPENDENCIES: SPM25, MATLAB R2022b+, steps 02–03 complete
%  =========================================================================

clear; clc; close all;
fprintf('=== Step 06: Functional Connectivity (Both Datasets) ===\n\n');

%% -----------------------------------------------------------------------
%  PARAMETERS
%  -----------------------------------------------------------------------
ROOT_DIR = 'C:\fmri_project';
TR       = 2.5;

% Seed regions in MNI coordinates (mm)
% DS1: Primary Motor Cortex (M1), hand area
SEED_DS1_MNI = [-38, -26, 60];   % L M1 hand [Biswal 1995 analog]
SEED_RADIUS  = 6;                  % mm radius sphere

% DS2: Primary Visual Cortex (V1), calcarine
SEED_DS2_MNI = [0, -88, 2];       % bilateral V1 calcarine

FIG_DIR  = fullfile(ROOT_DIR,'figures');
LOG_FILE = fullfile(ROOT_DIR,'project_log.txt');
if ~exist(FIG_DIR,'dir'), mkdir(FIG_DIR); end

%% -----------------------------------------------------------------------
%  DATASET CONFIGURATIONS
%  -----------------------------------------------------------------------
ds(1).name       = 'ds000114';
ds(1).subject    = 'sub-01';
ds(1).seed_mni   = SEED_DS1_MNI;
ds(1).seed_label = 'M1_hand';
ds(1).func_dir   = fullfile(ROOT_DIR,'data','ds000114','func');
ds(1).fn_bold    = 'sub-01_ses-test_task-fingerfootlips_run-1_bold';
ds(1).conn_dir   = fullfile(ROOT_DIR,'results','dataset1','connectivity');
ds(1).fig_prefix = 'ds1';
ds(1).fig_offset = 21;    % fig21 for DS1 connectivity

ds(2).name       = 'ds000105';
ds(2).subject    = 'sub-1';
ds(2).seed_mni   = SEED_DS2_MNI;
ds(2).seed_label = 'V1_calcarine';
ds(2).func_dir   = fullfile(ROOT_DIR,'data','ds000105','func');
ds(2).fn_bold    = 'sub-1_task-objectviewing_run-1_bold';
ds(2).conn_dir   = fullfile(ROOT_DIR,'results','dataset2','connectivity');
ds(2).fig_prefix = 'ds2';
ds(2).fig_offset = 22;    % fig22 for DS2 connectivity

spm('defaults','FMRI');
spm_jobman('initcfg');

conn_results = struct();

for d = 1:2
    fprintf('\n--- Dataset %d: %s (seed: %s) ---\n', d, ds(d).name, ds(d).seed_label);
    if ~exist(ds(d).conn_dir,'dir'), mkdir(ds(d).conn_dir); end

    % Locate smoothed functional
    swra_func = fullfile(ds(d).func_dir, ['swra' ds(d).fn_bold '.nii']);
    if ~exist(swra_func,'file')
        % Try results folder
        swra_func2 = fullfile(ROOT_DIR,'results', ...
            sprintf('dataset%d',d),'preprocessing', ...
            sprintf('swra_bold_ds%d.nii',d));
        if exist(swra_func2,'file')
            swra_func = swra_func2;
        else
            fprintf('  ERROR: Smoothed file not found: %s\n', swra_func);
            fprintf('  Run step0%d_preprocessing first.\n', d+1);
            continue;
        end
    end
    fprintf('  Functional: %s\n', swra_func);

    %% -------------------------------------------------------------------
    %  6.A — LOAD DATA AND EXTRACT SEED TIME SERIES
    %  -------------------------------------------------------------------
    fprintf('  Loading data...\n');
    V = spm_vol(swra_func);
    num_vols = numel(V);
    ref_vol  = V(1);
    dims     = ref_vol.dim;

    % Build seed mask (sphere in MNI space)
    seed_mni   = ds(d).seed_mni(:)';   % 1x3
    seed_mask  = make_sphere_mask(ref_vol, seed_mni, SEED_RADIUS);
    n_seed_vox = sum(seed_mask(:));
    fprintf('  Seed region: %d voxels (radius=%dmm, MNI=[%g,%g,%g])\n', ...
        n_seed_vox, SEED_RADIUS, seed_mni(1), seed_mni(2), seed_mni(3));

    % Load all volumes into memory (4D matrix)
    fprintf('  Loading %d volumes (may take ~30s)...\n', num_vols);
    data4d = zeros([dims, num_vols], 'single');
    for v = 1:num_vols
        data4d(:,:,:,v) = single(spm_read_vols(V(v)));
    end
    data2d = reshape(data4d, [], num_vols);  % [voxels x time]
    clear data4d;
    fprintf('  Data loaded: %d voxels x %d timepoints\n', size(data2d,1), size(data2d,2));

    % Extract seed time series (mean across seed voxels)
    seed_flat  = reshape(seed_mask, [], 1);
    seed_ts    = mean(data2d(seed_flat, :), 1)';  % [num_vols x 1]

    %% -------------------------------------------------------------------
    %  6.B — BAND-PASS FILTER (optional low-level)
    %  -------------------------------------------------------------------
    % High-pass filter (0.01 Hz) — remove slow drifts
    % We use SPM's discrete cosine transform approach
    fprintf('  High-pass filtering seed time series (0.01 Hz)...\n');
    K.HParam = 100;  % 100s high-pass cutoff
    K.RT     = TR;
    K.row    = 1:num_vols;
    K = spm_filter(K);
    seed_ts_hp = spm_filter(K, seed_ts);

    % Also filter the whole-brain data
    data_hp = spm_filter(K, data2d');  % spm_filter expects [time x voxels]
    data_hp = data_hp';                % back to [voxels x time]

    %% -------------------------------------------------------------------
    %  6.C — COMPUTE WHOLE-BRAIN CORRELATION
    %  -------------------------------------------------------------------
    fprintf('  Computing voxel-wise correlations...\n');

    % Demean
    seed_z = seed_ts_hp - mean(seed_ts_hp);
    seed_z = seed_z / (std(seed_z) + eps);

    data_dm = data_hp - mean(data_hp, 2);           % [vox x time] demean
    data_std = std(data_hp, 0, 2) + eps;            % [vox x 1]
    data_norm = data_dm ./ data_std;

    % Pearson r for all voxels simultaneously
    r_map = (data_norm * seed_z) / (num_vols - 1);  % [vox x 1]

    % Fisher z-transform
    r_map = min(max(r_map, -0.9999), 0.9999);       % clip for atanh
    z_map = atanh(r_map);

    % Reshape to volume
    r_vol = reshape(r_map, dims);
    z_vol = reshape(z_map, dims);

    %% -------------------------------------------------------------------
    %  6.D — SAVE CONNECTIVITY MAPS
    %  -------------------------------------------------------------------
    % r-map
    r_fname = fullfile(ds(d).conn_dir, ...
        sprintf('conn_r_%s_%s.nii', ds(d).fig_prefix, ds(d).seed_label));
    Vout = ref_vol;
    Vout.fname  = r_fname;
    Vout.dt     = [spm_type('float32') spm_platform('bigend')];
    Vout.dim    = dims;
    Vout.pinfo  = [1;0;0];
    spm_write_vol(Vout, r_vol);

    % z-map
    z_fname = fullfile(ds(d).conn_dir, ...
        sprintf('conn_z_%s_%s.nii', ds(d).fig_prefix, ds(d).seed_label));
    Vout.fname = z_fname;
    spm_write_vol(Vout, z_vol);

    fprintf('  Saved r-map: %s\n', r_fname);
    fprintf('  Saved z-map: %s\n', z_fname);

    % Store for comparison later
    conn_results(d).z_vol      = z_vol;
    conn_results(d).r_vol      = r_vol;
    conn_results(d).dims       = dims;
    conn_results(d).mat        = ref_vol.mat;
    conn_results(d).seed_mni   = seed_mni;
    conn_results(d).seed_mask  = seed_mask;
    conn_results(d).name       = ds(d).name;
    conn_results(d).seed_label = ds(d).seed_label;
    conn_results(d).num_vols   = num_vols;
    conn_results(d).z_fname    = z_fname;

    %% -------------------------------------------------------------------
    %  FIGURE — CONNECTIVITY MAP
    %  -------------------------------------------------------------------
    fig_name = sprintf('fig%02d_%s_connectivity_%s', ...
        ds(d).fig_offset, ds(d).fig_prefix, ds(d).seed_label);

    if ~exist(fullfile(FIG_DIR,[fig_name '.png']),'file')
        fprintf('  Generating connectivity map figure...\n');
        R_THRESH = 0.3;  % display threshold
        z_thresh = atanh(R_THRESH);

        fig = figure('Visible','off','Position',[100 100 1400 500]);

        % Show 4 axial slices at evenly spaced z-levels
        nz = dims(3);
        sl_list = round(linspace(nz*0.25, nz*0.75, 4));

        for p = 1:4
            sl = max(1, min(sl_list(p), nz));
            subplot(1,4,p);

            z_sl = z_vol(:,:,sl)';
            bg   = sum(data_dm(:,1:min(10,num_vols)),2)/(min(10,num_vols));
            bg   = reshape(bg, dims); bg_sl = bg(:,:,sl)';
            bg_n = (bg_sl - min(bg_sl(:)))/(max(bg_sl(:))-min(bg_sl(:))+eps);

            % Positive connectivity (warm colours)
            pos_mask = z_sl > z_thresh;
            % Negative connectivity (cool colours)
            neg_mask = z_sl < -z_thresh;

            rgb = repmat(bg_n, [1 1 3]);
            % Add warm (positive)
            rgb(:,:,1) = min(rgb(:,:,1) + pos_mask * 0.9, 1);
            rgb(:,:,2) = min(rgb(:,:,2) + pos_mask * 0.4, 1);
            rgb(:,:,3) = max(rgb(:,:,3) - pos_mask * 0.5, 0);
            % Add cool (negative)
            rgb(:,:,1) = max(rgb(:,:,1) - neg_mask * 0.5, 0);
            rgb(:,:,2) = min(rgb(:,:,2) + neg_mask * 0.4, 1);
            rgb(:,:,3) = min(rgb(:,:,3) + neg_mask * 0.9, 1);

            imagesc(rgb); axis image off;
            title(sprintf('z=%d', sl), 'FontSize',9);
        end

        sgtitle(sprintf('Connectivity Map — %s | Seed: %s (MNI [%g,%g,%g]) | r>%.1f', ...
            ds(d).name, ds(d).seed_label, ...
            seed_mni(1), seed_mni(2), seed_mni(3), R_THRESH), 'FontSize',10);

        saveFig(fig, FIG_DIR, fig_name);
        printLatex(fig_name, ...
            sprintf(['Seed-based functional connectivity map for Dataset %d (%s). ' ...
            'Seed region: %s (MNI [%g, %g, %g], %d mm radius). ' ...
            'Warm colours: positive connectivity ($r > %.1f$). ' ...
            'Cool colours: negative connectivity ($r < -%.1f$).'], ...
            d, ds(d).name, strrep(ds(d).seed_label,'_',' '), ...
            seed_mni(1), seed_mni(2), seed_mni(3), SEED_RADIUS, R_THRESH, R_THRESH), ...
            sprintf('%s:connectivity', ds(d).fig_prefix));
    end

    %% Save peak connectivity table
    csv_file = fullfile(ds(d).conn_dir, ...
        sprintf('table_%s_connectivity_peaks.csv', ds(d).fig_prefix));
    r_thresh_csv = 0.3;
    [~, top_idx] = sort(r_map, 'descend');
    fid = fopen(csv_file,'w');
    fprintf(fid, 'Rank,r_value,z_value,MNI_x,MNI_y,MNI_z\n');
    count = 0;
    for i = 1:numel(top_idx)
        if r_map(top_idx(i)) < r_thresh_csv, break; end
        [vx,vy,vz] = ind2sub(dims, top_idx(i));
        mni_c = ref_vol.mat * [vx;vy;vz;1];
        count = count + 1;
        fprintf(fid,'%d,%.4f,%.4f,%.1f,%.1f,%.1f\n', ...
            count, r_map(top_idx(i)), z_map(top_idx(i)), ...
            mni_c(1), mni_c(2), mni_c(3));
        if count >= 20, break; end
    end
    fclose(fid);
    fprintf('  Saved connectivity peaks: %s\n', csv_file);

    %% Log step
    fid = fopen(LOG_FILE,'a');
    fprintf(fid,'STEP 06 [DS%d]: %s | seed=%s | nvox=%d | time=%s\n', ...
        d, ds(d).name, ds(d).seed_label, n_seed_vox, datestr(now));
    fclose(fid);
end

%% -----------------------------------------------------------------------
%  FIGURE 23 — SIDE-BY-SIDE CONNECTIVITY COMPARISON
%  -----------------------------------------------------------------------
fig23_name = 'fig23_ds1_vs_ds2_connectivity_comparison';
if ~exist(fullfile(FIG_DIR,[fig23_name '.png']),'file') && ...
   isfield(conn_results,'z_vol') && numel(conn_results) == 2 && ...
   ~isempty(conn_results(1).z_vol) && ~isempty(conn_results(2).z_vol)

    fprintf('\n--- Generating cross-dataset connectivity comparison figure...\n');
    Z_THRESH_DISP = atanh(0.3);

    fig = figure('Visible','off','Position',[100 100 1400 600]);
    titles = {'DS1 (Motor — M1 seed)', 'DS2 (Visual — V1 seed)'};

    for d = 1:2
        zv  = conn_results(d).z_vol;
        nz  = conn_results(d).dims(3);
        sl  = round(nz * 0.50);

        z_sl = zv(:,:,sl)';
        pos  = z_sl > Z_THRESH_DISP;
        neg  = z_sl < -Z_THRESH_DISP;

        bg_n = (z_sl - min(z_sl(:))) / (max(z_sl(:)) - min(z_sl(:)) + eps);
        rgb  = repmat(bg_n,[1 1 3]);
        rgb(:,:,1) = min(rgb(:,:,1) + pos*0.9, 1);
        rgb(:,:,2) = min(rgb(:,:,2) + pos*0.4, 1);
        rgb(:,:,3) = max(rgb(:,:,3) - pos*0.5, 0);
        rgb(:,:,1) = max(rgb(:,:,1) - neg*0.5, 0);
        rgb(:,:,2) = min(rgb(:,:,2) + neg*0.4, 1);
        rgb(:,:,3) = min(rgb(:,:,3) + neg*0.9, 1);

        subplot(2,4, (d-1)*4 + (1:4));
        imagesc(rgb); axis image off;
        title(titles{d}, 'FontSize',11);
    end
    sgtitle(sprintf('Cross-Dataset Connectivity Comparison (r > 0.30)'), 'FontSize',12);
    saveFig(fig, FIG_DIR, fig23_name);
    printLatex(fig23_name, ...
        ['Side-by-side comparison of seed-based connectivity maps for Dataset 1 ' ...
         '(motor cortex seed, top) and Dataset 2 (visual cortex seed, bottom). ' ...
         'Warm: positive connectivity; cool: negative. Threshold $r > 0.30$.'], ...
        'comparison:connectivity');
end

%% -----------------------------------------------------------------------
%  SAVE SUMMARY
%  -----------------------------------------------------------------------
fprintf('\n--- Saving combined connectivity summary...\n');
summary_file = fullfile(ROOT_DIR,'results','comparison','result_connectivity_summary.mat');
if ~exist(fullfile(ROOT_DIR,'results','comparison'),'dir')
    mkdir(fullfile(ROOT_DIR,'results','comparison'));
end
save(summary_file,'conn_results','SEED_DS1_MNI','SEED_DS2_MNI','SEED_RADIUS');
fprintf('  Saved: %s\n', summary_file);

%% -----------------------------------------------------------------------
%  VALIDATION
%  -----------------------------------------------------------------------
fprintf('\n--- [6.Z] Validation...\n');
all_pass = true;
for d = 1:2
    fname = fullfile(ds(d).conn_dir, ...
        sprintf('conn_z_%s_%s.nii', ds(d).fig_prefix, ds(d).seed_label));
    if exist(fname,'file')
        info = dir(fname);
        fprintf('  PASS: DS%d z-map (%.1f MB)\n', d, info.bytes/1e6);
    else
        fprintf('  FAIL: DS%d z-map not found\n', d);
        all_pass = false;
    end
end

status_str = iif(all_pass,'PASS','WARNINGS');
fprintf('\n=== CONNECTIVITY COMPLETE | Status: %s ===\n\n', status_str);

fid = fopen(LOG_FILE,'a');
fprintf(fid,'========================================================\n');
fprintf(fid,'STEP 06 : Functional Connectivity\n');
fprintf(fid,'Date    : %s\n', datestr(now));
fprintf(fid,'DS1 seed: %s at MNI [%g,%g,%g]\n', ...
    ds(1).seed_label, ds(1).seed_mni(1), ds(1).seed_mni(2), ds(1).seed_mni(3));
fprintf(fid,'DS2 seed: %s at MNI [%g,%g,%g]\n', ...
    ds(2).seed_label, ds(2).seed_mni(1), ds(2).seed_mni(2), ds(2).seed_mni(3));
fprintf(fid,'STATUS  : %s\n', status_str);
fprintf(fid,'========================================================\n\n');
fclose(fid);

fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', fullfile(ROOT_DIR,'github'));
fprintf('git add .\n');
fprintf('git commit -m "Step 06 complete: Functional connectivity maps for both datasets"\n');
fprintf('git push origin main\n');
fprintf('\n=== END OF STEP 06 ===\n');

%% -----------------------------------------------------------------------
%  LOCAL HELPER FUNCTIONS
%  -----------------------------------------------------------------------

function mask = make_sphere_mask(Vref, mni_centre, radius_mm)
    % Build a binary sphere mask in voxel space given an MNI centre
    dims = Vref.dim;
    [xi,yi,zi] = ndgrid(1:dims(1), 1:dims(2), 1:dims(3));
    vox_coords  = [xi(:)';  yi(:)';  zi(:)';  ones(1,numel(xi))];
    mni_all     = Vref.mat * vox_coords;  % 4 x nvox
    d2 = sum((mni_all(1:3,:) - mni_centre(:)).^2, 1);
    mask = reshape(d2 <= radius_mm^2, dims);
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
