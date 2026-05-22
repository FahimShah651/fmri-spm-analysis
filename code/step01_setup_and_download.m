%% =========================================================================
%  SETUP AND DATA DOWNLOAD — fMRI Brain Data Analysis using SPM and MATLAB
%  =========================================================================
%  Title      : Preprocessing, GLM, Statistical Mapping and Functional
%               Connectivity
%  Author     : Fahim Ur Rehman Shah
%  Reg No     : EE2629
%  Course     : CSE532 — Signal and Image Processing (MS Level)
%  Supervisor : Dr. Adnan Shah, FCSE, GIKI
%  Institution: GIK Institute of Engineering Sciences & Technology
%  Date       : May 2026
%
%  STEP NUMBER : 01
%  DESCRIPTION : Check environment (MATLAB, SPM, Git), create full
%                directory structure, verify datasets exist or download
%                them, initialize git repo, create README, log everything.
%
%  INPUTS     : None (datasets downloaded from OpenNeuro S3)
%  OUTPUTS    : Verified directory tree, project_log.txt, README.md
%  DEPENDENCIES: SPM25, MATLAB R2022b+, Git, Python 3 (for downloads)
%  =========================================================================

clear; clc; close all;
fprintf('=== fMRI Project Setup — Step 01 ===\n\n');

%% ------------------------------------------------------------------------
%  PARAMETERS
%  ------------------------------------------------------------------------
ROOT_DIR    = 'C:\fmri_project';
DATASET_1   = 'ds000114';
DATASET_2   = 'ds000105';
LOG_FILE    = fullfile(ROOT_DIR, 'project_log.txt');
README_FILE = fullfile(ROOT_DIR, 'README.md');

%% ------------------------------------------------------------------------
%  1. CHECK MATLAB VERSION
%  ------------------------------------------------------------------------
fprintf('---[1] Checking MATLAB Version...\n');
v = ver('MATLAB');
matlab_ok = str2double(v.Version) >= 9.3;
if matlab_ok
    fprintf('  PASS: MATLAB %s (%s) detected.\n', v.Version, v.Release);
else
    fprintf('  WARNING: MATLAB %s (%s) found. R2022b+ recommended.\n', ...
            v.Version, v.Release);
end

%% ------------------------------------------------------------------------
%  2. CHECK SPM INSTALLATION
%  ------------------------------------------------------------------------
fprintf('---[2] Checking SPM Installation...\n');
spm_ok = false;
try
    spm_path = fileparts(which('spm'));
    if ~isempty(spm_path)
        spm_ver = spm('Ver');
        fprintf('  PASS: %s found at:\n    %s\n', spm_ver, spm_path);
        spm_ok = true;
    else
        fprintf('  FAIL: spm() not on MATLAB path.\n');
    end
catch ME
    fprintf('  FAIL: %s\n', ME.message);
end

%% ------------------------------------------------------------------------
%  3. CHECK GIT INSTALLATION
%  ------------------------------------------------------------------------
fprintf('---[3] Checking Git Installation...\n');
[git_ok, git_out] = system('git --version');
if git_ok == 0
    fprintf('  PASS: %s', strtrim(git_out));
else
    fprintf('  WARNING: Git not found. Install from https://git-scm.com\n');
end

%% ------------------------------------------------------------------------
%  4. CREATE DIRECTORY TREE
%  ------------------------------------------------------------------------
fprintf('---[4] Creating Directory Structure...\n');
dirs = {
    'code\utils'
    'data\ds000114\func'
    'data\ds000114\anat'
    'data\ds000105\func'
    'data\ds000105\anat'
    'figures'
    'github'
    'latex\figures'
    'latex\sections'
    'results\dataset1\preprocessing'
    'results\dataset1\glm'
    'results\dataset1\connectivity'
    'results\dataset2\preprocessing'
    'results\dataset2\glm'
    'results\dataset2\connectivity'
    'results\comparison'
};
for i = 1:numel(dirs)
    d = fullfile(ROOT_DIR, dirs{i});
    if ~exist(d, 'dir'); mkdir(d); fprintf('  Created: %s\n', d);
    else; fprintf('  Exists : %s\n', d); end
end

%% ------------------------------------------------------------------------
%  5. VERIFY DATASET 1
%  ------------------------------------------------------------------------
fprintf('---[5] Verifying Dataset 1 (ds000114 — Motor Task)...\n');
ds1_anat  = fullfile(ROOT_DIR, 'data', DATASET_1, 'anat');
ds1_func  = fullfile(ROOT_DIR, 'data', DATASET_1, 'func');

ds1_files = {...
    'sub-01_ses-test_T1w.nii.gz' ...
    'sub-01_ses-test_task-fingerfootlips_run-1_bold.nii.gz' ...
    'sub-01_ses-test_task-fingerfootlips_run-1_events.tsv'};

ds1_ok = true;
for f = 1:numel(ds1_files)
    fullp = fullfile(ds1_func, ds1_files{f});
    if ~exist(fullp, 'file'); fullp = fullfile(ds1_anat, ds1_files{f}); end
    if exist(fullp, 'file')
        info = dir(fullp);
        if info.bytes > 1e6
            fprintf('  [%s] %.2f MB\n', ds1_files{f}, info.bytes/1e6);
        elseif info.bytes > 1000
            fprintf('  [%s] %.2f KB\n', ds1_files{f}, info.bytes/1e3);
        else
            fprintf('  [%s] %d bytes\n', ds1_files{f}, info.bytes);
        end
    else
        fprintf('  [%s] MISSING\n', ds1_files{f});
        ds1_ok = false;
    end
end

%% ------------------------------------------------------------------------
%  6. VERIFY DATASET 2
%  ------------------------------------------------------------------------
fprintf('---[6] Verifying Dataset 2 (ds000105 — Visual)...\n');
ds2_anat  = fullfile(ROOT_DIR, 'data', DATASET_2, 'anat');
ds2_func  = fullfile(ROOT_DIR, 'data', DATASET_2, 'func');

ds2_files = {...
    'sub-1_T1w.nii.gz' ...
    'sub-1_task-objectviewing_run-1_bold.nii.gz' ...
    'sub-1_task-objectviewing_run-1_events.tsv'};

ds2_ok = true;
for f = 1:numel(ds2_files)
    fullp = fullfile(ds2_func, ds2_files{f});
    if ~exist(fullp, 'file'); fullp = fullfile(ds2_anat, ds2_files{f}); end
    if exist(fullp, 'file')
        info = dir(fullp);
        if info.bytes > 1e6
            fprintf('  [%s] %.2f MB\n', ds2_files{f}, info.bytes/1e6);
        elseif info.bytes > 1000
            fprintf('  [%s] %.2f KB\n', ds2_files{f}, info.bytes/1e3);
        else
            fprintf('  [%s] %d bytes\n', ds2_files{f}, info.bytes);
        end
    else
        fprintf('  [%s] MISSING\n', ds2_files{f});
        ds2_ok = false;
    end
end

%% ------------------------------------------------------------------------
%  7. DOWNLOAD DATASET FILES (if missing or empty)
%  ------------------------------------------------------------------------
fprintf('---[7] Downloading Dataset Files...\n');

% Use the standalone Python download script (robust S3 downloader)
py_script = fullfile(ROOT_DIR, 'code', 'utils', 'download_datasets.py');
fprintf('  Calling: python %s\n', py_script);
[st, out] = system(sprintf('python "%s" "%s"', py_script, ROOT_DIR));
fprintf('%s\n', strtrim(out));

% Re-verify Dataset 1 (NIfTI > 10KB, TSV > 30 bytes)
fprintf('  Re-verifying datasets...\n');
ds1_ok = true;
for f = 1:numel(ds1_files)
    dest = fullfile(ds1_func, ds1_files{f});
    if ~exist(dest,'file'); dest = fullfile(ds1_anat, ds1_files{f}); end
    min_bytes = iif(endsWith(ds1_files{f}, '.nii.gz'), 10000, 30);
    if ~exist(dest,'file') || dir(dest).bytes < min_bytes
        ds1_ok = false;
        fprintf('  Dataset 1: %s MISSING or empty\n', ds1_files{f});
    else
        info = dir(dest);
        if info.bytes > 1e6
            fprintf('  Dataset 1: %s OK (%.2f MB)\n', ds1_files{f}, info.bytes/1e6);
        else
            fprintf('  Dataset 1: %s OK (%d bytes)\n', ds1_files{f}, info.bytes);
        end
    end
end

ds2_ok = true;
for f = 1:numel(ds2_files)
    dest = fullfile(ds2_func, ds2_files{f});
    if ~exist(dest,'file'); dest = fullfile(ds2_anat, ds2_files{f}); end
    min_bytes = iif(endsWith(ds2_files{f}, '.nii.gz'), 10000, 30);
    if ~exist(dest,'file') || dir(dest).bytes < min_bytes
        ds2_ok = false;
        fprintf('  Dataset 2: %s MISSING or empty\n', ds2_files{f});
    else
        info = dir(dest);
        if info.bytes > 1e6
            fprintf('  Dataset 2: %s OK (%.2f MB)\n', ds2_files{f}, info.bytes/1e6);
        else
            fprintf('  Dataset 2: %s OK (%d bytes)\n', ds2_files{f}, info.bytes);
        end
    end
end

%% ------------------------------------------------------------------------
%  8. INITIALIZE GIT REPOSITORY
%  ------------------------------------------------------------------------
fprintf('---[8] Initializing Git Repository...\n');
git_dir = fullfile(ROOT_DIR, 'github');
if ~exist(fullfile(git_dir, '.git'), 'dir')
    old = cd(git_dir);
    system('git init');
    fprintf('  Git repository initialized at: %s\n', git_dir);
    cd(old);
else
    fprintf('  Git repository already exists.\n');
end

%% ------------------------------------------------------------------------
%  9. CREATE .gitignore
%  ------------------------------------------------------------------------
fprintf('---[9] Creating .gitignore...\n');
gitignore_path = fullfile(ROOT_DIR, '.gitignore');
if ~exist(gitignore_path, 'file')
    fid = fopen(gitignore_path, 'w');
    fprintf(fid, '%% MATLAB\n*.mat\n*.asv\n*~\n\n');
    fprintf(fid, '%% Data (large NIfTI files)\n');
    fprintf(fid, 'data/*.nii\ndata/*.nii.gz\ndata/*.img\ndata/*.hdr\n\n');
    fprintf(fid, '%% Results (regenerated)\nresults/*.mat\nresults/*.nii\n\n');
    fprintf(fid, '%% System\nThumbs.db\n.DS_Store\n');
    fclose(fid);
    fprintf('  Created .gitignore\n');
else
    fprintf('  .gitignore already exists.\n');
end

%% ------------------------------------------------------------------------
%  10. CREATE README.md
%  ------------------------------------------------------------------------
fprintf('---[10] Creating README.md...\n');
if ~exist(README_FILE, 'file')
    fid = fopen(README_FILE, 'w');
    fprintf(fid, '# fMRI Brain Data Analysis using SPM and MATLAB\n\n');
    fprintf(fid, '**Author:** Fahim Ur Rehman Shah (EE2629)\n\n');
    fprintf(fid, '**Course:** CSE532 — Signal and Image Processing (MS Level)\n\n');
    fprintf(fid, '**Supervisor:** Dr. Adnan Shah, FCSE, GIKI\n\n');
    fprintf(fid, '---\n\n## Overview\n\n');
    fprintf(fid, 'End-to-end fMRI analysis pipeline:\n');
    fprintf(fid, '- Preprocessing (SliceTiming, Realign, Coreg, Segment, Norm, Smooth)\n');
    fprintf(fid, '- GLM estimation and contrast mapping\n');
    fprintf(fid, '- Seed-based functional connectivity\n');
    fprintf(fid, '- Cross-dataset comparative analysis\n\n');
    fprintf(fid, '## Pipeline Steps\n\n');
    fprintf(fid, '| Step | Script | Description |\n');
    fprintf(fid, '|------|--------|-------------|\n');
    fprintf(fid, '| 01 | step01_* | Setup and data download |\n');
    fprintf(fid, '| 02 | step02_* | Preprocessing — Dataset 1 |\n');
    fprintf(fid, '| 03 | step03_* | Preprocessing — Dataset 2 |\n');
    fprintf(fid, '| 04 | step04_* | GLM — Dataset 1 |\n');
    fprintf(fid, '| 05 | step05_* | GLM — Dataset 2 |\n');
    fprintf(fid, '| 06 | step06_* | Connectivity — Both |\n');
    fprintf(fid, '| 07 | step07_* | Cross-dataset comparison |\n\n');
    fprintf(fid, '---\n*Generated: %s*\n', datestr(now));
    fclose(fid);
    fprintf('  Created README.md\n');
else
    fprintf('  README.md already exists.\n');
end

%% ------------------------------------------------------------------------
%  11. UPDATE PROJECT LOG
%  ------------------------------------------------------------------------
fprintf('---[11] Updating Project Log...\n');
fid = fopen(LOG_FILE, 'a');
fprintf(fid, '========================================================\n');
fprintf(fid, 'STEP 01: Setup and Data Verification\n');
fprintf(fid, 'Date/Time: %s\n', datestr(now));
fprintf(fid, 'MATLAB    : %s (%s)\n', v.Version, v.Release);
if spm_ok; fprintf(fid, 'SPM       : %s\n', spm_ver);
else; fprintf(fid, 'SPM       : NOT FOUND\n'); end
fprintf(fid, 'Dataset 1 : %s\n', iif(ds1_ok,'PASS','MISSING FILES'));
fprintf(fid, 'Dataset 2 : %s\n', iif(ds2_ok,'PASS','MISSING FILES'));
fprintf(fid, 'STATUS    : %s\n', iif(ds1_ok&&ds2_ok,'PASS','WARNING'));
fprintf(fid, '========================================================\n\n');
fclose(fid);
fprintf('  Appended to project_log.txt\n');

%% ------------------------------------------------------------------------
%  12. SUMMARY
%  ------------------------------------------------------------------------
fprintf('\n=============== SETUP COMPLETE ===============\n');
fprintf('  MATLAB   : %s\n', iif(matlab_ok, 'OK', 'WARN'));
fprintf('  SPM      : %s\n', iif(spm_ok, 'OK', 'FAIL'));
fprintf('  Git      : %s\n', iif(git_ok==0, 'OK', 'WARN'));
fprintf('  Dataset1 : %s\n', iif(ds1_ok, 'OK', 'FAIL'));
fprintf('  Dataset2 : %s\n', iif(ds2_ok, 'OK', 'FAIL'));
fprintf('=============================================\n');

%% ------------------------------------------------------------------------
%  HELPER FUNCTION
%  ------------------------------------------------------------------------
function s = iif(cond, tval, fval)
    if cond; s = tval; else; s = fval; end
end

%% =========================================================================
%  GIT COMMANDS
%  =========================================================================
fprintf('\n--- Git Commit Commands ---\n');
fprintf('cd %s\n', git_dir);
fprintf('git add .\n');
fprintf('git commit -m "Step 01: Setup, env check, dataset download, README"\n');
fprintf('git remote add origin https://github.com/<username>/fmri-spm-analysis.git\n');
fprintf('git push -u origin main\n');
fprintf('\n=== END OF STEP 01 ===\n');
%% =========================================================================