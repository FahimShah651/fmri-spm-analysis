function status = check_environment()
%CHECK_ENVIRONMENT  Verify MATLAB, SPM, and Git are properly configured.
%
%  Author     : Fahim Ur Rehman Shah (EE2629)
%  Course     : CSE532 — Signal and Image Processing (MS Level)
%  Supervisor : Dr. Adnan Shah, FCSE, GIKI
%  Date       : May 2026
%
%  Returns status struct with fields: matlab_ok, spm_ok, git_ok, all_ok

status.matlab_ok = false;
status.spm_ok    = false;
status.git_ok    = false;
status.all_ok    = false;

fprintf('=== Environment Check ===\n');

% MATLAB version
v = ver('MATLAB');
matlab_ver = str2double(v.Version);
status.matlab_ok = matlab_ver >= 9.3;  % R2022b = 9.13
fprintf('  MATLAB : %s %s — %s\n', v.Version, v.Release, ...
    iif(status.matlab_ok,'OK','WARN (R2022b+ recommended)'));

% SPM
try
    spm_path = fileparts(which('spm'));
    if ~isempty(spm_path)
        spm_ver = spm('Ver');
        status.spm_ok = true;
        fprintf('  SPM    : %s at %s — OK\n', spm_ver, spm_path);
    else
        fprintf('  SPM    : Not on MATLAB path — FAIL\n');
    end
catch
    fprintf('  SPM    : Error checking SPM — FAIL\n');
end

% Git
[rc, out] = system('git --version');
status.git_ok = (rc == 0);
if status.git_ok
    fprintf('  Git    : %s — OK\n', strtrim(out));
else
    fprintf('  Git    : Not found — WARN\n');
end

status.all_ok = status.matlab_ok && status.spm_ok;
fprintf('=== %s ===\n', iif(status.all_ok,'ALL GOOD','ISSUES FOUND'));

function s = iif(c,t,f)
    if c; s = t; else; s = f; end
end
end
