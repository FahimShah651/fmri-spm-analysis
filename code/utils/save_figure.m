function save_figure(fig, fig_dir, fig_name, dpi)
%SAVE_FIGURE  Save a figure as both PNG (raster) and EPS (vector) and
%             print the corresponding LaTeX \begin{figure}...\end{figure}
%             snippet to the console.
%
%  Usage:
%    save_figure(fig, fig_dir, fig_name)
%    save_figure(fig, fig_dir, fig_name, dpi)
%
%  Inputs:
%    fig      - figure handle
%    fig_dir  - output directory (must exist)
%    fig_name - filename without extension (e.g. 'fig01_ds1_raw_epi')
%    dpi      - PNG resolution in dpi (default: 150)
%
%  Author     : Fahim Ur Rehman Shah (EE2629)
%  Course     : CSE532 — Signal and Image Processing
%  Date       : May 2026

if nargin < 4 || isempty(dpi), dpi = 150; end

png_path = fullfile(fig_dir, [fig_name '.png']);
eps_path = fullfile(fig_dir, [fig_name '.eps']);

print(fig, png_path, '-dpng',  sprintf('-r%d', dpi));
print(fig, eps_path, '-depsc');
close(fig);

fprintf('  Saved: %s.png / .eps\n', fig_name);

% Print LaTeX snippet
caption = strrep(strrep(fig_name,'_',' '), 'fig', 'Figure ');
label   = strrep(fig_name, '_', ':');
fprintf('\n%%%% LaTeX snippet for %s.png\n', fig_name);
fprintf('\\begin{figure}[htbp]\n');
fprintf('\\centering\n');
fprintf('\\includegraphics[width=\\columnwidth]{figures/%s}\n', fig_name);
fprintf('\\caption{%s. ADD CAPTION HERE.}\n', caption);
fprintf('\\label{fig:%s}\n', label);
fprintf('\\end{figure}\n\n');
end
