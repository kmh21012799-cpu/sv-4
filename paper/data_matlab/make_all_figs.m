function make_all_figs()
% MAKE_ALL_FIGS  Run every figure script in this folder.
%   Run from inside paper/matlab_figs/ (the .mat data and viridis.m are here).
%   Each script writes a PNG at 300 dpi. Swap print(...) for
%   exportgraphics(gcf,'name.pdf','ContentType','vector') if you want vector PDF.

here = fileparts(mfilename('fullpath'));
cd(here);

fig1_threeaxis();    % 4x3 core panel (t_c / dig / chi)
fig2_transport();    % V_PD & DeltaT vs kperp
fig3_correlation();  % r(cKAM,V_PD) vs kperp (Paul's expectation refuted)
fig4_wba_tconv();    % WBA dig T-convergence (curves cross)
fig5_poincare();     % Poincare sections, four fields
fig6_isotherms();    % isotherms + parallel-dominated region
fig_val1_island();   % validation 1 (Paul Fig. 1)
fig_val2_chaos();    % validation 2 (Paul Fig. 3)

fprintf('\nAll figures written.\n');
end
