MATLAB figure scripts — converse-KAM / WBA / V_PD three-axis comparison
=======================================================================

Self-contained: this folder has the scripts, the colormap helper (viridis.m),
and all the .mat data. Run from inside this folder.

    >> cd paper/matlab_figs
    >> make_all_figs        % runs everything, writes PNGs at 300 dpi

or run any script on its own, e.g.  >> fig1_threeaxis

Scripts
-------
  fig1_threeaxis.m   Core figure. 4 fields (rows m=4/12/20/36) x 3 columns:
                     col1 converse-KAM t_c, col2 WBA dig (both 24x24 core
                     grid, rho in [0.25,0.75]), col3 Temperature T (129x32,
                     FULL domain rho in [0,1], zeta=0 slice, kperp=1e-6).
                     t_c/dig share a colour scale per column; T shares 0..1
                     across the four fields. The core band [0.25,0.75] is
                     dashed on the T panels. Saves vector PDF + PNG.
  fig2_transport.m   V_PD and DeltaT vs kperp (four curves each). Transport
                     separates the fields: V_PD largest at m=4, DeltaT at m=36.
  fig3_correlation.m r(converse-KAM, V_PD) vs kperp. Stays weak, does not grow
                     as kperp -> 0 (Paul's expectation not observed).
  fig4_wba_tconv.m   WBA dig vs integration time T; curves cross (finite-time
                     artifact), so no stable ordering.
  fig5_poincare.m    Poincare sections (zeta=0) of the four fields.
  fig6_isotherms.m   Isotherms (T contours) + parallel-dominated region (chi=1),
                     kperp=1e-6 slice (cf. Paul Fig. 6).
  fig_val1_island.m  Validation 1: single-island transition at eps_crit (Paul Fig.1)
  fig_val2_chaos.m   Validation 2: chaotic-layer transition ~ D_QL (Paul Fig.3)

Requirements
------------
  Base MATLAB only (no toolboxes). viridis.m is included and must stay in this
  folder (it is on the path when you run from here).
  yline / xline / sgtitle need R2018b+. If you are on an older release, replace
  yline(y) with plot(xlim,[y y]) and sgtitle with a title on the layout.
  For vector output, replace the print(...) line in each script with:
      exportgraphics(gcf, 'figX.pdf', 'ContentType','vector')

Data conventions (see also paper/data_matlab/README.md)
-------------------------------------------------------
  rho   in [0.25, 0.75] (core), common to all fields.
  theta in [0, 2*pi/m): RANGE DIFFERS PER FIELD. The panels normalise it as
        theta/max(theta) so the four rows share a horizontal axis.
  chi   is binary (0/1); V_PD is the fraction of chi=1. Spatial chi and T are
        at kperp = 1e-6 only.
  tc    contains Inf at undetected points; the scripts mask them with the
        'detected' array (shown grey).
  T     is a single zeta=0 slice, not a 3-D field.
  Poincare theta spans the full [0, 2*pi); kind = 0 regular seed, 1 chaotic.

Provenance
----------
  Corrected Paul-envelope field. Repo kmh21012799-cpu/sv-4, branch main.
  Data produced in C3b v2 (commit 8a52ca1); converse-KAM/WBA are the original
  tools/ code (393328e / bf5ca06). Poincare arrays recomputed and committed
  under paper/data/. wba_tconv values are from RECORD_C3a (bf5ca06).
