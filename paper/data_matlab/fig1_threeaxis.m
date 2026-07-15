function fig1_threeaxis()
% FIG1  The core figure (layout A): a 4x3 map panel + a companion line profile.
%   Rows m = 4,12,20,36.
%   Col 1: converse-KAM t_c   (stage3_m{m}.mat, 24x24, core rho in [0.25,0.75])
%   Col 2: WBA dig            (stage3_m{m}.mat, 24x24, core)
%   Col 3: Temperature T      (vpd_m{m}.mat, 129x32, RESTRICTED to the core
%                              rho in [0.25,0.75], zeta=0 slice, kperp=1e-6)
%   All three columns share the same rho axis [0.25,0.75]; theta is normalised
%   to [0,1] (its range 2*pi/m differs per field). t_c and dig share a colour
%   scale per column; T shares 0..1 across the four fields (no rescaling).
%   Bottom row (spanning): angle-averaged profile T(rho)=<T>_theta for the four
%   fields, with DeltaT computed from the data shown in the legend.
%   Requires R2019b+ (tiledlayout/yline). Saved as vector PDF + PNG.

fields = [4 12 20 36];
ls = {'-','--','-.',':'};
col = lines(4);

% ---- shared colour limits (t_c/dig from data; T fixed 0..1) ----
tc_all = [];
for m = fields
    S = load(sprintf('stage3_m%d.mat', m));
    v = S.tc(:); tc_all = [tc_all; v(isfinite(v))]; %#ok<AGROW>
end
tc_clim = [min(tc_all), 2*median(tc_all)];
dig_clim = [0 14];
T_clim = [0 1];

f = figure('Color','w','Position',[100 100 780 1050]);
tl = tiledlayout(5,3,'Padding','compact','TileSpacing','compact');

axcol = gobjects(1,3);   % keep one axes per column for the colorbars
for r = 1:4
    m = fields(r);
    S = load(sprintf('stage3_m%d.mat', m));
    V = load(sprintf('vpd_m%d.mat', m));
    tc = S.tc; tc(~logical(S.detected)) = NaN;
    thn  = S.theta / max(S.theta);
    thnT = V.theta / max(V.theta);
    coreT = (V.rho >= 0.25) & (V.rho <= 0.75);   % restrict T to the core

    % col 1: t_c
    ax = nexttile; drawmap(ax, thn, S.rho, tc, tc_clim, viridis(256), true);
    ylabel(sprintf('m = %d\n\\rho', m));
    if r==1, title('converse-KAM  t_c'); end
    if r==4, xlabel('\theta / (2\pi/m)'); axcol(1)=ax; end

    % col 2: dig
    ax = nexttile; drawmap(ax, thn, S.rho, S.dig, dig_clim, viridis(256), false);
    if r==1, title('WBA  dig'); end
    if r==4, xlabel('\theta / (2\pi/m)'); axcol(2)=ax; end

    % col 3: Temperature T (core only)
    ax = nexttile; drawmap(ax, thnT, V.rho(coreT), V.T_zeta0(coreT,:), T_clim, viridis(256), false);
    if r==1, title('Temperature  T'); end
    if r==4, xlabel('\theta / (2\pi/m)'); axcol(3)=ax; end
end

% per-column colorbars (each reflects its own axes caxis)
cb1 = colorbar(axcol(1),'southoutside'); cb1.Label.String='t_c'; cb1.FontSize=7;
cb2 = colorbar(axcol(2),'southoutside'); cb2.Label.String='dig'; cb2.FontSize=7;
cb3 = colorbar(axcol(3),'southoutside'); cb3.Label.String='T';   cb3.FontSize=7;

% ---- bottom row: companion angle-averaged profile ----
axp = nexttile([1 3]); hold(axp,'on'); box(axp,'on');
for i = 1:4
    m = fields(i);
    V = load(sprintf('vpd_m%d.mat', m));
    Tbar = mean(V.T_zeta0, 2);                     % <T>_theta, 129x1
    rho = V.rho(:); core = (rho >= 0.25) & (rho <= 0.75);
    dT = interp1(rho, Tbar, 0.75) - interp1(rho, Tbar, 0.25);   % from the data
    plot(axp, rho(core), Tbar(core), ls{i}, 'Color',col(i,:), 'LineWidth',1.8, ...
         'DisplayName', sprintf('m=%d  (\\DeltaT=%.3f)', m, dT));
end
set(axp,'FontSize',10); xlim(axp,[0.25 0.75]);
xlabel(axp,'\rho'); ylabel(axp,'\langle T\rangle_\theta');
title(axp,'Angle-averaged temperature across the core');
legend(axp,'Location','northwest'); grid(axp,'on');

exportgraphics(f, 'fig1_threeaxis.pdf', 'ContentType','vector');
print(f, 'fig1_threeaxis.png', '-dpng', '-r300');
fprintf('wrote fig1_threeaxis.pdf (vector) + .png\n');
end

function drawmap(ax, x, y, C, clim, cmap, showNaN)
axes(ax); %#ok<LAXES>
h = imagesc(x, y, C);
if showNaN
    set(h, 'AlphaData', ~isnan(C));
    set(ax, 'Color', [0.8 0.8 0.8]);
end
axis(ax,'xy'); axis(ax,'tight'); ylim(ax,[0.25 0.75]);
caxis(ax, clim); colormap(ax, cmap);
set(ax,'FontSize',8,'XTick',[0 1],'YTick',[0.25 0.5 0.75]);
end
