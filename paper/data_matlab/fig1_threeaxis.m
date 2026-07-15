function fig1_threeaxis()
% FIG1  The core figure: 4 fields (rows) x 3 indicators (columns).
%   Rows: m = 4, 12, 20, 36.
%   Col 1: converse-KAM t_c   (stage3_m{m}.mat, 24x24, core rho in [0.25,0.75])
%   Col 2: WBA dig            (stage3_m{m}.mat, 24x24, core)
%   Col 3: Temperature T      (vpd_m{m}.mat, 129x32, FULL domain rho in [0,1],
%                              zeta=0 slice, kperp=1e-6)
%   theta is normalised to [0,1] in every panel (its range 2*pi/m differs per m).
%   Col 3 y-axis spans the full [0,1]; the core band [0.25,0.75] used by cols 1-2
%   is marked with dashed lines so the rows can be read together.
%   Saved as vector PDF (exportgraphics) + PNG preview.

fields = [4 12 20 36];

% ---- shared colour limits per column ----
tc_all = [];
for m = fields
    S = load(sprintf('stage3_m%d.mat', m));
    v = S.tc(:); tc_all = [tc_all; v(isfinite(v))]; %#ok<AGROW>
end
tc_clim  = [min(tc_all), 2*median(tc_all)];   % robust, toolbox-free (~[10 40])
dig_clim = [0 14];                            % gate: regular ~14, chaotic ~1
T_clim   = [0 1];                             % temperature, shared across fields

f = figure('Color','w','Position',[100 100 780 900]);
for r = 1:4
    m = fields(r);
    S  = load(sprintf('stage3_m%d.mat', m));   % t_c, dig (core grid)
    V  = load(sprintf('vpd_m%d.mat', m));      % T_zeta0 (full-domain grid)
    tc = S.tc; tc(~logical(S.detected)) = NaN; % undetected -> NaN
    thn  = S.theta / max(S.theta);             % core theta, normalised
    thnT = V.theta / max(V.theta);             % T-slice theta, normalised

    % --- column 1: t_c ---
    ax = subplot(4,3,(r-1)*3+1);
    drawmap(ax, thn, S.rho, tc, tc_clim, viridis(256), true);
    ylabel(sprintf('m = %d\n\\rho', m));
    if r==1, title('converse-KAM  t_c'); end

    % --- column 2: dig ---
    ax = subplot(4,3,(r-1)*3+2);
    drawmap(ax, thn, S.rho, S.dig, dig_clim, viridis(256), false);
    if r==1, title('WBA  dig'); end

    % --- column 3: Temperature T (full domain) ---
    ax = subplot(4,3,(r-1)*3+3);
    axes(ax); %#ok<LAXES>
    imagesc(thnT, V.rho, V.T_zeta0);
    axis(ax,'xy'); ylim(ax,[0 1]);
    caxis(ax, T_clim); colormap(ax, viridis(256));
    yline(0.25,'w--','LineWidth',0.6); yline(0.75,'w--','LineWidth',0.6); % core band
    set(ax,'FontSize',8,'XTick',[0 1],'YTick',[0 0.25 0.5 0.75 1]);
    if r==1, title('Temperature  T'); end

    if r==4
        for c=1:3, subplot(4,3,9+c); xlabel('\theta / (2\pi/m)'); end
    end
end

% one colorbar per column (attached to the bottom-row axes)
addcbar(subplot(4,3,10), tc_clim);
addcbar(subplot(4,3,11), dig_clim);
addcbar(subplot(4,3,12), T_clim);

exportgraphics(f, 'fig1_threeaxis.pdf', 'ContentType','vector');
print(f, 'fig1_threeaxis.png', '-dpng', '-r300');   % preview
fprintf('wrote fig1_threeaxis.pdf (vector) + .png\n');
end

function drawmap(ax, x, y, C, clim, cmap, showNaN)
axes(ax); %#ok<LAXES>
h = imagesc(x, y, C);
if showNaN
    set(h, 'AlphaData', ~isnan(C));
    set(ax, 'Color', [0.8 0.8 0.8]);   % undetected shown grey
end
axis(ax,'xy'); axis(ax,'tight');
caxis(ax, clim);
colormap(ax, cmap);
set(ax,'FontSize',8,'XTick',[0 1],'YTick',[0.25 0.5 0.75]);
end

function addcbar(ax, clim)
c = colorbar(ax,'southoutside');
caxis(ax, clim);
c.FontSize = 7;
end
