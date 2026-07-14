function fig1_threeaxis()
% FIG1  The core figure: 4 fields (rows) x 3 indicators (columns).
%   Rows: m = 4, 12, 20, 36.  Columns: converse-KAM t_c, WBA dig, V_PD chi.
%   Same (rho,theta) grid; shared colour scale down each column.
%   Data: stage3_m{m}.mat (24x24).  chi is at kperp = 1e-6.

fields = [4 12 20 36];

% ---- shared colour limits per column (computed from the data) ----
tc_all = [];
for m = fields
    S = load(sprintf('stage3_m%d.mat', m));
    v = S.tc(:); tc_all = [tc_all; v(isfinite(v))]; %#ok<AGROW>
end
tc_clim  = [min(tc_all), 2*median(tc_all)];   % robust, toolbox-free (~[10 40])
dig_clim = [0 14];                            % gate: regular ~14, chaotic ~1
chi_clim = [0 1];                             % binary indicator

binmap = [0.20 0.20 0.28; 0.98 0.90 0.14];    % chi: perp-dominated / parallel

f = figure('Color','w','Position',[100 100 780 900]);
for r = 1:4
    m = fields(r);
    S = load(sprintf('stage3_m%d.mat', m));
    tc = S.tc; tc(~logical(S.detected)) = NaN;   % undetected -> NaN
    thn = S.theta / max(S.theta);                % common x-axis, [0,1]

    % --- column 1: t_c ---
    ax = subplot(4,3,(r-1)*3+1);
    drawmap(ax, thn, S.rho, tc, tc_clim, viridis(256), true);
    ylabel(sprintf('m = %d\n\\rho', m));
    if r==1, title('converse-KAM  t_c'); end

    % --- column 2: dig ---
    ax = subplot(4,3,(r-1)*3+2);
    drawmap(ax, thn, S.rho, S.dig, dig_clim, viridis(256), false);
    if r==1, title('WBA  dig'); end

    % --- column 3: chi ---
    ax = subplot(4,3,(r-1)*3+3);
    drawmap(ax, thn, S.rho, S.chi, chi_clim, binmap, false);
    if r==1, title('V_{PD} indicator  \chi'); end

    if r==4
        for c=1:3, subplot(4,3,9+c); xlabel('\theta / (2\pi/m)'); end
    end
end

% one colorbar per column (attached to the bottom-row axes)
addcbar(subplot(4,3,10), tc_clim);
addcbar(subplot(4,3,11), dig_clim);
addcbar(subplot(4,3,12), chi_clim);

print(f, 'fig1_threeaxis.png', '-dpng', '-r300');
fprintf('wrote fig1_threeaxis.png\n');
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
