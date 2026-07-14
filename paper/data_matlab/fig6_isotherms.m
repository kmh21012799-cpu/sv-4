function fig6_isotherms()
% FIG6  Isotherms (T contours) and the parallel-dominated region (chi=1),
%   zeta = 0 slice, kperp = 1e-6 (cf. Paul Fig. 6).
%   Data: vpd_m{m}.mat  (T_zeta0, chi_zeta0 129x32; rho 129x1; theta 32x1).

fields = [4 12 20 36];
f = figure('Color','w','Position',[100 100 900 760]);
for i=1:4
    m = fields(i);
    S = load(sprintf('vpd_m%d.mat', m));
    thn = S.theta / max(S.theta);              % normalise theta to [0,1]
    ax = subplot(2,2,i); hold on; box on;
    % chi = 1 region (parallel-dominated) as a light background
    imagesc(thn, S.rho, S.chi_zeta0); axis xy;
    colormap(ax, [1 1 1; 0.80 0.88 1.0]);      % 0 -> white, 1 -> light blue
    caxis(ax,[0 1]);
    % isotherms of T on top
    [C,h] = contour(thn, S.rho, S.T_zeta0, 0:0.1:1, 'k', 'LineWidth',0.7);
    clabel(C,h,0:0.2:1,'FontSize',6);
    xlim([0 1]); ylim([0 1]);
    set(ax,'FontSize',9,'XTick',[0 0.5 1],'YTick',0:0.25:1);
    xlabel('\theta / (2\pi/m)'); ylabel('\rho');
    title(sprintf('m = %d', m));
end
sgtitle('Isotherms (black) + parallel-dominated region \chi=1 (blue),  \kappa_\perp=10^{-6}');
print(f, 'fig6_isotherms.png', '-dpng', '-r300');
fprintf('wrote fig6_isotherms.png\n');
end
