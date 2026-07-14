function fig_val1_island()
% FIG (validation 1)  Single-island gate, cf. Paul Fig. 1.
%   V_PD turns on near eps_crit = sqrt(kperp)/2 and scales as sqrt(eps).
%   Data: validation1_paul.mat  eps(16x3), V_PD(16x3), kperp[1e-2;1e-4;1e-6],
%   eps_crit(3x1).  NOTE: the eps grid differs per kperp (column j <-> kperp(j)).

S = load('validation1_paul.mat');
col = lines(3);
f = figure('Color','w','Position',[100 100 560 430]);
hold on; box on;
for j=1:3
    x = S.eps(:,j) / S.eps_crit(j);
    plot(x, S.V_PD(:,j), '-o', 'Color',col(j,:), 'MarkerFaceColor',col(j,:), ...
         'LineWidth',1.3, 'DisplayName',sprintf('\\kappa_\\perp = %g',S.kperp(j)));
end
xline(1,'k--','LineWidth',0.8,'HandleVisibility','off');
set(gca,'XScale','log','FontSize',10);
xlabel('\epsilon / \epsilon_{crit}   (\epsilon_{crit} = \surd\kappa_\perp/2)');
ylabel('V_{PD}');
title('Validation 1: island transition at \epsilon_{crit}');
legend('Location','northwest'); grid on;
print(f, 'fig_val1_island.png', '-dpng', '-r300');
fprintf('wrote fig_val1_island.png\n');
end
