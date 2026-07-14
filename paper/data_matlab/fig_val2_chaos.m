function fig_val2_chaos()
% FIG (validation 2)  Chaotic-layer gate, cf. Paul Fig. 3.
%   V_PD rises 0 -> 1 as kperp decreases; the transition tracks the
%   quasilinear D_QL (dotted verticals).
%   Data: validation2_paul.mat  kperp(9x1), V_PD(9x3), S[2;2.83;4],
%   D_QL(3x1), kperp_trans(3x1).

S = load('validation2_paul.mat');
col = lines(3);
f = figure('Color','w','Position',[100 100 560 430]);
hold on; box on;
for j=1:3
    plot(S.kperp, S.V_PD(:,j), '-o', 'Color',col(j,:), 'MarkerFaceColor',col(j,:), ...
         'LineWidth',1.3, 'DisplayName',sprintf('S = %.2f',S.S(j)));
    xline(S.D_QL(j), ':', 'Color',col(j,:), 'LineWidth',1.2, 'HandleVisibility','off');
end
yline(0.5,'k--','LineWidth',0.6,'HandleVisibility','off');
set(gca,'XScale','log','XDir','reverse','FontSize',10);
xlabel('\kappa_\perp');
ylabel('V_{PD}');
title('Validation 2: chaotic-layer transition \propto D_{QL} (dotted)');
legend('Location','northwest'); grid on;
print(f, 'fig_val2_chaos.png', '-dpng', '-r300');
fprintf('wrote fig_val2_chaos.png\n');
end
