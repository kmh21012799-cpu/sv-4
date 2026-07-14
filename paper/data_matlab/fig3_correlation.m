function fig3_correlation()
% FIG3  Paul's expectation refuted: r(converse-KAM, V_PD) does not grow as
%   kperp -> 0.  Graded t_c correlation r_tc_vpd (4x3: rows m, cols kperp).
%   Data: correlations_paul.mat.

S = load('correlations_paul.mat');
m = S.m(:); kperp = S.kperp(:);
col = lines(4);
mk = {'-o','-s','-^','-d'};

f = figure('Color','w','Position',[100 100 560 430]);
hold on; box on;
for i=1:4
    plot(kperp, S.r_tc_vpd(i,:), mk{i}, 'Color',col(i,:), ...
         'MarkerFaceColor',col(i,:), 'LineWidth',1.4, 'DisplayName',sprintf('m=%d',m(i)));
end
yline(0,'k-','LineWidth',0.6,'HandleVisibility','off');
set(gca,'XScale','log','XDir','reverse','FontSize',10);
xlabel('\kappa_\perp');
ylabel('r(  -t_c ,  \chi )   [Spearman]');
title({'converse-KAM  \leftrightarrow  V_{PD} correlation','does not grow as \kappa_\perp \rightarrow 0'});
ylim([-0.15 0.25]); legend('Location','northwest'); grid on;

print(f, 'fig3_correlation.png', '-dpng', '-r300');
fprintf('wrote fig3_correlation.png\n');
end
