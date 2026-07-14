function fig4_wba_tconv()
% FIG4  WBA dig T-convergence: the short-time ordering is a finite-time
%   artifact; the curves cross (m=36 falls, m=4 rises) as T grows.
%   Data: wba_tconv.mat  (dig_median 4x4: rows m, cols T).

S = load('wba_tconv.mat');
m = S.m(:); T = S.T(:);
col = lines(4);
mk = {'-o','-s','-^','-d'};

f = figure('Color','w','Position',[100 100 560 430]);
hold on; box on;
for i=1:4
    plot(T, S.dig_median(i,:), mk{i}, 'Color',col(i,:), ...
         'MarkerFaceColor',col(i,:), 'LineWidth',1.5, 'DisplayName',sprintf('m=%d',m(i)));
end
set(gca,'XScale','log','FontSize',10);
xlabel('integration time  T  (toroidal periods)');
ylabel('median core WBA dig');
title('WBA dig vs T: no stable ordering (curves cross)');
xlim([400 6000]); legend('Location','northeast'); grid on;

print(f, 'fig4_wba_tconv.png', '-dpng', '-r300');
fprintf('wrote fig4_wba_tconv.png\n');
end
