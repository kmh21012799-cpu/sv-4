function fig2_transport()
% FIG2  Transport separates the four fields.
%   (a) V_PD vs kperp,  (b) DeltaT vs kperp,  four curves each.
%   Data: stage2_paul.mat  (V_PD, DeltaT are 4x3: rows m, cols kperp).

S = load('stage2_paul.mat');
m = S.m(:); kperp = S.kperp(:);
col = lines(4);
mk = {'-o','-s','-^','-d'};

f = figure('Color','w','Position',[100 100 900 380]);

subplot(1,2,1); hold on; box on;
for i=1:4
    plot(kperp, S.V_PD(i,:), mk{i}, 'Color',col(i,:), ...
         'MarkerFaceColor',col(i,:), 'LineWidth',1.4, 'DisplayName',sprintf('m=%d',m(i)));
end
set(gca,'XScale','log','XDir','reverse','FontSize',10);
xlabel('\kappa_\perp'); ylabel('V_{PD}');
title('(a)  V_{PD}: largest at m=4'); legend('Location','northwest'); grid on;

subplot(1,2,2); hold on; box on;
for i=1:4
    plot(kperp, S.DeltaT(i,:), mk{i}, 'Color',col(i,:), ...
         'MarkerFaceColor',col(i,:), 'LineWidth',1.4, 'DisplayName',sprintf('m=%d',m(i)));
end
set(gca,'XScale','log','XDir','reverse','FontSize',10);
xlabel('\kappa_\perp'); ylabel('\Delta T');
title('(b)  \Delta T: largest at m=36 (best insulator)'); legend('Location','northeast'); grid on;

print(f, 'fig2_transport.png', '-dpng', '-r300');
fprintf('wrote fig2_transport.png\n');
end
