function fig5_poincare()
% FIG5  Poincare sections (section zeta = 0) of the four Paul fields.
%   Chaotic sea (black) fills [1/8, 7/8] in all four; regular seeds (blue)
%   trace surviving KAM curves near the boundaries + island O-points.
%   Data: poincare_m{m}.mat  (rho, theta, kind, res_psi).

fields = [4 12 20 36];
f = figure('Color','w','Position',[100 100 820 720]);
for i=1:4
    m = fields(i);
    S = load(sprintf('poincare_m%d.mat', m));
    reg = (S.kind==0); ch = (S.kind==1);
    ax = subplot(2,2,i); hold on; box on;
    % resonance reference lines
    for k=1:numel(S.res_psi)
        yline(S.res_psi(k), 'Color',[0.85 0.85 0.85], 'LineWidth',0.3);
    end
    plot(S.theta(ch),  S.rho(ch),  '.', 'MarkerSize',1, 'Color',[0 0 0]);
    plot(S.theta(reg), S.rho(reg), '.', 'MarkerSize',1, 'Color',[0.10 0.35 0.80]);
    yline(1/8,'r--','LineWidth',0.6); yline(7/8,'r--','LineWidth',0.6);
    xlim([0 2*pi]); ylim([0 1]);
    set(ax,'FontSize',9,'XTick',[0 pi 2*pi],'XTickLabel',{'0','\pi','2\pi'});
    xlabel('\theta'); ylabel('\psi = \rho');
    title(sprintf('m = %d   (%d chains)', m, numel(S.res_psi)));
end
sgtitle('Poincare sections (\zeta = 0), Paul critical-overlap fields');
print(f, 'fig5_poincare.png', '-dpng', '-r300');
fprintf('wrote fig5_poincare.png\n');
end
