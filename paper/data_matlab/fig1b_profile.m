function fig1b_profile()
% FIG1b  Standalone companion to Fig 1: angle-averaged temperature profile
%   T(rho) = <T(rho,theta)>_theta over the core rho in [0.25,0.75], four fields.
%   DeltaT = <T(0.75)> - <T(0.25)> is computed from the data and shown in the
%   legend. Data: vpd_m{m}.mat (T_zeta0 129x32, zeta=0, kperp=1e-6).
%   Use this if Fig 1 (maps) and the profile go as separate figures (option B).

fields = [4 12 20 36];
ls = {'-','--','-.',':'};
col = lines(4);

f = figure('Color','w','Position',[100 100 560 430]); hold on; box on;
for i = 1:4
    m = fields(i);
    V = load(sprintf('vpd_m%d.mat', m));
    Tbar = mean(V.T_zeta0, 2);
    rho = V.rho(:); core = (rho >= 0.25) & (rho <= 0.75);
    dT = interp1(rho, Tbar, 0.75) - interp1(rho, Tbar, 0.25);
    plot(rho(core), Tbar(core), ls{i}, 'Color',col(i,:), 'LineWidth',1.8, ...
         'DisplayName', sprintf('m=%d  (\\DeltaT=%.3f)', m, dT));
end
set(gca,'FontSize',11); xlim([0.25 0.75]);
xlabel('\rho'); ylabel('\langle T\rangle_\theta');
title('Angle-averaged core temperature: m=4 flat, m=36 graded');
legend('Location','northwest'); grid on;

exportgraphics(f, 'fig1b_profile.pdf', 'ContentType','vector');
print(f, 'fig1b_profile.png', '-dpng', '-r300');
fprintf('wrote fig1b_profile.pdf (vector) + .png\n');
end
