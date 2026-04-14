%% 18探头双弧阵声场仿真 (z=0)
% 目标：恢复之前的双弧设计，使用18个探头分布在两条圆弧上，计算z=0平面的声压分布。
% 该脚本与直线阵、单弧阵、3x6阵保持相同的功率标定、网格和输出格式，便于横向比较。

clear; clc; close all;
rng(20260414);

%% 1) 物理参数
c = 343;                % 声速 (m/s)
f = 40000;              % 工作频率 (Hz)
lambda = c / f;         % 波长 (m)
k = 2 * pi / lambda;    % 波数 (rad/m)

%% 2) 阵列参数：双弧18探头
N_inner = 6;
N_outer = 12;
N = N_inner + N_outer;

% 双弧几何参数（可根据需要微调）
r1 = 0.055;             % 内弧半径 (m)
r2 = 0.117;             % 外弧半径 (m)
a1 = 74.1;              % 内弧张角半角 (deg)
a2 = 72.7;              % 外弧张角半角 (deg)
y1 = -0.008;            % 内弧整体平移 (m)
y2 = -0.003;            % 外弧整体平移 (m)

ang1 = linspace(-a1, a1, N_inner).';
ang2 = linspace(-a2, a2, N_outer).';

% 两条弧均以 +y 方向为主轴展开
x1 = r1 * sind(ang1);
y1p = r1 * cosd(ang1) + y1;

x2 = r2 * sind(ang2);
y2p = r2 * cosd(ang2) + y2;

pos = [x1, y1p; x2, y2p];

%% 3) 辐射模型与功率标定
% 单个探头指向性：cos^n(theta), -6 dB全角70度
theta_half = 35;
dir_n = -6 / (20 * log10(cosd(theta_half)));

% 总功率设定（可在1~4 W调整）
P_ref_W = 4.0;
P_total_W = 4.0;
P_total_W = min(max(P_total_W, 1.0), 4.0);
SPL_ref_at_Pref_dB = 120;  % 参考功率下，1m正前方标定声压级
SPL_ref_dB = SPL_ref_at_Pref_dB + 10 * log10(P_total_W / P_ref_W);

fprintf('18探头双弧阵总探头数: %d\n', N);
fprintf('内弧: r1=%.3f m, a1=%.1f deg, y1=%.3f m\n', r1, a1, y1);
fprintf('外弧: r2=%.3f m, a2=%.1f deg, y2=%.3f m\n', r2, a2, y2);
fprintf('总功率设定: %.1f W (参考 %.1f W) -> 参考点标定 %.1f dB SPL\n', P_total_W, P_ref_W, SPL_ref_dB);

%% 4) 观察平面网格（z=0平面）
r_vec = 0.3:0.05:5.0;
theta_vec = -90:1:90;
[R, Theta] = meshgrid(r_vec, theta_vec);
X = R .* sind(Theta);
Y = R .* cosd(Theta);

%% 5) 计算声场（复声压相干叠加）
P = zeros(size(R));
for m = 1:N
    dx = X - pos(m, 1);
    dy = Y - pos(m, 2);
    r_m = sqrt(dx.^2 + dy.^2);

    % 与 +y 方向的夹角，后向瓣抑制
    theta_m = atan2d(dx, dy);
    D = max(0, cosd(theta_m)).^dir_n;

    P = P + (D ./ r_m) .* exp(-1j * k * r_m);
end

P_abs = abs(P);

% 参考点：x=0, y=1m
[~, idx_ref] = min(abs(X(:)) + abs(Y(:) - 1));
P_ref = P_abs(idx_ref);
SPL = SPL_ref_dB + 20 * log10(P_abs / P_ref);
SPL_vis = min(max(SPL, 60), 140);

%% 6) 绘图：平面声压分布
figure('Name', '18探头双弧阵声压分布', 'Color', 'w');
pcolor(X, Y, SPL_vis);
shading interp;
colormap(jet);
colorbar;
caxis([60 140]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('18探头双弧阵在 z=0 平面的声压级分布');
xlim([-5, 5]);
ylim([0, 5]);
hold on;
plot(pos(:, 1), pos(:, 2), 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', '探头位置');
legend('SPL', '探头位置');

%% 7) 指向性曲线与轴向曲线
figure('Name', '不同距离处的角度声压变化', 'Color', 'w');
distances = [1, 2, 3, 4, 5];
cc = lines(numel(distances));
hold on;
for i = 1:numel(distances)
    [~, col_r] = min(abs(r_vec - distances(i)));
    plot(theta_vec, SPL(:, col_r), 'LineWidth', 1.5, 'Color', cc(i, :), ...
        'DisplayName', sprintf('r = %.1f m', distances(i)));
end
xlabel('角度 (deg)');
ylabel('SPL (dB)');
title('18探头双弧阵在不同距离处的角度声压变化');
grid on;
legend('Location', 'best');
xlim([-90, 90]);

[~, idx_theta0] = min(abs(theta_vec - 0));
figure('Name', '轴向声压曲线', 'Color', 'w');
plot(r_vec, SPL(idx_theta0, :), 'r-', 'LineWidth', 1.8);
grid on;
xlabel('距离 r (m)');
ylabel('SPL (dB)');
title('18探头双弧阵轴向（theta = 0°）声压随距离变化');
xlim([0.3, 5]);

%% 8) 输出结果
[max_SPL, idx_max] = max(SPL(:));
fprintf('最大声压级: %.1f dB SPL\n', max_SPL);
fprintf('最大值位置: x = %.2f m, y = %.2f m\n', X(idx_max), Y(idx_max));

disp('18探头双弧阵仿真完成。');
