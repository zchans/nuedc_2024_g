%% 18探头直线阵声场仿真 (z=0)
% 目标：作为统一对比基准，使用18个探头沿直线排列，计算z=0平面的声压分布。
% 该脚本与单弧阵列、3x6阵列保持相同的功率标定、网格和输出格式，便于横向比较。

clear; clc; close all;
rng(20260414);

%% 1) 物理参数
c = 343;                % 声速 (m/s)
f = 40000;              % 工作频率 (Hz)
lambda = c / f;         % 波长 (m)
k = 2 * pi / lambda;    % 波数 (rad/m)

%% 2) 阵列参数：18探头直线阵
N = 18;                 % 探头数量
pitch = 0.016;          % 相邻探头中心距 (m)
aperture = (N - 1) * pitch;

% 18个探头沿x轴等间距排布，中心对称，全部位于y=0平面
x_array = linspace(-aperture / 2, aperture / 2, N).';
y_array = zeros(N, 1);
pos = [x_array, y_array];

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

fprintf('18探头直线阵总探头数: %d\n', N);
fprintf('总孔径: %.3f m, 单个间距: %.3f m\n', aperture, pitch);
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
figure('Name', '18探头直线阵声压分布', 'Color', 'w');
pcolor(X, Y, SPL_vis);
shading interp;
colormap(jet);
colorbar;
caxis([60 140]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('18探头直线阵在 z=0 平面的声压级分布');
xlim([-5, 5]);
ylim([0, 5]);
hold on;
plot(x_array, y_array, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 5, 'DisplayName', '探头位置');
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
title('18探头直线阵在不同距离处的角度声压变化');
grid on;
legend('Location', 'best');
xlim([-90, 90]);

[~, idx_theta0] = min(abs(theta_vec - 0));
figure('Name', '轴向声压曲线', 'Color', 'w');
plot(r_vec, SPL(idx_theta0, :), 'r-', 'LineWidth', 1.8);
grid on;
xlabel('距离 r (m)');
ylabel('SPL (dB)');
title('18探头直线阵轴向（theta = 0°）声压随距离变化');
xlim([0.3, 5]);

%% 8) 输出结果
[max_SPL, idx_max] = max(SPL(:));
fprintf('最大声压级: %.1f dB SPL\n', max_SPL);
fprintf('最大值位置: x = %.2f m, y = %.2f m\n', X(idx_max), Y(idx_max));

disp('18探头直线阵仿真完成。');
