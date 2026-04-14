%% 单弧18探头声场仿真 (z=0)
% 目标：将18个探头均匀排布在同一条圆弧上，使用相干叠加计算z=0平面的声压分布。
% 用途：对比单弧阵列与直线阵/双弧阵列的声场覆盖效果。

clear; clc; close all;
rng(20260414);

%% 1) 物理参数
c = 343;                % 声速 (m/s)
f = 40000;              % 工作频率 (Hz)
lambda = c / f;         % 波长 (m)
k = 2 * pi / lambda;    % 波数 (rad/m)

%% 2) 阵列参数：单弧18探头
N = 18;                 % 探头数量
arc_radius = 0.095;     % 圆弧半径 (m)
arc_span_deg = 150;     % 圆弧张角 (deg)
arc_shift_y = 0.000;    % 圆弧整体平移 (m)

% 探头在圆弧上的角度，中心对称分布
theta_arc = linspace(-arc_span_deg/2, arc_span_deg/2, N).';

% 单弧坐标：以 +y 方向为主轴展开
x_array = arc_radius * sind(theta_arc);
y_array = arc_radius * cosd(theta_arc) + arc_shift_y;
pos = [x_array, y_array];

% 探头最小间距检查
min_d = inf;
for i = 1:N
    for j = i+1:N
        d_ij = hypot(pos(i,1) - pos(j,1), pos(i,2) - pos(j,2));
        if d_ij < min_d
            min_d = d_ij;
        end
    end
end

%% 3) 辐射模型与功率标定
% 单个探头指向性：cos^n(theta), -6 dB全角70度
theta_half = 35;
dir_n = -6 / (20 * log10(cosd(theta_half)));

% 总功率设定（可在 1~4 W 调整）
P_ref_W = 4.0;
P_total_W = 4.0;
P_total_W = min(max(P_total_W, 1.0), 4.0);
SPL_ref_at_Pref_dB = 120;  % 参考功率下，1m正前方标定声压级
SPL_ref_dB = SPL_ref_at_Pref_dB + 10 * log10(P_total_W / P_ref_W);

fprintf('单弧阵列总探头数: %d\n', N);
fprintf('圆弧半径: %.3f m, 张角: %.1f deg, 最小中心距: %.3f m\n', arc_radius, arc_span_deg, min_d);
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
figure('Name', '单弧18探头声压分布', 'Color', 'w');
pcolor(X, Y, SPL_vis);
shading interp;
colormap(jet);
colorbar;
caxis([60 140]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('单弧18探头在 z=0 平面的声压级分布');
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
title('单弧阵列在不同距离处的角度声压变化');
grid on;
legend('Location', 'best');
xlim([-90, 90]);

[~, idx_theta0] = min(abs(theta_vec - 0));
figure('Name', '轴向声压曲线', 'Color', 'w');
plot(r_vec, SPL(idx_theta0, :), 'r-', 'LineWidth', 1.8);
grid on;
xlabel('距离 r (m)');
ylabel('SPL (dB)');
title('单弧阵列轴向（theta = 0°）声压随距离变化');
xlim([0.3, 5]);

%% 8) 输出结果
[max_SPL, idx_max] = max(SPL(:));
fprintf('最大声压级: %.1f dB SPL\n', max_SPL);
fprintf('最大值位置: x = %.2f m, y = %.2f m\n', X(idx_max), Y(idx_max));

disp('单弧18探头仿真完成。');
