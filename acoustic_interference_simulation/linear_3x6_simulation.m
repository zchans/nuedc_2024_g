%% 3x6超声波探头板声场仿真
% 6个探头位置沿直线排列，每个位置叠加3个同相探头，共18个探头。
% 目标：观察“3个探头叠加为1个位置”的3x6板在z=0平面的声压分布。

clear; clc; close all;
rng(20260414);

%% 1) 物理参数
c = 343;                % 声速 (m/s)
f = 40000;              % 工作频率 (Hz)
lambda = c / f;         % 波长 (m)
k = 2 * pi / lambda;    % 波数 (rad/m)

%% 2) 阵列参数
num_positions = 6;      % 线阵位置数
num_per_position = 3;   % 每个位置叠加探头数
N = num_positions * num_per_position;

pitch = 0.016;          % 相邻位置间距 (m)
line_length = (num_positions - 1) * pitch;

% 6个位置沿x轴等间距排布，中心对称，全部位于y=0平面
x_positions = linspace(-line_length / 2, line_length / 2, num_positions).';
y_positions = zeros(num_positions, 1);

% 将每个位置复制3次，表示同一点叠加3个探头
x_array = repelem(x_positions, num_per_position);
y_array = repelem(y_positions, num_per_position);
pos = [x_array, y_array];

%% 3) 探头辐射与功率标定
% 单个探头指向性：cos^n(theta)，-6 dB全角70度
theta_half = 35;
dir_n = -6 / (20 * log10(cosd(theta_half)));

% 总功率设定（可改为1~4 W）
P_ref_W = 4.0;
P_total_W = 4.0;
P_total_W = min(max(P_total_W, 1.0), 4.0);
SPL_ref_at_Pref_dB = 120;  % 参考功率下，1m正前方标定声压级
SPL_ref_dB = SPL_ref_at_Pref_dB + 10 * log10(P_total_W / P_ref_W);

fprintf('3x6阵列总探头数: %d\n', N);
fprintf('总功率设定: %.1f W (参考 %.1f W) -> 参考点标定 %.1f dB SPL\n', P_total_W, P_ref_W, SPL_ref_dB);

%% 4) 观察平面网格（z=0平面）
r_vec = 0.3:0.05:5.0;
theta_vec = -90:1:90;
[R, Theta] = meshgrid(r_vec, theta_vec);
X = R .* sind(Theta);
Y = R .* cosd(Theta);

%% 5) 计算声场
P = zeros(size(R));
for m = 1:N
    dx = X - pos(m, 1);
    dy = Y - pos(m, 2);
    r_m = sqrt(dx.^2 + dy.^2);

    % 与+Y轴夹角，超出前向半空间后贡献迅速衰减
    theta_m = atan2d(dx, dy);
    D = max(0, cosd(theta_m)).^dir_n;

    % 同一点叠加的3个探头视为同相同幅辐射，直接叠加复压
    P = P + (D ./ r_m) .* exp(-1j * k * r_m);
end

P_abs = abs(P);

% 参考点：x=0, y=1m
[~, idx_ref] = min(abs(X(:)) + abs(Y(:) - 1));
P_ref = P_abs(idx_ref);
SPL = SPL_ref_dB + 20 * log10(P_abs / P_ref);
SPL_vis = min(max(SPL, 60), 140);

%% 6) 绘图：平面声压分布
figure('Name', '3x6探头板声压分布', 'Color', 'w');
pcolor(X, Y, SPL_vis);
shading interp;
colormap(jet);
colorbar;
caxis([60 140]);
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('3x6超声波探头板在 z=0 平面的声压级分布');
xlim([-5, 5]);
ylim([0, 5]);
hold on;
plot(x_positions, y_positions, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 6, 'DisplayName', '6个位置');
legend('SPL', '探头位置');

%% 7) 轴向声压曲线
[~, idx_theta0] = min(abs(theta_vec - 0));
figure('Name', '轴向声压曲线', 'Color', 'w');
plot(r_vec, SPL(idx_theta0, :), 'r-', 'LineWidth', 1.8);
grid on;
xlabel('距离 r (m)');
ylabel('SPL (dB)');
title('3x6阵列轴向（theta=0°）声压随距离变化');
xlim([0.3, 5]);

%% 8) 输出结果
[max_SPL, idx_max] = max(SPL(:));
fprintf('最大声压级: %.1f dB SPL\n', max_SPL);
fprintf('最大值位置: x = %.2f m, y = %.2f m\n', X(idx_max), Y(idx_max));

disp('3x6超声波探头板仿真完成。');
