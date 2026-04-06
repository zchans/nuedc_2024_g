import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib

# 设置中文字体支持
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

# 参数设置
SAMPLE_RATE = 44100      # 采样率 (Hz)
DURATION = 0.05          # 信号时长 (秒)
NUM_SAMPLES = int(SAMPLE_RATE * DURATION)

# 频率设置
f1 = 700.0               # 频率1 (Hz)
f2 = 1100.0              # 频率2 (Hz)
amplitude = 0.8          # 信号幅度
nonlinear_gain = 0.5     # 非线性增益系数 (失真强度)

def generate_clean_signal(f1, f2, A, t):
    """生成纯净的双音信号"""
    return A * np.sin(2 * np.pi * f1 * t) + A * np.sin(2 * np.pi * f2 * t)

def generate_distorted_signal(f1, f2, A, alpha, t):
    """生成包含互调失真的信号 (三阶非线性模型)"""
    x = A * np.sin(2 * np.pi * f1 * t) + A * np.sin(2 * np.pi * f2 * t)
    # 非线性系统: y = x + α * x³
    y = x + alpha * x**3
    return y

def compute_spectrum(signal, sample_rate):
    """计算信号的功率谱"""
    n = len(signal)
    freq = np.fft.fftfreq(n, 1/sample_rate)
    spectrum = np.abs(np.fft.fft(signal)) / n
    # 只取正频率部分
    return freq[:n//2], spectrum[:n//2]

def plot_waveform(t, clean_signal, distorted_signal, f1, f2, alpha):
    """绘制时域波形对比图"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # 上半图：纯净信号
    ax1.plot(t, clean_signal, 'b-', linewidth=1.5, label='纯净信号')
    ax1.set_ylabel('幅度', fontsize=12)
    ax1.set_title(f'纯净双音信号 (f₁={f1}Hz, f₂={f2}Hz)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.legend(loc='upper right')
    ax1.set_xlim([t[0], t[-1]])
    ax1.set_ylim([-1.5, 1.5])
    
    # 下半图：失真信号
    ax2.plot(t, distorted_signal, 'r-', linewidth=1.5, label=f'失真信号 (α={alpha})')
    ax2.set_xlabel('时间 (秒)', fontsize=12)
    ax2.set_ylabel('幅度', fontsize=12)
    ax2.set_title('经过非线性系统后的信号（包含互调失真）', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--')
    ax2.legend(loc='upper right')
    ax2.set_xlim([t[0], t[-1]])
    ax2.set_ylim([-1.5, 1.5])
    
    plt.tight_layout()
    plt.savefig('waveform_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("✓ 波形图已保存: waveform_comparison.png")

def plot_spectrum(clean_signal, distorted_signal, sample_rate, f1, f2, alpha):
    """绘制频谱对比图"""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    
    # 计算频谱
    freq_clean, spec_clean = compute_spectrum(clean_signal, sample_rate)
    freq_dist, spec_dist = compute_spectrum(distorted_signal, sample_rate)
    
    # 限制频率范围到 2500Hz (包含所有重要的互调产物)
    max_freq = 2500
    idx_clean = freq_clean <= max_freq
    idx_dist = freq_dist <= max_freq
    
    # 上半图：纯净信号频谱
    ax1.stem(freq_clean[idx_clean], spec_clean[idx_clean], 
             basefmt=" ", linefmt='b-', markerfmt='bo', markersize=4)
    ax1.set_ylabel('幅度', fontsize=12)
    ax1.set_title('纯净信号频谱（只有两个基频分量）', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlim([0, max_freq])
    ax1.set_ylim([0, max(spec_clean) * 1.1])
    
    # 标注基频
    ax1.annotate(f'{f1}Hz', xy=(f1, spec_clean[int(f1 * len(clean_signal)/sample_rate)]),
                xytext=(f1-50, spec_clean[int(f1 * len(clean_signal)/sample_rate)]+0.05),
                fontsize=10, color='blue', fontweight='bold')
    ax1.annotate(f'{f2}Hz', xy=(f2, spec_clean[int(f2 * len(clean_signal)/sample_rate)]),
                xytext=(f2+20, spec_clean[int(f2 * len(clean_signal)/sample_rate)]+0.05),
                fontsize=10, color='blue', fontweight='bold')
    
    # 下半图：失真信号频谱
    # 使用条形图更清晰地显示互调产物
    ax2.bar(freq_dist[idx_dist], spec_dist[idx_dist], width=8, alpha=0.7, color='red', label='失真信号频谱')
    ax2.set_xlabel('频率 (Hz)', fontsize=12)
    ax2.set_ylabel('幅度', fontsize=12)
    ax2.set_title(f'失真信号频谱（出现新的互调产物, α={alpha}）', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, linestyle='--', axis='y')
    ax2.set_xlim([0, max_freq])
    ax2.set_ylim([0, max(spec_dist) * 1.1])
    
    # 计算并标注互调产物
    im3_1 = 2*f1 - f2  # 300 Hz
    im3_2 = 2*f2 - f1  # 1500 Hz
    im5_1 = 2*f1 + f2  # 2500 Hz
    im5_2 = 2*f2 + f1  # 2900 Hz
    
    # 查找互调产物的实际幅度
    def find_amplitude_at_freq(freq, freq_array, spec_array):
        idx = np.argmin(np.abs(freq_array - freq))
        return spec_array[idx] if idx < len(spec_array) else 0
    
    amp_im1 = find_amplitude_at_freq(im3_1, freq_dist, spec_dist)
    amp_im2 = find_amplitude_at_freq(im3_2, freq_dist, spec_dist)
    
    # 标注互调产物
    if amp_im1 > 0.01:
        ax2.annotate(f'2f₁-f₂\n{im3_1:.0f}Hz', xy=(im3_1, amp_im1),
                    xytext=(im3_1-40, amp_im1+0.08),
                    fontsize=9, color='red', fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1))
    
    if amp_im2 > 0.01:
        ax2.annotate(f'2f₂-f₁\n{im3_2:.0f}Hz', xy=(im3_2, amp_im2),
                    xytext=(im3_2+20, amp_im2+0.08),
                    fontsize=9, color='red', fontweight='bold',
                    arrowprops=dict(arrowstyle='->', color='red', lw=1))
    
    # 标注基频
    ax2.annotate(f'f₁={f1}Hz', xy=(f1, find_amplitude_at_freq(f1, freq_dist, spec_dist)),
                xytext=(f1-50, find_amplitude_at_freq(f1, freq_dist, spec_dist)+0.05),
                fontsize=10, color='blue', fontweight='bold')
    ax2.annotate(f'f₂={f2}Hz', xy=(f2, find_amplitude_at_freq(f2, freq_dist, spec_dist)),
                xytext=(f2+20, find_amplitude_at_freq(f2, freq_dist, spec_dist)+0.05),
                fontsize=10, color='blue', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('spectrum_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("✓ 频谱图已保存: spectrum_comparison.png")

def plot_combined_view(t, clean_signal, distorted_signal, sample_rate, f1, f2, alpha):
    """绘制组合视图（波形+频谱）"""
    fig = plt.figure(figsize=(16, 12))
    
    # 波形对比
    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(t, clean_signal, 'b-', linewidth=1.5, label='纯净信号')
    ax1.set_ylabel('幅度')
    ax1.set_title('纯净信号波形', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    ax1.set_xlim([t[0], t[-1]])
    
    ax2 = plt.subplot(2, 2, 2)
    ax2.plot(t, distorted_signal, 'r-', linewidth=1.5, label=f'失真信号 (α={alpha})')
    ax2.set_title('失真信号波形', fontsize=12, fontweight='bold')
    ax2.set_xlabel('时间 (秒)')
    ax2.set_ylabel('幅度')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim([t[0], t[-1]])
    
    # 频谱对比
    freq_clean, spec_clean = compute_spectrum(clean_signal, sample_rate)
    freq_dist, spec_dist = compute_spectrum(distorted_signal, sample_rate)
    
    max_freq = 2500
    idx_clean = freq_clean <= max_freq
    idx_dist = freq_dist <= max_freq
    
    ax3 = plt.subplot(2, 2, 3)
    ax3.stem(freq_clean[idx_clean], spec_clean[idx_clean], 
             basefmt=" ", linefmt='b-', markerfmt='bo', markersize=3)
    ax3.set_title('纯净信号频谱', fontsize=12, fontweight='bold')
    ax3.set_ylabel('幅度')
    ax3.set_xlabel('频率 (Hz)')
    ax3.grid(True, alpha=0.3, axis='y')
    ax3.set_xlim([0, max_freq])
    
    ax4 = plt.subplot(2, 2, 4)
    ax4.bar(freq_dist[idx_dist], spec_dist[idx_dist], width=8, alpha=0.7, color='red')
    ax4.set_title('失真信号频谱（互调产物可见）', fontsize=12, fontweight='bold')
    ax4.set_xlabel('频率 (Hz)')
    ax4.set_ylabel('幅度')
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.set_xlim([0, max_freq])
    
    # 在失真频谱上标注互调产物
    im3_1, im3_2 = 2*f1 - f2, 2*f2 - f1
    for im_freq in [im3_1, im3_2]:
        if 0 < im_freq < max_freq:
            idx = np.argmin(np.abs(freq_dist - im_freq))
            if spec_dist[idx] > 0.01:
                ax4.annotate(f'{im_freq:.0f}Hz', xy=(im_freq, spec_dist[idx]),
                            xytext=(im_freq-30, spec_dist[idx]+0.05),
                            fontsize=8, color='darkred',
                            arrowprops=dict(arrowstyle='->', color='red', lw=0.8))
    
    plt.suptitle(f'互调失真分析 (f₁={f1}Hz, f₂={f2}Hz, α={alpha})', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig('combined_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("✓ 组合分析图已保存: combined_analysis.png")

def analyze_intermodulation_products(f1, f2, alpha):
    """分析互调产物的理论值"""
    print("\n" + "="*60)
    print("互调失真理论分析")
    print("="*60)
    print(f"输入信号: {f1}Hz 和 {f2}Hz 双音信号")
    print(f"非线性系数 α = {alpha}")
    print("\n产生的频率成分:")
    print("-"*40)
    
    # 基频
    print(f"✓ 基频分量: {f1}Hz, {f2}Hz")
    
    # 二阶产物（由平方项产生）
    print(f"\n二阶互调产物 (x²项):")
    print(f"  - 和频: {f1+f2:.1f}Hz")
    print(f"  - 差频: {abs(f1-f2):.1f}Hz")
    print(f"  - 二次谐波: {2*f1:.1f}Hz, {2*f2:.1f}Hz")
    
    # 三阶产物（由立方项产生，最重要）
    print(f"\n三阶互调产物 (x³项, 最重要的失真产物):")
    print(f"  ★ 2f₁ - f₂ = {2*f1 - f2:.1f}Hz  (三阶互调，通常落在带内)")
    print(f"  ★ 2f₂ - f₁ = {2*f2 - f1:.1f}Hz  (三阶互调，通常落在带内)")
    print(f"  - 2f₁ + f₂ = {2*f1 + f2:.1f}Hz  (三阶互调，通常在带外)")
    print(f"  - 2f₂ + f₁ = {2*f2 + f1:.1f}Hz  (三阶互调，通常在带外)")
    print(f"  - 三次谐波: {3*f1:.1f}Hz, {3*f2:.1f}Hz")
    
    print("\n" + "="*60)
    print("实际影响: 2f₁-f₂ 和 2f₂-f₁ 是最有害的互调产物,")
    print("因为它们会落在工作频带内，造成邻道干扰。")
    print("="*60)

def interactive_demo():
    """交互式演示，允许用户调整参数"""
    print("\n" + "="*60)
    print("互调失真可视化程序")
    print("="*60)
    
    # 用户输入参数
    try:
        f1_input = input(f"\n请输入第一个频率 [默认 {f1}Hz]: ").strip()
        if f1_input:
            f1_val = float(f1_input)
        else:
            f1_val = f1
        
        f2_input = input(f"请输入第二个频率 [默认 {f2}Hz]: ").strip()
        if f2_input:
            f2_val = float(f2_input)
        else:
            f2_val = f2
        
        alpha_input = input(f"请输入非线性系数 [默认 {nonlinear_gain}]: ").strip()
        if alpha_input:
            alpha_val = float(alpha_input)
        else:
            alpha_val = nonlinear_gain
            
    except ValueError:
        print("输入无效，使用默认值")
        f1_val, f2_val, alpha_val = f1, f2, nonlinear_gain
    
    # 生成时间轴
    t = np.linspace(0, DURATION, NUM_SAMPLES)
    
    # 生成信号
    clean_signal = generate_clean_signal(f1_val, f2_val, amplitude, t)
    distorted_signal = generate_distorted_signal(f1_val, f2_val, amplitude, alpha_val, t)
    
    # 分析理论互调产物
    analyze_intermodulation_products(f1_val, f2_val, alpha_val)
    
    # 生成图像
    print("\n正在生成图像...")
    plot_waveform(t, clean_signal, distorted_signal, f1_val, f2_val, alpha_val)
    plot_spectrum(clean_signal, distorted_signal, SAMPLE_RATE, f1_val, f2_val, alpha_val)
    plot_combined_view(t, clean_signal, distorted_signal, SAMPLE_RATE, f1_val, f2_val, alpha_val)
    
    print("\n✨ 所有图像生成完成！")
    print("生成的图片文件:")
    print("  1. waveform_comparison.png - 时域波形对比")
    print("  2. spectrum_comparison.png - 频谱对比")
    print("  3. combined_analysis.png - 综合分析图")

if __name__ == "__main__":
    # 运行交互式演示
    interactive_demo()
    
    # 可选：直接运行默认参数版本
    # t = np.linspace(0, DURATION, NUM_SAMPLES)
    # clean_signal = generate_clean_signal(f1, f2, amplitude, t)
    # distorted_signal = generate_distorted_signal(f1, f2, amplitude, nonlinear_gain, t)
    # analyze_intermodulation_products(f1, f2, nonlinear_gain)
    # plot_waveform(t, clean_signal, distorted_signal, f1, f2, nonlinear_gain)
    # plot_spectrum(clean_signal, distorted_signal, SAMPLE_RATE, f1, f2, nonlinear_gain)
    # plot_combined_view(t, clean_signal, distorted_signal, SAMPLE_RATE, f1, f2, nonlinear_gain)