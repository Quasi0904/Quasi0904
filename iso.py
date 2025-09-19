import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import emcee
from scipy.spatial import cKDTree
from scipy.interpolate import Akima1DInterpolator
from matplotlib.ticker import AutoMinorLocator
import os, re
import corner
import warnings
warnings.filterwarnings('ignore')

# ---------------- Z与[Fe/H]转换 ----------------
def feh_to_z(feh, z_sun=0.0134):
    return z_sun * 10**feh

def av_to_ebv(av, rv=3.1):
    return av / rv

# ---------------- 数据读取 ----------------
def load_observation_data(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True, comment='#',
                     names=['RA','DEC','MAG_F475W','ERR_F475W', 
                            'MAG_F814W','ERR_F814W'])
    df['Color'] = df['MAG_F475W'] - df['MAG_F814W']
    mask = (np.isfinite(df['Color']) & (df['MAG_F475W']<99) & (df['MAG_F814W']<99) &
            (df['ERR_F475W']<1) & (df['ERR_F814W']<1))
    return {
        'color': df['Color'][mask].values,
        'mag': df['MAG_F814W'][mask].values,
        'color_err': np.sqrt(df['ERR_F475W'][mask]**2 + df['ERR_F814W'][mask]**2),
        'mag_err': df['ERR_F814W'][mask].values
    }

# ---------------- MIST加载与过滤 ----------------
def parse_feh(fname):
    m = re.search(r'feh_([mp])(\d+\.?\d*)', fname)
    return -float(m.group(2)) if m.group(1)=='m' else float(m.group(2)) if m else 0.0

def load_all_isochrones(mist_files):
    grid = []
    for f in mist_files:
        feh = parse_feh(os.path.basename(f))
        lines = open(f).readlines()
        for i, l in enumerate(lines):
            if 'log10_isochrone_age_yr' in l:
                start = i + 1
                break
        data = [list(map(float, l.split()[:23])) for l in lines[start:] if l.strip() and not l.startswith('#')]
        cols = ['EEP','log_age','IM','SM','logTe','logg','logL',
                'FeH_init','FeH','F435W','F475W','F502N','F550M','F555W','F606W',
                'F625W','F658N','F660N','F775W','F814W','F850LP','F892N','phase']
        df = pd.DataFrame(data, columns=cols)
        for age in df['log_age'].unique():
            s = df[df['log_age'] == age]
            eep = s['EEP'].values
            f475 = s['F475W'].values
            f814 = s['F814W'].values
            # 严格过滤：无NaN，EEP唯一且≥5个点
            if (not np.any(np.isnan(f475)) and not np.any(np.isnan(f814)) and
                len(np.unique(eep)) >= 5):
                grid.append({
                    'log_age': age,
                    'FeH': feh,
                    'EEP': eep,
                    'F475W': f475,
                    'F814W': f814
                })
    return grid

# ---------------- EEP对齐+Akima插值（安全版） ----------------
def interpolate_isochrone(grid, log_age, FeH, n_pts=400):
    ages = sorted(set(g['log_age'] for g in grid))
    fehs = sorted(set(g['FeH'] for g in grid))
    
    # 检查是否在网格范围内
    if not (ages[0] <= log_age <= ages[-1] and fehs[0] <= FeH <= fehs[-1]):
        return np.array([]), np.array([])
    
    a0 = max(a for a in ages if a <= log_age)
    a1 = min(a for a in ages if a >= log_age)
    f0 = max(f for f in fehs if f <= FeH)
    f1 = min(f for f in fehs if f >= FeH)
    combos = [(a0, f0), (a0, f1), (a1, f0), (a1, f1)]
    
    EEP_common = np.linspace(200, 800, n_pts)
    interp_data = []
    
    for (a, f) in combos:
        iso = next((g for g in grid if g['log_age'] == a and g['FeH'] == f), None)
        if iso is None:
            return np.array([]), np.array([])
        # 再次检查（以防过滤遗漏，但正常load_all_isochrones已处理）
        if (len(np.unique(iso['EEP'])) < 5 or
            np.any(np.isnan(iso['F475W'])) or np.any(np.isnan(iso['F814W']))):
            return np.array([]), np.array([])
        # 插值
        f475_fun = Akima1DInterpolator(iso['EEP'], iso['F475W'])
        f814_fun = Akima1DInterpolator(iso['EEP'], iso['F814W'])
        interp_data.append((f475_fun(EEP_common), f814_fun(EEP_common)))
    
    # 双线性插值
    tx = (log_age - a0) / (a1 - a0 + 1e-9)
    ty = (FeH - f0) / (f1 - f0 + 1e-9)
    F475_interp = (
        (1 - tx) * (1 - ty) * interp_data[0][0] +
        (1 - tx) * ty * interp_data[1][0] +
        tx * (1 - ty) * interp_data[2][0] +
        tx * ty * interp_data[3][0]
    )
    F814_interp = (
        (1 - tx) * (1 - ty) * interp_data[0][1] +
        (1 - tx) * ty * interp_data[1][1] +
        tx * (1 - ty) * interp_data[2][1] +
        tx * ty * interp_data[3][1]
    )
    
    # 最终检查：无NaN/空
    if (F475_interp.size == 0 or F814_interp.size == 0 or
        np.any(~np.isfinite(F475_interp)) or np.any(~np.isfinite(F814_interp))):
        return np.array([]), np.array([])
    
    return F475_interp, F814_interp

# ---------------- 消光&距离 ----------------
def apply_reddening(F475, F814, Av):
    A475, A814 = Av * 1.21179, Av * 0.59696
    app475 = F475 + A475 + 5 * np.log10(780000) - 5  # DM=5*log10(780000)-5
    app814 = F814 + A814 + 5 * np.log10(780000) - 5
    return app475 - app814, app814
    
def av_to_ebv_final(av):
    return av / 3.634
# ---------------- chi2（防御版） ----------------
def chi2(obs, model_color, model_mag):
    model_color = np.asarray(model_color)
    model_mag = np.asarray(model_mag)
    
    # 严格检查
    if (model_color.size == 0 or model_mag.size == 0 or
        np.any(~np.isfinite(model_color)) or np.any(~np.isfinite(model_mag))):
        return np.inf
    
    tree = cKDTree(np.c_[model_color, model_mag])
    _, idx = tree.query(np.c_[obs['color'], obs['mag']], k=1)
    
    color_err = np.maximum(obs['color_err'], 0.01)
    mag_err = np.maximum(obs['mag_err'], 0.01)
    
    return np.sum(((model_color[idx] - obs['color']) / color_err)**2 +
                  ((model_mag[idx] - obs['mag']) / mag_err)**2)

# ---------------- MCMC模型 ----------------
class IsoMCMC:
    def __init__(self, obs, grid):
        self.obs = obs
        self.grid = grid
    
    def ln_prior(self, theta):
        age, feh, Av = theta
        grid_ages = sorted(set(g['log_age'] for g in self.grid))
        grid_fehs = sorted(set(g['FeH'] for g in self.grid))
        if (grid_ages[0] <= age <= grid_ages[-1] and
            grid_fehs[0] <= feh <= grid_fehs[-1] and
            0 <= Av < 3):
            return 0.0
        return -np.inf
    
    def ln_likelihood(self, theta):
        age, feh, Av = theta
        F475, F814 = interpolate_isochrone(self.grid, age, feh)
        if F475.size == 0 or F814.size == 0:
            return -np.inf
        c, m = apply_reddening(F475, F814, Av)
        return -0.5 * chi2(self.obs, c, m)
    
    def ln_prob(self, theta):
        lp = self.ln_prior(theta)
        return lp + self.ln_likelihood(theta) if np.isfinite(lp) else -np.inf

# ---------------- 绘制CMD ----------------
def plot_cmd(obs, model_c, model_m, best, best_err, APID):
    color = np.asarray(obs.get('color', []))
    mag814 = np.asarray(obs.get('mag', []))

    # 如果观测数据为空，给出提示并尝试用模型确定坐标范围或直接跳过绘图
    if color.size == 0 or mag814.size == 0:
        print(f"[Warning] No observed stars for APID {APID}. Skipping observed CMD plotting.")
        # 如果有模型等龄线，则用模型范围绘制等龄线
        if model_c.size > 0 and model_m.size > 0:
            x_min, x_max = np.nanmin(model_c), np.nanmax(model_c)
            y_min, y_max = np.nanmin(model_m), np.nanmax(model_m)
        else:
            print("[Info] No model isochrone either. Nothing to plot.")
            return
    else:
        # 正常情况：根据观测数据计算范围（防 NaN）
        # 使用 nanmin/nanmax 更稳健
        x_min, x_max = np.nanmin(color), np.nanmax(color)
        y_min, y_max = np.nanmin(mag814), np.nanmax(mag814)
        x_range = x_max - x_min if (x_max > x_min) else 0.5
        y_range = y_max - y_min if (y_max > y_min) else 0.5
        x_margin = x_range * 0.1
        y_margin = y_range * 0.1
        x_min -= x_margin
        x_max += x_margin
        y_min -= y_margin
        y_max += y_margin  

    # 创建图形与坐标轴
    fig, ax = plt.subplots(figsize=(8, 10), dpi=150)

    # 绘制简单圆圈数据点
    ax.scatter(color, mag814, 
               s=15,              
               edgecolors='black',
               facecolors='none', 
               alpha=0.7,
               linewidths=0.5)

    # 绘制最佳拟合等龄线
    ax.plot(model_c, model_m, 'r-', lw=2, label='Best-fit Isochrone')

    # 设置坐标范围
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # 反转Y轴
    ax.invert_yaxis()

    # 设置坐标轴标签与标题
    ax.set_xlabel('F475W - F814W', fontsize=12, fontweight='bold')
    ax.set_ylabel('F814W', fontsize=12, fontweight='bold')
    ax.set_title(f'J{APID} Color-Magnitude Diagram', fontsize=14, fontweight='bold', pad=20)

    # 主次刻度与网格（科研风格）
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which='major', alpha=0.5, linestyle='--')
    ax.grid(which='minor', alpha=0.3, linestyle=':')

    # 计算物理量
    z_value = feh_to_z(best[1])
    ebv_value = av_to_ebv(best[2])

    # 左下角文本框（最佳拟合参数）
    textstr = '\n'.join((
        f'log(age) = {best[0]:.3f} ± {best_err[0]:.3f}',
        f'Z = {z_value:.4f} ± {feh_to_z(best_err[1]):.4f}',
        f'E(B-V) = {ebv_value:.3f} ± {av_to_ebv(best_err[2]):.3f}'))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.05, 0.05, textstr, transform=ax.transAxes, fontsize=10,
            verticalalignment='bottom', bbox=props)

    ax.legend()

    # 保存并显示
    plt.savefig(f"/home/quasi/桌面/m31/JHONSON/J{APID}/{APID}_CMD_fit.png",
                dpi=300, bbox_inches='tight')
    plt.show()

# ---------------- 主程序 ----------------
def main():
    APID=input("APID: ").strip()
    obs_file=f"/home/quasi/桌面/m31/JHONSON/J{APID}/{APID}.txt"
    obs=load_observation_data(obs_file)
    mist_dir='/home/quasi/桌面/m31/MIST_v1.2_vvcrit0.0_HST_ACSWF'
    mist_files=glob.glob(os.path.join(mist_dir,'*.iso.cmd'))
    grid=load_all_isochrones(mist_files)
    model=IsoMCMC(obs, grid)
    ndim,nwalkers=3,64
    p0=np.array([9.0,0.0,0.5])+1e-3*np.random.randn(nwalkers,ndim)
    sampler=emcee.EnsembleSampler(nwalkers, ndim, model.ln_prob)
    print("Burn-in...")
    p0,_,_=sampler.run_mcmc(p0, 3000, progress=True); sampler.reset()
    print("Production...")
    sampler.run_mcmc(p0, 8000, progress=True)
    tau=sampler.get_autocorr_time(tol=0)
    print("Autocorr times:",tau)
    samples=sampler.get_chain(flat=True)
    best=np.median(samples,axis=0)
    q16,q84=np.percentile(samples,[16,84],axis=0)
    best_err=(q84-q16)/2
    print(f"Best-fit: log_age={best[0]:.3f}±{best_err[0]:.3f}, [Fe/H]={best[1]:.3f}±{best_err[1]:.3f}, Av={best[2]:.3f}±{best_err[2]:.3f}")
    print(f"Derived: Z={feh_to_z(best[1]):.6f}, E(B-V)={av_to_ebv_final(best[2]):.3f}")
    corner.corner(samples, labels=["log_age","[Fe/H]","Av"], truths=best)
    plt.savefig(f"/home/quasi/桌面/m31/JHONSON/J{APID}/{APID}_corner.png", dpi=300, bbox_inches='tight')
    plt.show()
    F475,F814=interpolate_isochrone(grid, best[0], best[1])
    c,m=apply_reddening(F475, F814, best[2])
    plot_cmd(obs, c, m, best, best_err, APID)

if __name__=="__main__":
    main()

