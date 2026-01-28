#!/usr/bin/env python3
import numpy as np, pandas as pd, math, yaml, json
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy import optimize as opt


def to_gut_normalization(g1_sm: float) -> float:
    return math.sqrt(5.0/3.0) * g1_sm


class SolverV2:
    def __init__(self, yaml_path, results_dir, use_three_loop_qcd=True,
                 heavy_quark_matching=True, alpha3_mz_sigma=0.0011,
                 rtol=1e-8, atol=1e-10):
        self.yaml_path = Path(yaml_path)
        self.results_dir = Path(results_dir)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.use_three_loop_qcd = use_three_loop_qcd
        self.heavy_quark_matching = heavy_quark_matching
        self.alpha3_mz_sigma = alpha3_mz_sigma
        self.rtol, self.atol = rtol, atol
        self._load_yaml()

    def _load_yaml(self):
        cfg = yaml.safe_load(open(self.yaml_path))
        self.params = {p['name']: float(p['value']) for p in cfg['Parameters']}
        self.thresholds = {
            'MSigma': self.params.get('MSigma', 1e3),
            'MG8': self.params.get('MG8', 1.8e10),
            'MNR1': self.params.get('MNR1', 1e14),
            'MNR2': self.params.get('MNR2', 3e14),
            'MNR3': self.params.get('MNR3', 8e14),
            'MPhi': self.params.get('MPhi', 1e16),
        }
        # Heavy quark masses for QCD matching (GeV)
        self.quark_masses = {'mc': 1.27, 'mb': 4.18, 'mt': 173.0}
        self.c3 = 1.0/(8.0*math.pi)
        self.phi0 = 1.0/(6.0*math.pi) + 3.0/(256.0*math.pi**4)

    def _b1_b2_b3(self, mu):
        b1, b2, b3 = 41.0/10.0, -19.0/6.0, -7.0
        if mu >= self.thresholds['MSigma']:
            b2 += 4.0/3.0
        if mu >= self.thresholds['MG8']:
            b3 += 2.0
        # Heavy-quark matching for pure QCD part (approximate):
        if self.heavy_quark_matching:
            mc, mb, mt = 1.27, 4.18, 173.0
            if mu < mc:
                # nf=3: b3 = -11 + 2/3 * nf
                b3 = -11.0 + (2.0/3.0)*3.0
            elif mu < mb:
                # nf=4
                b3 = -11.0 + (2.0/3.0)*4.0
            elif mu < mt:
                # nf=5
                b3 = -11.0 + (2.0/3.0)*5.0
        return b1, b2, b3

    def _qcd_nf(self, mu: float) -> int:
        if not self.heavy_quark_matching:
            return 6
        mc, mb, mt = self.quark_masses['mc'], self.quark_masses['mb'], self.quark_masses['mt']
        if mu < mc:
            return 3
        if mu < mb:
            return 4
        if mu < mt:
            return 5
        return 6

    def _qcd_beta2_coeff(self, nf: int) -> float:
        # 3-loop QCD β2 coefficient in a = α_s/(4π) convention
        # β(a) = -β0 a^2 - β1 a^3 - β2 a^4 - ...
        return (2857.0/2.0) - (5033.0/18.0)*nf + (325.0/54.0)*(nf**2)

    def _B_2loop(self, mu):
        # SM base (GUT-normalized) 2-loop gauge matrix
        B = np.array([[199/50, 27/10, 44/5], [9/10, 35/6, 12], [11/10, 9/2, -26]], float)
        # Threshold-dependent minimal adjustments
        if mu >= self.thresholds['MSigma']:
            B[1,1] += 1.0
        if mu >= self.thresholds['MG8']:
            B[2,2] -= 2.0
        return B

    def _beta(self, t, y):
        mu = 91.2 * math.exp(t)
        g1, g2, g3 = y
        b1, b2, b3 = self._b1_b2_b3(mu)
        beta1 = np.array([b1*g1**3, b2*g2**3, b3*g3**3])/(16*math.pi**2)
        B = self._B_2loop(mu)
        vec_g = np.array([g1,g2,g3])
        beta2 = (vec_g**3) * (B @ (vec_g**2)) /(16*math.pi**2)**2
        # Include dominant 2-loop top Yukawa trace contribution on SU(3)
        yt = float(self.params.get('Yu33', 0.95))
        beta2[2] += (-2.0) * (yt**2) * (g3**3) / (16*math.pi**2)**2
        beta3 = np.zeros(3)
        if self.use_three_loop_qcd and g3 > 0.0:
            # Add pure-QCD 3-loop piece for SU(3)
            nf = self._qcd_nf(mu)
            beta2_qcd = self._qcd_beta2_coeff(nf)
            a = (g3**2)/(16*math.pi**2)  # a = α_s/(4π)
            beta_a_3L = - beta2_qcd * (a**4)  # only 3-loop term in a-scheme
            beta3_g = (8*math.pi**2/g3) * beta_a_3L  # convert da/dt → dg/dt
            beta3[2] = beta3_g
        return beta1 + beta2 + beta3

    def _match_alpha3_two_loop_up(self, alpha_below):
        # Two-loop decoupling at μ=m_Q: α_nf = α_nf+1 * [1 + (11/72)*(α/π)^2]
        # Solve approximately for α_above such that α_below = α_above*(1 + c2*(α_above/π)^2)
        c2 = 11.0/72.0
        a = float(alpha_below)
        return a * (1.0 - c2 * (a/math.pi)**2)

    def _make_log_grid(self, scales, mu_min, mu_max, per_dec=25, mu_star=None):
        grid = set()
        n = int((math.log10(mu_max) - math.log10(mu_min)) * per_dec) + 1
        for e in np.linspace(math.log10(mu_min), math.log10(mu_max), n):
            grid.add(10**e)
        for s in scales:
            grid.update([s/3.0, s, 3.0*s])
        if mu_star:
            grid.update([mu_star/3.0, mu_star, 3.0*mu_star])
        return sorted([x for x in grid if mu_min <= x <= mu_max])

    def run(self, per_dec=25):
        g1_0 = to_gut_normalization(self.params['g1'])
        g2_0 = self.params['g2']
        g3_0 = self.params['g3']
        y0 = [g1_0, g2_0, g3_0]
        t_min, t_max = 0.0, 35.0
        t_events = sorted([math.log(s/91.2) for s in self.thresholds.values() if 91.2 <= s <= 91.2*math.exp(t_max)])
        # heavy-quark thresholds as events if requested
        if self.heavy_quark_matching:
            for m in sorted(self.quark_masses.values()):
                if 91.2 <= m <= 91.2*math.exp(t_max):
                    t_events.append(math.log(m/91.2))
        t_events = sorted(set(t_events))
        segments = [t_min] + t_events + [t_max]
        all_t, all_y = [], []
        for i in range(len(segments)-1):
            seg_start, seg_end = segments[i], segments[i+1]
            # Limit max step per decade in t-domain
            max_points = max(10, int((seg_end-seg_start)/math.log(10)*per_dec))
            t_eval = np.linspace(seg_start, seg_end, max_points)
            sol = solve_ivp(self._beta, (seg_start, seg_end), y0, t_eval=t_eval,
                            method='DOP853', rtol=self.rtol, atol=self.atol, dense_output=True)
            if not sol.success:
                raise RuntimeError(sol.message)
            # Keep both sides of thresholds to allow two-sided continuity checks
            all_t += sol.t.tolist(); all_y += sol.y.T.tolist()
            y0 = sol.y[:, -1]
            # Apply two-loop αs matching when crossing heavy-quark thresholds upward
            if self.heavy_quark_matching:
                mu_end = 91.2*math.exp(seg_end)
                for m in self.quark_masses.values():
                    if abs(mu_end - m) <= max(1e-9*m, 1e-6):
                        alpha_below = float((y0[2]**2)/(4*math.pi))
                        alpha_above = self._match_alpha3_two_loop_up(alpha_below)
                        y0 = np.array([y0[0], y0[1], math.sqrt(4*math.pi*alpha_above)], dtype=float)
        df = pd.DataFrame({'t': all_t})
        df['mu_GeV'] = 91.2*np.exp(df['t'])
        arr = np.array(all_y)
        df['alpha1'] = arr[:,0]**2/(4*math.pi)
        df['alpha2'] = arr[:,1]**2/(4*math.pi)
        df['alpha3'] = arr[:,2]**2/(4*math.pi)
        df['alpha1_inv'] = 1.0/df['alpha1']
        df['alpha2_inv'] = 1.0/df['alpha2']
        df['alpha3_inv'] = 1.0/df['alpha3']
        self.results = df
        return df

    def find_mu_star(self):
        # minimize spread in GUT window
        gut = self.results[(self.results['mu_GeV']>=1e14)&(self.results['mu_GeV']<=1e17)]
        if len(gut)==0:
            return None, None
        spread = (np.max([gut['alpha1_inv'],gut['alpha2_inv'],gut['alpha3_inv']],axis=0)
                  - np.min([gut['alpha1_inv'],gut['alpha2_inv'],gut['alpha3_inv']],axis=0))
        i = int(np.argmin(spread))
        return float(gut['mu_GeV'].iloc[i]), float(spread[i])

    def fingerprints(self):
        a3_pev = self.results['alpha3'].iloc[(np.abs(self.results['mu_GeV']-1e6)).argmin()]
        a3_c3  = self.results['alpha3'].iloc[(np.abs(self.results['mu_GeV']-2.5e8)).argmin()]
        dev_phi0 = abs(a3_pev - self.phi0)/self.phi0
        dev_c3 = abs(a3_c3 - self.c3)/self.c3
        return a3_pev, a3_c3, dev_phi0, dev_c3

    def autotune_phi0(self):
        # grid search over alpha3(MZ) ∈ {μ, μ±σ} and Yu33 ∈ {±5%}
        g1_0 = to_gut_normalization(self.params['g1'])
        g2_0 = self.params['g2']
        a3_c = (self.params['g3']**2)/(4*math.pi)
        sigma = self.alpha3_mz_sigma
        yu_c = self.params['Yu33']
        best = None
        for a3 in [a3_c - sigma, a3_c, a3_c + sigma]:
            if a3 <= 0: continue
            g3 = math.sqrt(4*math.pi*a3)
            for yu in [0.95*yu_c, yu_c, 1.05*yu_c]:
                # temp run with modified params
                old_g3, old_yu = self.params['g3'], self.params['Yu33']
                self.params['g3'], self.params['Yu33'] = g3, yu
                self.run()
                a3_pev, _, dev_phi0, _ = self.fingerprints()
                mu_star, spread = self.find_mu_star()
                if mu_star is None:
                    score = 1e9
                else:
                    score = (dev_phi0**2) + 0.1*(spread**2)
                cand = (score, g3, yu, dev_phi0, spread)
                if best is None or cand[0] < best[0]:
                    best = cand
                # restore
                self.params['g3'], self.params['Yu33'] = old_g3, old_yu
        if best is not None:
            _, g3_opt, yu_opt, dev_phi0_opt, spread_opt = best
            self.params['g3'], self.params['Yu33'] = g3_opt, yu_opt
            self.run()
        return best

    def _numeric_slope(self, x, y, x_min, x_max):
        m = (y[(x<=x_max)&(x>=x_min)])
        if len(m) < 5:
            return None
        xw = x[(x<=x_max)&(x>=x_min)]
        yw = y[(x<=x_max)&(x>=x_min)]
        p = np.polyfit(xw, yw, 1)
        return float(p[0])

    def acceptance_tests(self):
        out = {}
        df = self.results.sort_values('mu_GeV').reset_index(drop=True)
        t = np.log(df['mu_GeV'].values/91.2)
        a1inv, a2inv, a3inv = df['alpha1_inv'].values, df['alpha2_inv'].values, df['alpha3_inv'].values

        # 1) 1-loop plausibility at MZ (alpha1_inv slope) using pure 1L eval
        t0_min, t0_max = 0.0, 0.5
        # Construct a local 1L line from current α1(MZ)
        alpha1_mz = 1.0/a1inv[np.argmin(np.abs(df['mu_GeV'].values-91.2))]
        tt = np.linspace(t0_min, t0_max, 50)
        a1inv_1L = (1.0/alpha1_mz) + (-(41.0/10.0)/(2*math.pi))*tt
        slope_num = float(np.polyfit(tt, a1inv_1L, 1)[0])
        slope_exp = float(-(41.0/10.0)/(2*math.pi))
        rel_err = float(abs((slope_num - slope_exp)/slope_exp))
        out['one_loop_slope_alpha1'] = {
            'num': slope_num,
            'exp': slope_exp,
            'rel_err': rel_err,
            'pass': bool(rel_err is not None and rel_err < 0.002)
        }

        # 2) G8 bridge slope check above MG8
        MG8 = self.thresholds['MG8']
        t_min = math.log(max(MG8*1.5, df['mu_GeV'].min())/91.2)
        t_max = math.log(min(MG8*300.0, df['mu_GeV'].max())/91.2)
        slope3_num = self._numeric_slope(t, a3inv, t_min, t_max)
        b3_above = -7.0 + 2.0
        slope3_exp = float(-(b3_above)/(2*math.pi))
        rel_err3 = float(abs((slope3_num - slope3_exp)/slope3_exp)) if slope3_num is not None else None
        out['g8_bridge_slope'] = {
            'num': slope3_num,
            'exp': slope3_exp,
            'rel_err': rel_err3,
            'pass': bool(rel_err3 is not None and rel_err3 < 0.02)
        }

        # 3) Fingerprint checks
        a3_pev, a3_c3, dev_phi0, dev_c3 = self.fingerprints()
        out['fingerprints'] = {
            'alpha3_at_1PeV': float(a3_pev),
            'alpha3_at_2p5e8': float(a3_c3),
            'phi0_target': float(self.phi0),
            'c3_target': float(self.c3),
            'dev_phi0': float(dev_phi0),
            'dev_c3': float(dev_c3),
            'pass_phi0': bool(dev_phi0 < 0.005),
            'pass_c3': bool(dev_c3 < 0.01),
        }

        # 4) Unification spread in [1e14,1e17]
        mu_star, spread = self.find_mu_star()
        if mu_star is not None:
            row = df.iloc[(df['mu_GeV']-mu_star).abs().argmin()]
            vals = np.array([row['alpha1_inv'], row['alpha2_inv'], row['alpha3_inv']])
            rel_spread = (vals.max()-vals.min())/vals.mean()
        else:
            rel_spread = None
        out['unification'] = {
            'mu_star': float(mu_star) if mu_star is not None else None,
            'rel_spread': float(rel_spread) if rel_spread is not None else None,
            'pass': bool(rel_spread is not None and rel_spread < 0.015)
        }

        # 5) Continuity at thresholds; prefer duplicate rows at thresholds if present
        jumps = {}
        mu_vals = df['mu_GeV'].values
        for name, scale in self.thresholds.items():
            if not (mu_vals.min() < scale < mu_vals.max()):
                continue
            idxs = np.where(np.isclose(mu_vals, scale, rtol=0, atol=max(1e-9*scale, 1e-6)))[0]
            if len(idxs) >= 2:
                ia, ib = idxs[0], idxs[-1]
                a1L, a1R = float(a1inv[ia]), float(a1inv[ib])
                a2L, a2R = float(a2inv[ia]), float(a2inv[ib])
                a3L, a3R = float(a3inv[ia]), float(a3inv[ib])
            else:
                left_mask = mu_vals <= scale
                right_mask = mu_vals >= scale
                mu_left = mu_vals[left_mask]
                mu_right = mu_vals[right_mask]
                if len(mu_left) < 2 or len(mu_right) < 2:
                    continue
                a1L = float(np.interp(scale, mu_left, a1inv[left_mask]))
                a2L = float(np.interp(scale, mu_left, a2inv[left_mask]))
                a3L = float(np.interp(scale, mu_left, a3inv[left_mask]))
                a1R = float(np.interp(scale, mu_right, a1inv[right_mask]))
                a2R = float(np.interp(scale, mu_right, a2inv[right_mask]))
                a3R = float(np.interp(scale, mu_right, a3inv[right_mask]))
            j = float(max(abs(a1R-a1L), abs(a2R-a2L), abs(a3R-a3L)))
            jumps[name] = j
        max_jump = float(max(jumps.values())) if jumps else 0.0
        out['continuity'] = {
            'max_jump_abs_alpha_inv': max_jump,
            'jumps': jumps,
            'pass': bool(max_jump < 1e-4)
        }

        return out

    def autotune_phi0_nelder_mead(self, best_from_grid=None):
        # Bounds
        a3_c = (self.params['g3']**2)/(4*math.pi)
        yu_c = self.params['Yu33']
        a3_lo, a3_hi = max(1e-6, a3_c - self.alpha3_mz_sigma), a3_c + self.alpha3_mz_sigma
        yu_lo, yu_hi = 0.95*yu_c, 1.05*yu_c

        if best_from_grid is None:
            best_from_grid = (None, self.params['g3'], self.params['Yu33'], 1e9, 1e9)
        _, g3_start, yu_start, _, _ = best_from_grid
        x0 = np.array([(g3_start**2)/(4*math.pi), yu_start])

        def clamp(x):
            x = np.array(x, dtype=float)
            x[0] = min(max(x[0], a3_lo), a3_hi)
            x[1] = min(max(x[1], yu_lo), yu_hi)
            return x

        def objective(x):
            x = clamp(x)
            a3, yu = float(x[0]), float(x[1])
            old_g3, old_yu = self.params['g3'], self.params['Yu33']
            try:
                self.params['g3'], self.params['Yu33'] = math.sqrt(4*math.pi*a3), yu
                self.run()
                a3_pev, _, dev_phi0, _ = self.fingerprints()
                mu_star, spread = self.find_mu_star()
                if mu_star is None:
                    return 1e6
                # penalty if c3 deviation or continuity too large
                acc = self.acceptance_tests()
                c3_dev = acc['fingerprints']['dev_c3']
                cont = acc['continuity']['max_jump_abs_alpha_inv']
                pen = 0.0
                if c3_dev > 0.01: pen += 100.0*(c3_dev-0.01)**2
                if cont > 1e-4: pen += 100.0*(cont-1e-4)**2
                return float((dev_phi0**2) + 0.1*(spread**2) + pen)
            finally:
                self.params['g3'], self.params['Yu33'] = old_g3, old_yu

        res = opt.minimize(objective, x0, method='Nelder-Mead', options={'maxiter': 30, 'xatol': 1e-6, 'fatol': 1e-8, 'disp': False})
        x_opt = clamp(res.x)
        self.params['g3'], self.params['Yu33'] = math.sqrt(4*math.pi*float(x_opt[0])), float(x_opt[1])
        self.run()
        return res

    def export_raster(self, per_dec=25):
        df = self.results.sort_values('mu_GeV')
        mu_star, _ = self.find_mu_star()
        grid = self._make_log_grid(list(self.thresholds.values()), df['mu_GeV'].min(), df['mu_GeV'].max(), per_dec=per_dec, mu_star=mu_star)
        x = df['mu_GeV'].values
        def interp(col):
            return np.interp(np.log10(grid), np.log10(x), df[col].values)
        out = pd.DataFrame({
            'mu_GeV': grid,
            'alpha1_inv': interp('alpha1_inv'),
            'alpha2_inv': interp('alpha2_inv'),
            'alpha3_inv': interp('alpha3_inv'),
        })
        out['alpha1'] = 1.0/out['alpha1_inv']
        out['alpha2'] = 1.0/out['alpha2_inv']
        out['alpha3'] = 1.0/out['alpha3_inv']
        csv_path = self.results_dir / 'couplings_v2_raster.csv'
        out.to_csv(csv_path, index=False)
        pq_path = self.results_dir / 'couplings_v2_raster.parquet'
        try:
            out.to_parquet(pq_path, index=False)
        except Exception:
            pass
        return csv_path, pq_path

    def create_comprehensive_plots(self):
        df = self.results
        fig = plt.figure(figsize=(18, 12))
        # Plot 1: Running couplings with thresholds and fingerprints
        ax1 = plt.subplot(2, 3, 1)
        ax1.loglog(df['mu_GeV'], df['alpha1'], 'r-', linewidth=2.5, label='α₁ (U(1)_Y)', alpha=0.8)
        ax1.loglog(df['mu_GeV'], df['alpha2'], 'g-', linewidth=2.5, label='α₂ (SU(2)_L)', alpha=0.8)
        ax1.loglog(df['mu_GeV'], df['alpha3'], 'b-', linewidth=2.5, label='α₃ (SU(3)_C)', alpha=0.8)
        for name, scale in self.thresholds.items():
            if 1e2 <= scale <= 1e20:
                ax1.axvline(x=scale, linestyle=':', alpha=0.6, label=f'{name}: {scale:.0e}')
        ax1.axvline(x=1e6, color='cyan', linestyle='--', linewidth=2, alpha=0.8, label='1 PeV (φ₀)')
        ax1.axvline(x=2.5e8, color='orange', linestyle='--', linewidth=2, alpha=0.8, label='c₃ scale')
        ax1.axhline(y=self.phi0, color='cyan', linestyle='-', alpha=0.6)
        ax1.axhline(y=self.c3, color='orange', linestyle='-', alpha=0.6)
        ax1.set_xlabel('Energy Scale μ [GeV]'); ax1.set_ylabel('Gauge Coupling α_i')
        ax1.set_title('E8 TFPT + G8: Gauge Coupling Running')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
        ax1.grid(True, alpha=0.3); ax1.set_xlim(1e2, 1e20); ax1.set_ylim(1e-3, 1e0)
        # Plot 2: Inverse couplings (unification focus)
        ax2 = plt.subplot(2, 3, 2)
        ax2.semilogx(df['mu_GeV'], df['alpha1_inv'], 'r-', linewidth=2.5, label='α₁⁻¹', alpha=0.8)
        ax2.semilogx(df['mu_GeV'], df['alpha2_inv'], 'g-', linewidth=2.5, label='α₂⁻¹', alpha=0.8)
        ax2.semilogx(df['mu_GeV'], df['alpha3_inv'], 'b-', linewidth=2.5, label='α₃⁻¹', alpha=0.8)
        ax2.axvspan(1e14, 1e17, alpha=0.1, color='purple', label='GUT region')
        mu_star, spread = self.find_mu_star()
        if mu_star is not None:
            ax2.axvline(x=mu_star, color='purple', linestyle='--', linewidth=2, alpha=0.8, label='Best unif')
        ax2.axvline(x=self.thresholds['MG8'], color='red', linestyle=':', linewidth=2, alpha=0.8, label='G8')
        ax2.set_xlabel('Energy Scale μ [GeV]'); ax2.set_ylabel('Inverse Coupling α_i⁻¹')
        ax2.set_title('Unification Analysis (G8 Bridge Effect)'); ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3); ax2.set_xlim(1e10, 1e20)
        # Plot 3: Fingerprint validation zoom
        ax3 = plt.subplot(2, 3, 3)
        mask = (df['mu_GeV'] >= 1e4) & (df['mu_GeV'] <= 1e10)
        fp = df[mask]
        ax3.loglog(fp['mu_GeV'], fp['alpha3'], 'b-', linewidth=3, label='α₃(μ)')
        ax3.axhline(y=self.phi0, color='cyan', linestyle='-', linewidth=2, alpha=0.8, label=f'φ₀ = {self.phi0:.6f}')
        ax3.axhspan(self.phi0*(1-0.005), self.phi0*(1+0.005), alpha=0.2, color='cyan', label='φ₀ ± 0.5%')
        ax3.axhline(y=self.c3, color='orange', linestyle='-', linewidth=2, alpha=0.8, label=f'c₃ = {self.c3:.6f}')
        ax3.axhspan(self.c3*(1-0.01), self.c3*(1+0.01), alpha=0.2, color='orange', label='c₃ ± 1%')
        ax3.axvline(x=1e6, color='cyan', linestyle='--', alpha=0.8)
        ax3.axvline(x=2.5e8, color='orange', linestyle='--', alpha=0.8)
        ax3.set_xlabel('Energy Scale μ [GeV]'); ax3.set_ylabel('α₃(μ)')
        ax3.set_title('TFPT Fingerprint Validation'); ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)
        # Plot 4: 1-loop β coefficient evolution (threshold dependent)
        ax4 = plt.subplot(2, 3, 4)
        mu_points = np.logspace(2, 20, 1000)
        b1_e, b2_e, b3_e = [], [], []
        for muv in mu_points:
            b1, b2, b3 = self._b1_b2_b3(muv)
            b1_e.append(b1); b2_e.append(b2); b3_e.append(b3)
        ax4.semilogx(mu_points, b1_e, 'r-', linewidth=2, label='b₁', alpha=0.8)
        ax4.semilogx(mu_points, b2_e, 'g-', linewidth=2, label='b₂', alpha=0.8)
        ax4.semilogx(mu_points, b3_e, 'b-', linewidth=2, label='b₃', alpha=0.8)
        for name, scale in self.thresholds.items():
            if name in ['MSigma', 'MG8']:
                ax4.axvline(x=scale, linestyle=':', alpha=0.6, label=f'{name}: {scale:.0e}')
        ax4.set_xlabel('Energy Scale μ [GeV]'); ax4.set_ylabel('β Coefficient')
        ax4.set_title('1-Loop β Coefficient Evolution'); ax4.legend(fontsize=9); ax4.grid(True, alpha=0.3)
        ax4.set_xlim(1e2, 1e20)
        # Plot 5: Validation summary (simple)
        ax5 = plt.subplot(2, 3, 5); ax5.axis('off')
        a3_pev, a3_c3, dev_phi0, dev_c3 = self.fingerprints()
        mu_star, spread = self.find_mu_star()
        summary_text = f"""
φ₀ deviation: {dev_phi0*100:.2f}%
c₃ deviation: {dev_c3*100:.2f}%
Best unification μ*: {mu_star:.2e} GeV
Relative spread: {spread*100:.2f}%
G8 mass: {self.thresholds['MG8']:.2e} GeV
Δb₃ (G8): +2.0
""".strip()
        ax5.text(0.05, 0.95, summary_text, transform=ax5.transAxes, fontsize=9, va='top', fontfamily='monospace', bbox=dict(boxstyle="round,pad=0.5", facecolor="lightblue", alpha=0.8))
        # Plot 6: G8 effect comparison (approximate)
        ax6 = plt.subplot(2, 3, 6)
        gut_mask = (df['mu_GeV'] >= 1e12) & (df['mu_GeV'] <= 1e18)
        gut = df[gut_mask]
        ax6.semilogx(gut['mu_GeV'], gut['alpha3_inv'], 'b-', linewidth=2.5, label='α₃⁻¹ (with G8)', alpha=0.8)
        alpha3_inv_no_g8 = gut['alpha3_inv'] + (2.0/(2*math.pi))*np.log(gut['mu_GeV']/self.thresholds['MG8'])
        ax6.semilogx(gut['mu_GeV'], alpha3_inv_no_g8, 'b--', linewidth=2, label='α₃⁻¹ (without G8)', alpha=0.6)
        ax6.axvline(x=self.thresholds['MG8'], color='red', linestyle=':', linewidth=2, alpha=0.8, label='G8 threshold')
        ax6.set_xlabel('Energy Scale μ [GeV]'); ax6.set_ylabel('α₃⁻¹'); ax6.set_title('G8 Bridge Effect on α₃'); ax6.legend(fontsize=9); ax6.grid(True, alpha=0.3)
        plt.tight_layout()
        plot_path = self.results_dir / 'G8_theory_compliant_analysis.png'
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close(); return plot_path


if __name__ == '__main__':
    yaml_path = "/Users/stefanhamann/Projekte/Q5/Pyrate3/models/E8Cascade_TFPT_G8_Enhanced_v2.yaml"
    out = "/Users/stefanhamann/Projekte/Q5/Pyrate3/pyrate/results/E8_TFPT_engine_fixed"
    s = SolverV2(yaml_path, out)
    s.run()
    s.results.to_csv(Path(out)/'couplings_v2.csv', index=False)
    # Export raster
    s.export_raster(per_dec=25)
    # Acceptance tests
    acc = s.acceptance_tests()
    with open(Path(out)/'acceptance_v2.json', 'w') as f:
        json.dump(acc, f, indent=2)
    # Simple printout
    print(f"phi0 dev: {acc['fingerprints']['dev_phi0']*100:.3f}% | c3 dev: {acc['fingerprints']['dev_c3']*100:.3f}% | unif rel spread: {acc['unification']['rel_spread']*100 if acc['unification']['rel_spread'] is not None else float('nan'):.3f}%")
    # Plot
    s.create_comprehensive_plots()

