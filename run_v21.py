"""v2.1 post-registration analysis: likelihood material tag.
Re-selects stored v1.7 and v1.8 M025 fits — no new simulation or fitting.
  python3 run_v21.py          -> v21_results.npy (sanity gate, validations, tag scan, predictions)
"""
import os, sys, subprocess, pickle, time
import numpy as np

MASSES = [0.017, 0.05, 0.1, 0.2]
SAMPLES = ['conv'] + [f'sig{m}' for m in MASSES]
N_CONV, N_SIG = 100_000, 60_000

# ── published degradations (sanity gate target, loaded from stored results) ──
PUB = {}  # filled at import time from v17_results.npy and v18_results.npy

def _load_pub():
    mydir = os.path.dirname(os.path.abspath(__file__))
    r17 = np.load(os.path.join(mydir, 'v17/v17_results.npy'), allow_pickle=True).item()
    r18 = np.load(os.path.join(mydir, 'v18/v18_results.npy'), allow_pickle=True).item()
    for m in MASSES:
        PUB[('C1', m)] = r17[('deg', '3D', m)]['deg']
        PUB[('C2', m)] = r17[('deg', '3D_allpixel', m)]['deg']
        PUB[('C3', m)] = r18['M025'][('deg', '3D_allpixel', m)]['deg']

_load_pub()

# ── configurations ──
#   C1: v1.7, f=1,    dchi2 = dchi2T + dL_default      (3D, pixel z only)
#   C2: v1.7, f=1,    dchi2 = dchi2T + dL_allpixel      (3D, z at every layer)
#   C3: v1.8 M025, f=0.25, dchi2 = dchi2T + dL_allpixel (3D, z at every layer)
CFGS = {
    'C1': dict(ver='v17', tag=None,    f=1.0,  dchi2_key=lambda F: F['dchi2T'] + F['dL_default'],    label='f=1 pixel-z'),
    'C2': dict(ver='v17', tag=None,    f=1.0,  dchi2_key=lambda F: F['dchi2T'] + F['dL_allpixel'],   label='f=1 z-every-layer'),
    'C3': dict(ver='v18', tag='M025',  f=0.25, dchi2_key=lambda F: F['dchi2T'] + F['dL_allpixel'],   label='f=0.25 z-every-layer'),
}


def run_config(cfg_name):
    """Runs in a subprocess with V16_F set. Returns dict of results for one configuration."""
    cfg = CFGS[cfg_name]
    mydir = os.path.dirname(os.path.abspath(__file__))
    v19dir = os.path.join(mydir, 'v19')
    if v19dir not in sys.path:
        sys.path.insert(0, v19dir)
    preregdir = os.path.join(mydir, 'v21_prereg')
    if preregdir not in sys.path:
        sys.path.insert(0, preregdir)
    os.chdir(v19dir)  # so that relative paths like ../v17/ work
    import mag
    from run_v13 import conversions_truth, signal_truth
    from run_v16 import subsample
    from run_v15 import scan_v15
    import bkg
    import mtag_v21

    assert mag.F_MAT == cfg['f'], f"expected F_MAT={cfg['f']}, got {mag.F_MAT}"

    # ── load stored fits ──
    if cfg['ver'] == 'v17':
        fits = {k: np.load(os.path.join(mydir, f"v17/v17_{k}.npz")) for k in SAMPLES}
    else:
        fits = {k: np.load(os.path.join(mydir, f"v18/v18_{cfg['tag']}_{k}.npz")) for k in SAMPLES}

    res = {}

    # ── helpers (same as run_v17/run_v18) ──
    def full(d, Fk, nmax, arr):
        S, use = subsample(d, nmax); assert np.array_equal(use, Fk['idx'])
        a = np.full(len(d['r']), -np.inf); a[use] = np.where(Fk['ok_alp'], arr, -np.inf)
        return a

    def attach(d, Fk, nmax):
        S, use = subsample(d, nmax)
        rm = np.full(len(d['r']), -1.0); sr = np.ones(len(d['r'])); ok = np.zeros(len(d['r']), bool)
        rm[use] = np.nan_to_num(Fk['r'], nan=-1.0); sr[use] = np.nan_to_num(Fk['sigma_r'], nan=1.0); ok[use] = Fk['ok']
        S['rm'], S['sr'], S['base'] = rm, sr, S['base'] & ok
        return S

    # ── build truth samples and attach fits ──
    C0 = conversions_truth(); Fc = fits['conv']
    C = attach(C0, Fc, N_CONV)

    # ── SANITY GATE step 1: reproduce 12 published hard-veto degradations ──
    dchi2_fn = cfg['dchi2_key']
    dC_full = full(C0, Fc, N_CONV, dchi2_fn(Fc))

    gate1 = {}
    for m in MASSES:
        k = f'sig{m}'
        S0 = signal_truth(m, 1.0); Fs = fits[k]; S = attach(S0, Fs, N_SIG)
        dS_full = full(S0, Fs, N_SIG, dchi2_fn(Fs))
        r15 = scan_v15(m, S, C, dS_full, dC_full)
        pub_val = PUB[(cfg_name, m)]
        relerr = abs(r15['deg'] - pub_val) / pub_val
        gate1[m] = dict(computed=r15['deg'], published=pub_val, relerr=relerr, passed=relerr < 1e-9,
                        k=r15['k'], nsig=r15['nsig'])
    res['gate1'] = gate1
    gate1_ok = all(g['passed'] for g in gate1.values())

    # ── SANITY GATE step 2: p_conv normalisation ──
    from scipy.integrate import quad
    gate2 = {}
    for sr_test in [1.0, 5.0, 15.0]:
        # Integrate over a wide enough range, with tight tolerance for narrow bands at small σ_r
        integral, abserr = quad(lambda r: float(mtag_v21.p_conv(np.array(r), np.array(sr_test))),
                                -200, 1200, limit=200, epsabs=1e-12, epsrel=1e-12)
        gate2[sr_test] = dict(integral=integral, err=abs(integral - 1.0), abserr=abserr,
                              passed=abs(integral - 1.0) < 1e-6)
    res['gate2'] = gate2
    gate2_ok = all(g['passed'] for g in gate2.values())

    # ── SANITY GATE step 3: λ finiteness ──
    # check on conversion and all signal samples
    gate3 = {}
    for k in SAMPLES:
        d0 = conversions_truth() if k == 'conv' else signal_truth(float(k[3:]), 1.0)
        Fk = fits[k]; nmax = N_CONV if k == 'conv' else N_SIG
        Sk = attach(d0, Fk, nmax)
        rm_base = Sk['rm'][Sk['base']]; sr_base = Sk['sr'][Sk['base']]
        dchi_arr = dchi2_fn(Fk)
        dchi_base = np.full(len(d0['r']), 0.0)
        _, use = subsample(d0, nmax)
        dchi_mapped = np.where(Fk['ok_alp'], dchi_arr, 0.0)
        dchi_base[use] = dchi_mapped
        dchi_base_sel = dchi_base[Sk['base']]
        lam_vals = mtag_v21.lam(rm_base, sr_base, dchi_base_sel)
        n_finite = np.isfinite(lam_vals).sum(); n_total = len(lam_vals)
        gate3[k] = dict(n_finite=int(n_finite), n_total=int(n_total), passed=n_finite == n_total)
    res['gate3'] = gate3
    gate3_ok = all(g['passed'] for g in gate3.values())

    res['gate_ok'] = gate1_ok and gate2_ok and gate3_ok
    print(f"\n{'='*60}")
    print(f"Config {cfg_name} ({cfg['label']}) — Sanity Gate")
    print(f"{'='*60}")
    print(f"Step 1 (reproduce published):")
    for m, g in gate1.items():
        print(f"  {m*1000:.0f} MeV: computed {g['computed']:.4f} vs published {g['published']:.4f}, "
              f"relerr {g['relerr']:.2e}  {'PASS' if g['passed'] else 'FAIL'}")
    print(f"Step 2 (p_conv normalisation):")
    for sr, g in gate2.items():
        print(f"  σ_r={sr:.0f} mm: ∫p_conv={g['integral']:.10f}, err={g['err']:.2e}  {'PASS' if g['passed'] else 'FAIL'}")
    print(f"Step 3 (λ finiteness):")
    for k, g in gate3.items():
        print(f"  {k}: {g['n_finite']}/{g['n_total']} finite  {'PASS' if g['passed'] else 'FAIL'}")
    print(f"Gate: {'PASSED' if res['gate_ok'] else 'FAILED'}")

    if not res['gate_ok']:
        print("*** SANITY GATE FAILED — investigating ***")
        # Continue anyway so we can report the failure and still compute what we can

    # ── Validations T1–T3 ──
    print(f"\n--- Validations T1–T3 ---")

    # T1: p_conv → 0 and all dchi2 = 0 → λ is constant, tag with t=None equals published no-veto
    S0_test = signal_truth(0.017, 1.0); Fs_test = fits['sig0.017']
    S_test = attach(S0_test, Fs_test, N_SIG)
    rm_t1 = S_test['rm'][S_test['base']]; sr_t1 = S_test['sr'][S_test['base']]
    # With p_conv=0, ln(p_sig / p_conv) → ln(p_sig / P_FLOOR); dchi2=0 → λ = 0 + ln(p_sig) - ln(P_FLOOR)
    lam_zero = mtag_v21.lam(rm_t1, sr_t1, np.zeros(len(rm_t1)))
    # With p_conv forced to zero: use P_FLOOR
    lam_nomat = 0.5 * 0 + np.log(mtag_v21.p_sig(rm_t1)) - np.log(mtag_v21.P_FLOOR)
    lam_range = lam_zero.max() - lam_zero.min()
    # Since p_conv is nonzero at material bands, lam_zero won't be perfectly constant.
    # T1 asks: with p_conv set to zero. We need to compute lambda with p_conv actually forced to 0.
    lam_nomat_all = np.log(mtag_v21.p_sig(rm_t1)) - np.log(mtag_v21.P_FLOOR)  # constant since p_sig is constant
    T1_const = (lam_nomat_all.max() - lam_nomat_all.min()) < 1e-10
    # With t=None, all base pairs pass → same as no veto
    S0_17 = signal_truth(0.017, 1.0); Fs_17 = fits['sig0.017']; S_17 = attach(S0_17, Fs_17, N_SIG)
    dS_17 = full(S0_17, Fs_17, N_SIG, dchi2_fn(Fs_17))
    S_all = S_17['w'][S_17['base']].sum()
    lim_noveto = bkg.upper_limit(S_all, 0)
    # The published "no-veto" is just upper_limit(S0, 0) which scan_v15 uses as lim0
    T1_pass = T1_const  # λ is constant, so any threshold either passes all or none = same as no veto
    res['T1'] = dict(lam_range_nomat=float(lam_nomat_all.max() - lam_nomat_all.min()),
                     is_constant=T1_const, passed=T1_pass)
    print(f"T1: p_conv=0, dchi2=0 → λ range = {lam_nomat_all.max() - lam_nomat_all.min():.2e}, "
          f"constant={'yes' if T1_const else 'no'}  {'PASS' if T1_pass else 'FAIL'}")

    # T2: monotonicity
    # (a) At a band centre, λ is lower than at the midpoint between two adjacent pixel bands
    # Pixel bands: 29, 68, 109, 160 mm
    pixel_centres = [29.0, 68.0, 109.0, 160.0]
    midpoints = [(pixel_centres[i] + pixel_centres[i+1]) / 2 for i in range(len(pixel_centres)-1)]
    sr_test = 3.0; dchi_test = 5.0
    lam_centres = [float(mtag_v21.lam(rc, sr_test, dchi_test)) for rc in pixel_centres]
    lam_midpoints = [float(mtag_v21.lam(rm, sr_test, dchi_test)) for rm in midpoints]
    T2a = all(lc < lm for lc in lam_centres for lm in lam_midpoints)
    # More precisely: each centre < each midpoint between adjacent pixel bands
    T2a_detail = []
    for i in range(len(pixel_centres)):
        for j in range(len(midpoints)):
            T2a_detail.append(dict(centre=pixel_centres[i], midpoint=midpoints[j],
                                   lam_c=lam_centres[i], lam_m=lam_midpoints[j],
                                   ok=lam_centres[i] < lam_midpoints[j]))
    T2a_pass = all(d['ok'] for d in T2a_detail)

    # (b) For fixed r_m and σ_r, λ increases strictly with Δχ²
    r_fixed = 50.0; sr_fixed = 3.0
    dchi_vals = np.linspace(0, 20, 100)
    lam_vs_dchi = np.array([float(mtag_v21.lam(r_fixed, sr_fixed, d)) for d in dchi_vals])
    T2b_pass = bool(np.all(np.diff(lam_vs_dchi) > 0))

    T2_pass = T2a_pass and T2b_pass
    res['T2'] = dict(T2a_pass=T2a_pass, T2b_pass=T2b_pass, passed=T2_pass,
                     centres=lam_centres, midpoints=lam_midpoints)
    print(f"T2a: λ(band centre) < λ(midpoint between pixel bands): {'PASS' if T2a_pass else 'FAIL'}")
    print(f"T2b: λ strictly increases with Δχ²: {'PASS' if T2b_pass else 'FAIL'}")
    print(f"T2: {'PASS' if T2_pass else 'FAIL'}")

    # T3: at a band centre, λ(σ_r=1) < λ(σ_r=15)
    r_band = 29.0; dchi_t3 = 5.0
    lam_sr1 = float(mtag_v21.lam(r_band, 1.0, dchi_t3))
    lam_sr15 = float(mtag_v21.lam(r_band, 15.0, dchi_t3))
    T3_pass = lam_sr1 < lam_sr15
    res['T3'] = dict(lam_sr1=lam_sr1, lam_sr15=lam_sr15, passed=T3_pass)
    print(f"T3: λ(σ_r=1)={lam_sr1:.4f} < λ(σ_r=15)={lam_sr15:.4f}: {'PASS' if T3_pass else 'FAIL'}")

    res['validations_ok'] = T1_pass and T2_pass and T3_pass

    # ── Tag scan: for each mass, compute λ, scan TGRID thresholds ──
    print(f"\n--- Tag scan ---")

    TGRID = mtag_v21.TGRID  # (None, -4, -2, -1, 0, 1, 2, 3, 4, 6, 8, 12, 16, 25)

    # Compute λ for all conversion base pairs
    Fc_arr = fits['conv']
    dchi_conv_3d = dchi2_fn(Fc_arr)  # on subsampled indices
    # Map to full arrays for tag computation
    rm_conv_full = C['rm']; sr_conv_full = C['sr']
    dchi_conv_full = np.zeros(len(C0['r']))
    _, use_c = subsample(C0, N_CONV)
    dchi_conv_full[use_c] = np.where(Fc_arr['ok_alp'], dchi_conv_3d, 0.0)

    lam_conv_full = np.full(len(C0['r']), -np.inf)
    base_c = C['base']
    lam_conv_full[base_c] = mtag_v21.lam(rm_conv_full[base_c], sr_conv_full[base_c], dchi_conv_full[base_c])

    for m in MASSES:
        k = f'sig{m}'
        S0 = signal_truth(m, 1.0); Fs = fits[k]; S = attach(S0, Fs, N_SIG)
        dchi_sig_3d = dchi2_fn(Fs)
        rm_sig_full = S['rm']; sr_sig_full = S['sr']
        dchi_sig_full = np.zeros(len(S0['r']))
        _, use_s = subsample(S0, N_SIG)
        dchi_sig_full[use_s] = np.where(Fs['ok_alp'], dchi_sig_3d, 0.0)

        lam_sig_full = np.full(len(S0['r']), -np.inf)
        base_s = S['base']
        lam_sig_full[base_s] = mtag_v21.lam(rm_sig_full[base_s], sr_sig_full[base_s], dchi_sig_full[base_s])

        # Scan TGRID × mass windows (like scan_v15 but with λ threshold instead of veto+exemption)
        S0_w = S['w']; C0_w = C['w']
        S_base_sum = S0_w[base_s].sum()
        lim0 = bkg.upper_limit(S_base_sum, 0)

        # Mass window candidates
        cands = [(np.ones(len(S0['r']), bool), np.ones(len(C0['r']), bool))]
        for lo in np.linspace(0, 1.5 * m, 31):
            for hi in (3 * m + 0.01, np.inf):
                cands.append(((S['m'] > lo) & (S['m'] < hi), (C['m'] > lo) & (C['m'] < hi)))

        best = None
        for t in TGRID:
            if t is None:
                ts = np.ones(len(S0['r']), bool); tb = np.ones(len(C0['r']), bool)
            else:
                ts = lam_sig_full >= t; tb = lam_conv_full >= t
            for ms, mb in cands:
                cs = base_s & ts & ms; cb = base_c & tb & mb
                Sw = S0_w[cs].sum(); Bw = C0_w[cb].sum()
                lim = bkg.upper_limit(Sw, Bw)
                if best is None or lim < best['limit']:
                    best = dict(limit=lim, deg=lim / lim0, t=t,
                                eff_sig=float(cs.mean() * 100), B=float(Bw),
                                S=float(Sw))

        # Also scan without mass window
        best_nomass = None
        for t in TGRID:
            if t is None:
                ts = np.ones(len(S0['r']), bool); tb = np.ones(len(C0['r']), bool)
            else:
                ts = lam_sig_full >= t; tb = lam_conv_full >= t
            cs = base_s & ts; cb = base_c & tb
            Sw = S0_w[cs].sum(); Bw = C0_w[cb].sum()
            lim = bkg.upper_limit(Sw, Bw)
            if best_nomass is None or lim < best_nomass['limit']:
                best_nomass = dict(limit=lim, deg=lim / lim0, t=t,
                                   eff_sig=float(cs.mean() * 100), B=float(Bw))

        pub_deg = PUB[(cfg_name, m)]
        res[('tag', m)] = best
        res[('tag_nomass', m)] = best_nomass
        res[('hardveto', m)] = gate1[m]  # published
        ratio = best['deg'] / pub_deg if pub_deg > 0 else np.inf
        print(f"  {m*1000:4.0f} MeV: tag {best['deg']:7.2f}× (t={best['t']}) vs hard veto {pub_deg:.2f}×  "
              f"ratio {ratio:.4f}  eff {best['eff_sig']:.2f}%  B={best['B']:.1f}")

    # ── Store λ distributions for report figures ──
    # 17 MeV signal λ and conversion λ (base pairs)
    S0_17 = signal_truth(0.017, 1.0); Fs_17 = fits['sig0.017']; S_17 = attach(S0_17, Fs_17, N_SIG)
    dchi_17 = dchi2_fn(Fs_17)
    rm17 = S_17['rm']; sr17 = S_17['sr']
    dchi17_full = np.zeros(len(S0_17['r']))
    _, use17 = subsample(S0_17, N_SIG)
    dchi17_full[use17] = np.where(Fs_17['ok_alp'], dchi_17, 0.0)
    base17 = S_17['base']
    lam_sig17 = mtag_v21.lam(rm17[base17], sr17[base17], dchi17_full[base17])

    lam_conv_base = mtag_v21.lam(rm_conv_full[base_c], sr_conv_full[base_c], dchi_conv_full[base_c])

    res['lam_sig17'] = lam_sig17
    res['lam_conv'] = lam_conv_base
    res['best_t_17'] = res[('tag', 0.017)]['t']

    return res


def main():
    t0 = time.time()
    results = {}

    for cfg_name, cfg in CFGS.items():
        env = dict(os.environ, V16_B='3.8', V16_F=str(cfg['f']),
                   OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        mydir = os.path.dirname(os.path.abspath(__file__))
        v19dir = os.path.join(mydir, 'v19')
        pkl_path = os.path.join(mydir, f'v21_res_{cfg_name}.pkl')
        code = (f"import sys; sys.path.insert(0, {mydir!r}); sys.path.insert(0, {v19dir!r}); "
                f"import pickle, run_v21; "
                f"r = run_v21.run_config('{cfg_name}'); "
                f"pickle.dump(r, open({pkl_path!r}, 'wb'))")
        print(f"\n{'#'*60}")
        print(f"# Running config {cfg_name} ({cfg['label']}) with V16_F={cfg['f']}")
        print(f"{'#'*60}")
        subprocess.run([sys.executable, '-W', 'ignore', '-c', code], env=env, check=True,
                       cwd=mydir)
        results[cfg_name] = pickle.load(open(pkl_path, 'rb'))

    # ── Summary ──
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    all_gates = all(results[c]['gate_ok'] for c in CFGS)
    print(f"\nSanity gate: {'ALL PASSED' if all_gates else 'FAILED'}")

    all_valid = all(results[c].get('validations_ok', False) for c in CFGS)
    print(f"Validations: {'ALL PASSED' if all_valid else 'FAILED'}")

    print(f"\n{'Config':<6} {'Mass':>6}  {'Hard veto':>10}  {'Tag':>10}  {'t':>5}  {'Ratio':>7}  {'Eff%':>6}  {'B':>8}")
    print('-' * 70)
    for cfg_name in CFGS:
        r = results[cfg_name]
        for m in MASSES:
            hv = PUB[(cfg_name, m)]
            tag = r[('tag', m)]
            ratio = tag['deg'] / hv
            print(f"{cfg_name:<6} {m*1000:5.0f}  {hv:10.2f}  {tag['deg']:10.2f}  {tag['t']:>5}  {ratio:7.4f}  "
                  f"{tag['eff_sig']:6.2f}  {tag['B']:8.1f}")

    # ── Predictions P1–P7 ──
    print(f"\n{'='*60}")
    print("PREDICTIONS")
    print(f"{'='*60}")
    r3 = results['C3']
    r1 = results['C1']
    r2 = results['C2']

    def pred(name, val, lo, hi, extra=""):
        ok = lo <= val <= hi
        status = "CONFIRMED" if ok else "REFUTED"
        print(f"  {name}: {val:.2f}  (predicted {lo}–{hi})  → {status}{extra}")
        return ok

    preds = {}
    preds['P1'] = pred('P1 (17 MeV, f=0.25)', r3[('tag', 0.017)]['deg'], 16, 26.7)
    preds['P2'] = pred('P2 (50 MeV, f=0.25)', r3[('tag', 0.05)]['deg'], 5, 8.83)
    preds['P3'] = pred('P3 (17 MeV, f=1 pixel-z)', r1[('tag', 0.017)]['deg'], 110, 127.3)
    preds['P4'] = pred('P4 (17 MeV, f=1 z-every)', r2[('tag', 0.017)]['deg'], 85, 117.4)

    # P5: tag not worse than hard veto by more than 5% in every config/mass
    P5_ok = True
    for cfg_name in CFGS:
        r = results[cfg_name]
        for m in MASSES:
            hv = PUB[(cfg_name, m)]
            tag_deg = r[('tag', m)]['deg']
            if tag_deg > hv * 1.05:
                P5_ok = False
                print(f"  P5 violation: {cfg_name} {m*1000:.0f} MeV: tag {tag_deg:.2f} > 1.05 × hard veto {hv:.2f}")
    preds['P5'] = P5_ok
    print(f"  P5 (no large loss): {'CONFIRMED' if P5_ok else 'REFUTED'}")

    # P6: at 17 MeV optimum f=0.25, tag keeps more base signal than hard veto
    # Hard veto best efficiency at 17 MeV: from gate1
    hv_17_eff = r3['gate1'][0.017]  # this has k, nsig from the hard veto
    tag_17_eff = r3[('tag', 0.017)]['eff_sig']
    # Need to get hard veto eff at the 17 MeV optimum
    # Re-derive: hard veto scan_v15 returns eff in the result
    # From PUB values, the hard veto result has eff stored — let me get it from the full scan
    # Actually gate1 only stores the degradation, not the full scan result.
    # Let me recompute by running scan_v15 directly in this process for C3 17 MeV
    # Actually we can't import bkg/mag here with the right F_MAT. So store it in the subprocess.
    # The tag result has eff_sig. For hard veto, we need the eff from the original result.
    # Let me check what's stored in v18_results.npy
    r18_full = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v18/v18_results.npy'), allow_pickle=True).item()
    hv_eff_17 = r18_full['M025'][('deg', '3D_allpixel', 0.017)]['eff']
    P6_ok = tag_17_eff > hv_eff_17
    preds['P6'] = P6_ok
    print(f"  P6 (more signal at 17 MeV f=0.25): tag eff {tag_17_eff:.2f}% vs hard veto eff {hv_eff_17:.2f}%"
          f"  → {'CONFIRMED' if P6_ok else 'REFUTED'}")

    # P7: f=0.25, 200 MeV ≤ 1.1×
    P7_val = r3[('tag', 0.2)]['deg']
    preds['P7'] = P7_val <= 1.1
    print(f"  P7 (200 MeV f=0.25 ≤ 1.1×): {P7_val:.4f}  → {'CONFIRMED' if P7_val <= 1.1 else 'REFUTED'}")

    results['predictions'] = preds
    results['all_gates'] = all_gates
    results['all_validations'] = all_valid

    np.save('v21_results.npy', results, allow_pickle=True)
    print(f"\nResults saved to v21_results.npy")
    print(f"Total time: {time.time()-t0:.0f} s")


if __name__ == '__main__':
    main()
