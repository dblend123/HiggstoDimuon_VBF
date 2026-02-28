#!/usr/bin/env python3

import numpy as np
import ROOT


def weighted_mean_std(values, weights):
    values = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if values.size == 0:
        return np.nan, np.nan
    wsum = weights.sum()
    if wsum <= 0:
        return np.nan, np.nan
    mean = float(np.sum(weights * values) / wsum)
    var = float(np.sum(weights * (values - mean) ** 2) / wsum)
    return mean, np.sqrt(max(var, 0.0))


def chebyshev_terms(x, order):
    x = np.asarray(x, dtype=np.float64)
    terms = []
    if order >= 1:
        terms.append(x)
    if order >= 2:
        terms.append(2.0 * x * x - 1.0)
    return terms


def smooth_histogram_bins(values, sigma_bins):
    vals = np.asarray(values, dtype=np.float64)
    sigma = float(abs(sigma_bins))
    if vals.size == 0 or sigma <= 1e-8:
        return vals.copy()
    radius = max(1, int(np.ceil(4.0 * sigma)))
    grid = np.arange(-radius, radius + 1, dtype=np.float64)
    ker = np.exp(-0.5 * (grid / sigma) ** 2)
    ker /= max(np.sum(ker), 1e-12)
    return np.convolve(vals, ker, mode="same")


def morph_signal_shape(shape, shift_bins=0.0, sigma_bins=0.0):
    arr = np.asarray(shape, dtype=np.float64)
    if arr.size == 0:
        return arr.copy()
    x = np.arange(arr.size, dtype=np.float64)
    shifted = np.interp(x - float(shift_bins), x, arr, left=0.0, right=0.0)
    if float(abs(sigma_bins)) > 1e-8:
        shifted = smooth_histogram_bins(shifted, sigma_bins)
    shifted = np.clip(shifted, 1e-12, None)
    shifted /= max(float(np.sum(shifted)), 1e-12)
    return shifted


def double_crystal_ball_pdf(x, mean, sigma, alpha_l, n_l, alpha_r, n_r):
    x = np.asarray(x, dtype=np.float64)
    sigma = max(float(sigma), 1e-6)
    alpha_l = max(float(alpha_l), 1e-3)
    alpha_r = max(float(alpha_r), 1e-3)
    n_l = max(float(n_l), 1.1)
    n_r = max(float(n_r), 1.1)

    t = (x - float(mean)) / sigma
    out = np.zeros_like(t, dtype=np.float64)

    core = (t > -alpha_l) & (t < alpha_r)
    out[core] = np.exp(-0.5 * t[core] ** 2)

    left = t <= -alpha_l
    if np.any(left):
        a = alpha_l
        A = (n_l / a) ** n_l * np.exp(-0.5 * a * a)
        B = n_l / a - a
        out[left] = A * np.power(np.clip(B - t[left], 1e-12, None), -n_l)

    right = t >= alpha_r
    if np.any(right):
        a = alpha_r
        A = (n_r / a) ** n_r * np.exp(-0.5 * a * a)
        B = n_r / a - a
        out[right] = A * np.power(np.clip(B + t[right], 1e-12, None), -n_r)

    return np.clip(out, 1e-12, None)


def build_dscb_hist_shape(mass_edges, params):
    centers = 0.5 * (mass_edges[:-1] + mass_edges[1:])
    vals = double_crystal_ball_pdf(
        centers,
        params["mean"],
        params["sigma"],
        params["alpha_l"],
        params["n_l"],
        params["alpha_r"],
        params["n_r"],
    )
    vals = np.clip(vals, 1e-12, None)
    vals /= max(float(np.sum(vals)), 1e-12)
    return vals


def estimate_dscb_params(sig_mass, sig_w):
    m = np.asarray(sig_mass, dtype=np.float64)
    w = np.asarray(sig_w, dtype=np.float64)
    if m.size == 0 or np.sum(w) <= 0:
        return {
            "mean": 125.0,
            "sigma": 1.8,
            "alpha_l": 1.5,
            "n_l": 4.0,
            "alpha_r": 1.7,
            "n_r": 4.0,
        }
    mean, rms = weighted_mean_std(m, np.clip(w, 0.0, None))
    sigma = max(0.6, min(float(rms), 4.0))
    return {
        "mean": float(np.clip(mean, 123.0, 127.0)),
        "sigma": sigma,
        "alpha_l": 1.5,
        "n_l": 4.0,
        "alpha_r": 1.7,
        "n_r": 4.0,
    }


def fit_polynomial_modifier(b_hist, core_shape, x, order):
    if int(order) <= 0:
        return np.zeros(0, dtype=np.float64)

    b = np.asarray(b_hist, dtype=np.float64)
    core = np.asarray(core_shape, dtype=np.float64)
    terms = chebyshev_terms(x, int(order))
    if b.size == 0 or core.size == 0 or not terms:
        return np.zeros(int(order), dtype=np.float64)

    ysum = float(np.sum(b))
    if ysum <= 0:
        return np.zeros(int(order), dtype=np.float64)

    target = np.divide(b, ysum * np.clip(core, 1e-12, None)) - 1.0
    design = np.vstack([np.asarray(t, dtype=np.float64) for t in terms]).T
    w = np.sqrt(np.clip(b, 1e-12, None))
    aw = design * w[:, None]
    bw = target * w
    try:
        coeffs, _, _, _ = np.linalg.lstsq(aw, bw, rcond=None)
    except np.linalg.LinAlgError:
        coeffs = np.zeros(int(order), dtype=np.float64)

    out = np.zeros(int(order), dtype=np.float64)
    out[:min(out.size, coeffs.size)] = coeffs[:min(out.size, coeffs.size)]
    return np.clip(out, -5.0, 5.0)


def fit_shared_background(b_templates, mass_edges, poly_order, use_lowmed=False, eps_prior_scale=1.0, keep_mask=None):
    b_hist = np.asarray(b_templates, dtype=np.float64)
    if b_hist.ndim != 2:
        raise ValueError("b_templates must be a 2D array [n_bdt, n_mass]")

    source = b_hist
    if bool(use_lowmed):
        if keep_mask is not None:
            km = np.asarray(keep_mask, dtype=bool)
            if km.size == b_hist.shape[0] and np.any(km):
                if int(np.count_nonzero(km)) >= 2:
                    source = b_hist[km]
                elif b_hist.shape[0] > 1:
                    source = b_hist[:-1]
            elif b_hist.shape[0] > 1:
                source = b_hist[:-1]
        elif b_hist.shape[0] > 1:
            source = b_hist[:-1]

    core = np.sum(source, axis=0)
    core = np.divide(core, max(float(np.sum(core)), 1e-12))

    centers = 0.5 * (mass_edges[:-1] + mass_edges[1:])
    xmin, xmax = float(mass_edges[0]), float(mass_edges[-1])
    x = 2.0 * (centers - xmin) / max(xmax - xmin, 1e-12) - 1.0

    coeffs = []
    eps_sigma = []
    for i in range(b_hist.shape[0]):
        coeffs.append(fit_polynomial_modifier(b_hist[i], core, x, poly_order))
        n_i = float(np.sum(b_hist[i]))
        sigma_i = 3.0 / np.sqrt(max(n_i, 1.0))
        sigma_i = float(np.clip(sigma_i * max(float(eps_prior_scale), 1e-3), 0.05, 0.8))
        eps_sigma.append(sigma_i)
    return core, np.asarray(coeffs, dtype=np.float64), np.asarray(eps_sigma, dtype=np.float64)


def build_shape_inputs(sig_mass, sig_score, sig_w,
                       bkg_mass, bkg_score, bkg_w,
                       bdt_edges, mass_min, mass_max, mass_nbins,
                       bkg_genjets=None,
                       dy_mass=None, dy_score=None, dy_w=None, dy_genjets=None,
                       core_use_lowmed=False, core_poly_order=2, eps_prior_scale=1.0,
                       core_sb_threshold=1e-3):
    mass_edges = np.linspace(float(mass_min), float(mass_max), int(mass_nbins) + 1)
    mass_centers = 0.5 * (mass_edges[:-1] + mass_edges[1:])
    x = 2.0 * (mass_centers - float(mass_min)) / max(float(mass_max - mass_min), 1e-12) - 1.0

    n_bdt = len(bdt_edges) - 1
    n_mass = len(mass_edges) - 1
    s_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)
    b_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)
    dy_ge2_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)
    dy_lt2_hist = np.zeros((n_bdt, n_mass), dtype=np.float64)

    for i, (lo, hi) in enumerate(zip(bdt_edges[:-1], bdt_edges[1:])):
        if i == n_bdt - 1:
            sm = (sig_score >= lo) & (sig_score <= hi)
            bm = (bkg_score >= lo) & (bkg_score <= hi)
        else:
            sm = (sig_score >= lo) & (sig_score < hi)
            bm = (bkg_score >= lo) & (bkg_score < hi)

        s_hist[i], _ = np.histogram(sig_mass[sm], bins=mass_edges, weights=sig_w[sm])
        b_hist[i], _ = np.histogram(bkg_mass[bm], bins=mass_edges, weights=bkg_w[bm])

        if dy_mass is not None and dy_score is not None and dy_w is not None and dy_genjets is not None:
            if i == n_bdt - 1:
                dm = (dy_score >= lo) & (dy_score <= hi)
            else:
                dm = (dy_score >= lo) & (dy_score < hi)
            dm_ge2 = dm & (dy_genjets >= 2)
            dm_lt2 = dm & (dy_genjets < 2)
            dy_ge2_hist[i], _ = np.histogram(dy_mass[dm_ge2], bins=mass_edges, weights=dy_w[dm_ge2])
            dy_lt2_hist[i], _ = np.histogram(dy_mass[dm_lt2], bins=mass_edges, weights=dy_w[dm_lt2])
        elif bkg_genjets is not None:
            if i == n_bdt - 1:
                gsel = (bkg_score >= lo) & (bkg_score <= hi)
            else:
                gsel = (bkg_score >= lo) & (bkg_score < hi)
            gsel_ge2 = gsel & (bkg_genjets >= 2)
            gsel_lt2 = gsel & (bkg_genjets < 2)
            dy_ge2_hist[i], _ = np.histogram(bkg_mass[gsel_ge2], bins=mass_edges, weights=bkg_w[gsel_ge2])
            dy_lt2_hist[i], _ = np.histogram(bkg_mass[gsel_lt2], bins=mass_edges, weights=bkg_w[gsel_lt2])

    s_y = s_hist.sum(axis=1)
    b_y = b_hist.sum(axis=1)
    dy_ge2_y = dy_ge2_hist.sum(axis=1)
    dy_lt2_y = dy_lt2_hist.sum(axis=1)
    other_y = np.clip(b_y - dy_ge2_y - dy_lt2_y, 0.0, None)

    dscb_params = estimate_dscb_params(sig_mass, sig_w)
    s_shape_dscb = build_dscb_hist_shape(mass_edges, dscb_params)

    sb_ratio = np.divide(s_y, np.maximum(b_y, 1e-12), out=np.full_like(s_y, np.inf, dtype=np.float64), where=b_y > 0)
    keep_mask = sb_ratio < float(core_sb_threshold)

    core, core_modifiers, eps_sigma = fit_shared_background(
        b_hist,
        mass_edges,
        poly_order=int(np.clip(core_poly_order, 0, 2)),
        use_lowmed=bool(core_use_lowmed),
        eps_prior_scale=float(eps_prior_scale),
        keep_mask=keep_mask,
    )

    asimov = s_hist + b_hist
    return {
        "mass_edges": mass_edges,
        "mass_centers": mass_centers,
        "x": x,
        "s_hist": s_hist,
        "b_hist": b_hist,
        "s_y": s_y,
        "b_y": b_y,
        "dy_ge2_y": dy_ge2_y,
        "dy_lt2_y": dy_lt2_y,
        "other_y": other_y,
        "s_shape_dscb": s_shape_dscb,
        "dscb_params": dscb_params,
        "core": core,
        "eps_anchor": core_modifiers,
        "eps_sigma": eps_sigma,
        "sb_ratio": sb_ratio,
        "core_keep_mask": keep_mask,
        "asimov": asimov,
        "n_bdt": n_bdt,
        "n_mass": n_mass,
    }


def unpack_shape_nuisance(vec, n_bdt, poly_order):
    vec = np.asarray(vec, dtype=np.float64)
    n_poly = max(int(poly_order), 0)
    n_eps = n_bdt * n_poly
    n_eta = n_bdt
    need = n_eps + n_eta + 5
    if vec.size < need:
        out = np.zeros(need, dtype=np.float64)
        out[:vec.size] = vec
        vec = out
    eps = vec[:n_eps].reshape((n_bdt, n_poly)) if n_poly > 0 else np.zeros((n_bdt, 0), dtype=np.float64)
    eta = vec[n_eps:n_eps + n_eta]
    tail = vec[n_eps + n_eta:n_eps + n_eta + 5]
    k_dy_ge2 = float(tail[0])
    k_dy_lt2 = float(tail[1])
    k_sig = float(tail[2])
    k_mass = float(tail[3])
    k_res = float(tail[4])
    return eps, eta, k_dy_ge2, k_dy_lt2, k_sig, k_mass, k_res


def build_bkg_shapes_from_nuisance(core_shape, x, eps_row, poly_order):
    mod = np.ones_like(core_shape, dtype=np.float64)
    terms = chebyshev_terms(x, poly_order)
    for j in range(min(len(terms), eps_row.size)):
        mod += float(eps_row[j]) * terms[j]
    mod = np.clip(mod, 1e-6, None)
    shaped = core_shape * mod
    shaped /= max(float(np.sum(shaped)), 1e-12)
    return shaped


def build_initial_nuisance(shape_inputs, poly_order):
    n_bdt = shape_inputs["n_bdt"]
    n_poly = max(int(poly_order), 0)
    return np.zeros(n_bdt * n_poly + n_bdt + 5, dtype=np.float64)


def nuisance_penalty(nuisance_vec, shape_inputs, poly_order, poly_unc):
    v = np.asarray(nuisance_vec, dtype=np.float64)
    n_bdt = shape_inputs["n_bdt"]
    n_poly = max(int(poly_order), 0)
    n_eps = n_bdt * n_poly
    n_eta = n_bdt
    penalty = 0.0

    if n_poly > 0:
        eps = v[:n_eps].reshape((n_bdt, n_poly))
        eps_sigma_base = np.asarray(shape_inputs.get("eps_sigma", np.ones(n_bdt, dtype=np.float64) * float(poly_unc)), dtype=np.float64)
        if eps_sigma_base.size != n_bdt:
            eps_sigma_base = np.ones(n_bdt, dtype=np.float64) * float(poly_unc)
        eps_sigma = np.repeat(np.clip(eps_sigma_base, 1e-6, None), n_poly).reshape((n_bdt, n_poly))
        pull = eps / eps_sigma
        penalty += 0.5 * float(np.sum(pull ** 2))

    eta = v[n_eps:n_eps + n_eta]
    if eta.size:
        penalty += 0.5 * float(np.sum(eta ** 2))

    tail = v[n_eps + n_eta:n_eps + n_eta + 5]
    if tail.size < 5:
        tmp = np.zeros(5, dtype=np.float64)
        tmp[:tail.size] = tail
        tail = tmp
    penalty += 0.5 * float(tail[0] ** 2)
    penalty += 0.0 * float(tail[1] ** 2)
    penalty += 0.5 * float(tail[2] ** 2)
    penalty += 0.5 * float(tail[3] ** 2)
    penalty += 0.5 * float(tail[4] ** 2)
    return penalty


def build_expected_counts(mu, nuisance_vec, shape_inputs,
                          poly_order=1, bkg_norm_unc=0.20, poly_unc=0.30,
                          sig_norm_unc=0.05, sig_mass_shift_unc=0.20, sig_resolution_unc=0.10,
                          dy_ge2_unc=0.05):
    n_bdt = shape_inputs["n_bdt"]
    eps, eta, k_dy_ge2, k_dy_lt2, k_sig, k_mass, k_res = unpack_shape_nuisance(nuisance_vec, n_bdt, poly_order)

    mass_edges = shape_inputs["mass_edges"]
    mass_binw = float(mass_edges[1] - mass_edges[0]) if mass_edges.size >= 2 else 1.0
    shift_bins = (k_mass * float(sig_mass_shift_unc)) / max(mass_binw, 1e-12)
    smear_bins = abs(k_res * float(sig_resolution_unc))

    lam = np.zeros((shape_inputs["n_bdt"], shape_inputs["n_mass"]), dtype=np.float64)

    for i in range(shape_inputs["n_bdt"]):
        s_shape = morph_signal_shape(shape_inputs["s_shape_dscb"], shift_bins=shift_bins, sigma_bins=smear_bins)
        s_norm = max(float(mu), 0.0) * np.exp(float(sig_norm_unc) * k_sig)
        s = s_norm * s_shape

        b_ge2_nom = max(float(shape_inputs["dy_ge2_y"][i]), 0.0)
        b_lt2_nom = max(float(shape_inputs["dy_lt2_y"][i]), 0.0)
        b_other_nom = max(float(shape_inputs["other_y"][i]), 0.0)
        b_nom_total = max(b_ge2_nom + b_lt2_nom + b_other_nom, 1e-12)

        b_ge2 = b_ge2_nom * np.exp(float(dy_ge2_unc) * k_dy_ge2)
        b_lt2 = b_lt2_nom * np.exp(float(k_dy_lt2))
        b_other = b_other_nom
        b_components_total = max(b_ge2 + b_lt2 + b_other, 1e-12)

        frac_ge2 = b_ge2_nom / b_nom_total
        if frac_ge2 > 0.8:
            bin_unc = 0.05
        elif frac_ge2 > 0.5:
            bin_unc = 0.10
        else:
            bin_unc = float(bkg_norm_unc)
        b_norm = b_components_total * np.exp(float(bin_unc) * float(eta[i] if i < eta.size else 0.0))

        b_shape = build_bkg_shapes_from_nuisance(
            shape_inputs["core"],
            shape_inputs["x"],
            eps[i] if i < eps.shape[0] else np.zeros(0, dtype=np.float64),
            poly_order,
        )
        b = b_norm * b_shape
        lam[i] = s + b
    return np.clip(lam, 1e-12, None)


def poisson_nll(data, model):
    d = np.asarray(data, dtype=np.float64)
    m = np.asarray(model, dtype=np.float64)
    return float(np.sum(m - d * np.log(np.clip(m, 1e-12, None))))


def total_shape_nll(mu, nuisance_vec, shape_inputs,
                    poly_order=1, bkg_norm_unc=0.20, poly_unc=0.30,
                    sig_norm_unc=0.05, sig_mass_shift_unc=0.20, sig_resolution_unc=0.10,
                    dy_ge2_unc=0.05,
                    data=None):
    lam = build_expected_counts(
        mu,
        nuisance_vec,
        shape_inputs,
        poly_order=poly_order,
        bkg_norm_unc=bkg_norm_unc,
        poly_unc=poly_unc,
        sig_norm_unc=sig_norm_unc,
        sig_mass_shift_unc=sig_mass_shift_unc,
        sig_resolution_unc=sig_resolution_unc,
        dy_ge2_unc=dy_ge2_unc,
    )
    obs = shape_inputs["asimov"] if data is None else np.asarray(data, dtype=np.float64)
    return poisson_nll(obs, lam) + nuisance_penalty(nuisance_vec, shape_inputs, poly_order, poly_unc)


def profile_nll_fixed_mu(mu, shape_inputs, poly_order, bkg_norm_unc, poly_unc,
                         sig_norm_unc, sig_mass_shift_unc, sig_resolution_unc,
                         dy_ge2_unc,
                         data,
                         max_rounds=7, init=None):
    n_bdt = shape_inputs["n_bdt"]
    n_poly = max(int(poly_order), 0)
    n_eps = n_bdt * n_poly
    n_eta = n_bdt
    n_params = n_eps + n_eta + 5
    if init is None:
        vec = build_initial_nuisance(shape_inputs, poly_order)
    else:
        vec = np.asarray(init, dtype=np.float64).copy()
        if vec.size != n_params:
            tmp = np.zeros(n_params, dtype=np.float64)
            tmp[:min(n_params, vec.size)] = vec[:min(n_params, vec.size)]
            vec = tmp

    step_schedule = [1.0, 0.6, 0.35, 0.20, 0.12, 0.07, 0.04]
    step_schedule = step_schedule[:max(1, min(max_rounds, len(step_schedule)))]

    best = total_shape_nll(
        mu, vec, shape_inputs,
        poly_order, bkg_norm_unc, poly_unc,
        sig_norm_unc, sig_mass_shift_unc, sig_resolution_unc,
        dy_ge2_unc,
        data,
    )

    bounds = [(-5.0, 5.0)] * n_params
    if n_poly > 0:
        for i in range(n_eps):
            bounds[i] = (-3.0, 3.0)
    idx_k_dy_lt2 = n_eps + n_eta + 1
    if 0 <= idx_k_dy_lt2 < n_params:
        bounds[idx_k_dy_lt2] = (-3.0, 3.0)

    for step in step_schedule:
        improved = False
        for ip in range(n_params):
            cur = vec[ip]
            cands = np.array([cur - step, cur, cur + step], dtype=np.float64)
            cands = np.clip(cands, bounds[ip][0], bounds[ip][1])
            vals = []
            for c in cands:
                test = vec.copy()
                test[ip] = c
                vals.append(total_shape_nll(
                    mu,
                    test,
                    shape_inputs,
                    poly_order,
                    bkg_norm_unc,
                    poly_unc,
                    sig_norm_unc,
                    sig_mass_shift_unc,
                    sig_resolution_unc,
                    dy_ge2_unc,
                    data,
                ))
            best_idx = int(np.argmin(vals))
            if vals[best_idx] + 1e-10 < best:
                vec[ip] = float(cands[best_idx])
                best = float(vals[best_idx])
                improved = True
        if not improved:
            continue

    return best, vec


def compute_shape_expected_significance(shape_inputs,
                                        poly_order=1,
                                        bkg_norm_unc=0.20,
                                        poly_unc=0.30,
                                        sig_norm_unc=0.05,
                                        sig_mass_shift_unc=0.20,
                                        sig_resolution_unc=0.10,
                                        dy_ge2_unc=0.05,
                                        mu_max=3.0,
                                        mu_points=121,
                                        data_override=None):
    mu_points = max(int(mu_points), 11)
    mu_grid = np.linspace(0.0, float(mu_max), mu_points)
    data = shape_inputs["asimov"] if data_override is None else np.asarray(data_override, dtype=np.float64)

    nll_mu0, vec_mu0 = profile_nll_fixed_mu(
        0.0,
        shape_inputs,
        poly_order,
        bkg_norm_unc,
        poly_unc,
        sig_norm_unc,
        sig_mass_shift_unc,
        sig_resolution_unc,
        dy_ge2_unc,
        data,
        init=None,
    )

    prof = []
    vec_init = vec_mu0
    best = None
    for mu in mu_grid:
        nll_mu, vec_init = profile_nll_fixed_mu(
            float(mu),
            shape_inputs,
            poly_order,
            bkg_norm_unc,
            poly_unc,
            sig_norm_unc,
            sig_mass_shift_unc,
            sig_resolution_unc,
            dy_ge2_unc,
            data,
            init=vec_init,
        )
        prof.append((float(mu), float(nll_mu)))
        if best is None or nll_mu < best[1]:
            best = (float(mu), float(nll_mu), vec_init.copy())

    mu_hat, nll_best, vec_best = best
    q0 = max(0.0, 2.0 * (nll_mu0 - nll_best))
    zexp = float(np.sqrt(q0))

    model_b = build_expected_counts(
        0.0,
        vec_mu0,
        shape_inputs,
        poly_order,
        bkg_norm_unc,
        poly_unc,
        sig_norm_unc,
        sig_mass_shift_unc,
        sig_resolution_unc,
        dy_ge2_unc,
    )
    model_sb = build_expected_counts(
        mu_hat,
        vec_best,
        shape_inputs,
        poly_order,
        bkg_norm_unc,
        poly_unc,
        sig_norm_unc,
        sig_mass_shift_unc,
        sig_resolution_unc,
        dy_ge2_unc,
    )

    return {
        "zexp": zexp,
        "q0": q0,
        "mu_hat": float(mu_hat),
        "nll_mu0": float(nll_mu0),
        "nll_best": float(nll_best),
        "vec_mu0": vec_mu0,
        "vec_best": vec_best,
        "profile": prof,
        "model_b": model_b,
        "model_sb": model_sb,
    }


def run_shape_toy_ensembles(shape_inputs,
                            n_toys,
                            injected_mus,
                            poly_order=1,
                            bkg_norm_unc=0.20,
                            poly_unc=0.30,
                            sig_norm_unc=0.05,
                            sig_mass_shift_unc=0.20,
                            sig_resolution_unc=0.10,
                            dy_ge2_unc=0.05,
                            mu_max=3.0,
                            mu_points=121,
                            seed=12345):
    n_toys = int(max(0, n_toys))
    if n_toys <= 0:
        return {}
    rng = np.random.default_rng(int(seed))
    results = {}
    nom_vec = build_initial_nuisance(shape_inputs, poly_order)

    for mu_inj in injected_mus:
        mu_inj = float(mu_inj)
        lam_true = build_expected_counts(
            mu_inj,
            nom_vec,
            shape_inputs,
            poly_order=poly_order,
            bkg_norm_unc=bkg_norm_unc,
            poly_unc=poly_unc,
            sig_norm_unc=sig_norm_unc,
            sig_mass_shift_unc=sig_mass_shift_unc,
            sig_resolution_unc=sig_resolution_unc,
            dy_ge2_unc=dy_ge2_unc,
        )

        mu_hats = []
        z_vals = []
        for _ in range(n_toys):
            toy = rng.poisson(np.clip(lam_true, 1e-12, None)).astype(np.float64)
            fit = compute_shape_expected_significance(
                shape_inputs,
                poly_order=poly_order,
                bkg_norm_unc=bkg_norm_unc,
                poly_unc=poly_unc,
                sig_norm_unc=sig_norm_unc,
                sig_mass_shift_unc=sig_mass_shift_unc,
                sig_resolution_unc=sig_resolution_unc,
                dy_ge2_unc=dy_ge2_unc,
                mu_max=mu_max,
                mu_points=mu_points,
                data_override=toy,
            )
            mu_hats.append(float(fit["mu_hat"]))
            z_vals.append(float(fit["zexp"]))

        mu_arr = np.asarray(mu_hats, dtype=np.float64)
        z_arr = np.asarray(z_vals, dtype=np.float64)
        results[mu_inj] = {
            "mu_hat_mean": float(np.mean(mu_arr)) if mu_arr.size else np.nan,
            "mu_hat_std": float(np.std(mu_arr)) if mu_arr.size else np.nan,
            "mu_hat_bias": float(np.mean(mu_arr) - mu_inj) if mu_arr.size else np.nan,
            "z_mean": float(np.mean(z_arr)) if z_arr.size else np.nan,
            "z_std": float(np.std(z_arr)) if z_arr.size else np.nan,
            "n_toys": int(n_toys),
        }
    return results


def _safe_sigmoid(x):
    x = np.clip(float(x), -40.0, 40.0)
    return 1.0 / (1.0 + np.exp(-x))


def _bernstein2(t, b1, b2):
    t = np.clip(np.asarray(t, dtype=np.float64), 0.0, 1.0)
    c0 = (1.0 - t) ** 2
    c1 = 2.0 * t * (1.0 - t)
    c2 = t ** 2
    w0 = 1.0
    w1 = np.exp(float(b1))
    w2 = np.exp(float(b2))
    out = w0 * c0 + w1 * c1 + w2 * c2
    return np.clip(out, 1e-12, None)


def _bg_model_shape(name, mass_centers, params, mass_min, mass_max):
    m = np.asarray(mass_centers, dtype=np.float64)
    t = (m - float(mass_min)) / max(float(mass_max - mass_min), 1e-12)
    x = 2.0 * t - 1.0
    mz = 91.1876

    p = np.asarray(params, dtype=np.float64)
    if name == "BWZRedux":
        p1, p2, g = p
        gamma = max(0.5, np.exp(g))
        bw = 1.0 / (((m - mz) ** 2) + gamma ** 2)
        y = np.exp(p1 * x + p2 * x * x) * bw
    elif name == "SumExp":
        p1, p2, a = p
        w = _safe_sigmoid(a)
        y = w * np.exp(p1 * x) + (1.0 - w) * np.exp(p2 * x)
    elif name == "PowerLaw":
        p0, p1, p2 = p
        base = np.power(np.clip(m / max(float(mass_min), 1e-6), 1e-6, None), p0)
        y = base * np.exp(p1 * x + p2 * x * x)
    elif name == "FEWZxBernstein":
        p0, b1, b2 = p
        fewz_like = np.power(np.clip(m / max(float(mz), 1e-6), 1e-6, None), p0)
        y = fewz_like * _bernstein2(t, b1, b2)
    elif name == "LandaUxBernstein":
        mpv, scale, b1 = p
        landau = np.array([ROOT.TMath.Landau(float(mi), float(mpv), max(0.5, np.exp(scale)), False) for mi in m], dtype=np.float64)
        y = landau * np.exp(b1 * x)
    elif name == "BWZxBernstein":
        g, b1, b2 = p
        gamma = max(0.5, np.exp(g))
        bw = 1.0 / (((m - mz) ** 2) + gamma ** 2)
        y = bw * _bernstein2(t, b1, b2)
    elif name == "BWZGamma":
        g, p1, p2 = p
        gamma = max(0.5, np.exp(g))
        bw = 1.0 / (((m - mz) ** 2) + gamma ** 2)
        y = bw * np.power(np.clip(m / max(float(mz), 1e-6), 1e-6, None), p1) * np.exp(p2 * x)
    else:
        raise ValueError(f"Unknown background model: {name}")

    y = np.clip(y, 1e-12, None)
    y /= max(float(np.sum(y)), 1e-12)
    return y


def _bg_model_specs(mass_min, mass_max):
    center = 0.5 * (float(mass_min) + float(mass_max))
    return {
        "BWZRedux": {
            "init": np.array([-1.0, -0.3, np.log(2.0)], dtype=np.float64),
            "bounds": [(-8.0, 2.0), (-5.0, 5.0), (np.log(0.4), np.log(20.0))],
        },
        "SumExp": {
            "init": np.array([-1.2, -0.2, 0.0], dtype=np.float64),
            "bounds": [(-8.0, 2.0), (-8.0, 2.0), (-8.0, 8.0)],
        },
        "PowerLaw": {
            "init": np.array([-1.5, 0.0, 0.0], dtype=np.float64),
            "bounds": [(-10.0, 4.0), (-5.0, 5.0), (-5.0, 5.0)],
        },
        "FEWZxBernstein": {
            "init": np.array([-1.0, 0.0, 0.0], dtype=np.float64),
            "bounds": [(-10.0, 4.0), (-5.0, 5.0), (-5.0, 5.0)],
        },
        "LandaUxBernstein": {
            "init": np.array([center, np.log(2.0), 0.0], dtype=np.float64),
            "bounds": [(float(mass_min), float(mass_max)), (np.log(0.4), np.log(10.0)), (-6.0, 6.0)],
        },
        "BWZxBernstein": {
            "init": np.array([np.log(2.0), 0.0, 0.0], dtype=np.float64),
            "bounds": [(np.log(0.4), np.log(20.0)), (-5.0, 5.0), (-5.0, 5.0)],
        },
        "BWZGamma": {
            "init": np.array([np.log(2.0), -1.0, -0.2], dtype=np.float64),
            "bounds": [(np.log(0.4), np.log(20.0)), (-8.0, 5.0), (-8.0, 5.0)],
        },
    }


def _fit_mu_with_bg_model(obs_hist, s_shape, s_y, b_y, mass_centers, mass_min, mass_max,
                          model_name, mu_start, param_grid=21):
    specs = _bg_model_specs(mass_min, mass_max)
    cfg = specs[model_name]
    pars = np.asarray(cfg["init"], dtype=np.float64).copy()
    bounds = cfg["bounds"]

    mu_hat = max(0.0, float(mu_start))
    ln_b = 0.0

    def nll_val(mu, ln_bkg, pars_vec):
        b_shape = _bg_model_shape(model_name, mass_centers, pars_vec, mass_min, mass_max)
        lam = max(float(mu), 0.0) * float(s_y) * s_shape + float(b_y) * np.exp(float(ln_bkg)) * b_shape
        lam = np.clip(lam, 1e-12, None)
        d = np.asarray(obs_hist, dtype=np.float64)
        return float(np.sum(lam - d * np.log(lam)))

    best = nll_val(mu_hat, ln_b, pars)
    step_schedule = [1.0, 0.5, 0.25, 0.1, 0.05]
    for step in step_schedule:
        improved = True
        while improved:
            improved = False

            for dmu in [-step, 0.0, step]:
                cand = max(0.0, mu_hat + dmu)
                val = nll_val(cand, ln_b, pars)
                if val + 1e-10 < best:
                    mu_hat = cand
                    best = val
                    improved = True

            for db in [-step, 0.0, step]:
                cand = np.clip(ln_b + db, -4.0, 4.0)
                val = nll_val(mu_hat, cand, pars)
                if val + 1e-10 < best:
                    ln_b = cand
                    best = val
                    improved = True

            for ip in range(pars.size):
                lo, hi = bounds[ip]
                delta = (hi - lo) * 0.05 * step
                cands = [pars[ip] - delta, pars[ip], pars[ip] + delta]
                for c in cands:
                    test = pars.copy()
                    test[ip] = np.clip(c, lo, hi)
                    val = nll_val(mu_hat, ln_b, test)
                    if val + 1e-10 < best:
                        pars = test
                        best = val
                        improved = True

    return float(mu_hat), float(ln_b), pars, float(best)


def run_bias_test_matrix(shape_inputs, mu_true=1.0, n_toys=30, threshold=0.20, seed=12345, param_grid=21):
    n_toys = max(1, int(n_toys))
    mu_true = float(mu_true)
    threshold = float(abs(threshold))

    mass_edges = shape_inputs["mass_edges"]
    mass_centers = shape_inputs["mass_centers"]
    mass_min = float(mass_edges[0])
    mass_max = float(mass_edges[-1])

    s_hist_tot = np.sum(shape_inputs["s_hist"], axis=0)
    b_hist_tot = np.sum(shape_inputs["b_hist"], axis=0)
    s_y = float(np.sum(s_hist_tot))
    b_y = float(np.sum(b_hist_tot))
    s_shape = np.clip(np.asarray(s_hist_tot, dtype=np.float64), 1e-12, None)
    s_shape /= max(float(np.sum(s_shape)), 1e-12)

    model_names = list(_bg_model_specs(mass_min, mass_max).keys())
    rng = np.random.default_rng(int(seed))

    matrix = {g: {} for g in model_names}

    for gen_name in model_names:
        gen_shape = _bg_model_shape(gen_name, mass_centers, _bg_model_specs(mass_min, mass_max)[gen_name]["init"], mass_min, mass_max)
        lam_truth = np.clip(mu_true * s_y * s_shape + b_y * gen_shape, 1e-12, None)

        for fit_name in model_names:
            mu_hats = []
            for _ in range(n_toys):
                obs = rng.poisson(lam_truth).astype(np.float64)
                mu_hat, _, _, _ = _fit_mu_with_bg_model(
                    obs,
                    s_shape,
                    s_y,
                    b_y,
                    mass_centers,
                    mass_min,
                    mass_max,
                    fit_name,
                    mu_start=max(mu_true, 0.0),
                    param_grid=param_grid,
                )
                mu_hats.append(mu_hat)

            mu_arr = np.asarray(mu_hats, dtype=np.float64)
            mean_mu = float(np.mean(mu_arr)) if mu_arr.size else np.nan
            std_mu = float(np.std(mu_arr)) if mu_arr.size else np.nan
            if abs(mu_true) > 1e-12:
                bias = float((mean_mu - mu_true) / mu_true)
            else:
                bias = float(mean_mu - mu_true)
            passed = bool(abs(bias) < threshold)
            matrix[gen_name][fit_name] = {
                "mu_hat_mean": mean_mu,
                "mu_hat_std": std_mu,
                "bias": bias,
                "pass": passed,
                "n_toys": int(n_toys),
            }

    return {
        "models": model_names,
        "matrix": matrix,
        "mu_true": mu_true,
        "threshold": threshold,
        "n_toys": int(n_toys),
    }
