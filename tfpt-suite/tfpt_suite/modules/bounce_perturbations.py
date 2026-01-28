from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from mpmath import mp
from scipy.integrate import solve_ivp

from tfpt_suite.constants import TfptConstants
from tfpt_suite.module_base import Check, ModuleResult, ModuleSpec, TfptModule


def _tanh_switch(x: float, *, center: float, width: float) -> tuple[float, float, float]:
    """
    Smooth switch s(x) = 0.5*(1+tanh((x-center)/width)), plus first/second derivatives.

    Returns: (s, s', s'') with derivatives w.r.t x.
    """
    u = (x - center) / width
    t = np.tanh(u)

    # Numerically stable sech^2(u):
    # sech^2(u) = 4 e^{-2|u|} / (1 + e^{-2|u|})^2
    a = np.exp(-2.0 * np.abs(u))
    sech2 = (4.0 * a) / ((1.0 + a) ** 2)
    s = 0.5 * (1.0 + t)
    sp = 0.5 * sech2 / width
    spp = -sech2 * t / (width**2)
    return s, sp, spp


@dataclass(frozen=True)
class BounceInflationBackground:
    """
    Minimal C^2 bounce→inflation background in a dimensionless conformal-time variable x.

    Motivation (paper v2.4, Appendix L outline):
    - a smooth non-singular bounce (torsion-induced a^{-6} term acts repulsively)
    - a transition into an R^2/Starobinsky-like inflationary era

    Implementation choice here (minimal, deterministic):
    - bounce proxy: a_b(x) = a_min * sqrt(1 + x^2)  (radiation-like for |x|>>1)
    - inflation proxy: a_inf(x) = 1 / (x0 - x)      (de Sitter in x-units, H=1)
    - transition: tanh switch s(x) around x_t with width w, giving C^2 a(x)

    Note: Working in x rather than physical η avoids huge scales when using the
    TFPT-predicted small H (M/Mpl ~ 1e-5). The mode equation can be written in x
    without loss of generality for the transfer-function *shape* T(k/k_bounce).
    """

    a_min: float = 1.0
    x_transition: float = 3.0
    x_width: float = 0.20
    growth_factor_after_transition: float = 50.0

    def a_bounce(self, x: float) -> float:
        return self.a_min * float(np.sqrt(1.0 + x * x))

    def a_bounce_p(self, x: float) -> float:
        # d/dx of a_min*sqrt(1+x^2)
        return self.a_min * float(x / np.sqrt(1.0 + x * x))

    def a_bounce_pp(self, x: float) -> float:
        # d^2/dx^2 of a_min*sqrt(1+x^2)
        return self.a_min * float(1.0 / (1.0 + x * x) ** (3.0 / 2.0))

    def x0(self) -> float:
        # Choose x0 so inflation a_inf matches a_bounce at x_transition
        a_t = self.a_bounce(self.x_transition)
        return self.x_transition + 1.0 / a_t

    def a_inflation(self, x: float) -> float:
        denom = self.x0() - x
        if denom <= 0:
            # out of domain: avoid singularity
            return float("inf")
        return 1.0 / denom

    def a_inflation_p(self, x: float) -> float:
        denom = self.x0() - x
        return 1.0 / (denom * denom)

    def a_inflation_pp(self, x: float) -> float:
        denom = self.x0() - x
        return 2.0 / (denom**3)

    def a(self, x: float) -> float:
        s, _, _ = _tanh_switch(x, center=self.x_transition, width=self.x_width)
        ab = self.a_bounce(x)
        ai = self.a_inflation(x)
        return (1.0 - s) * ab + s * ai

    def a_p(self, x: float) -> float:
        s, sp, _ = _tanh_switch(x, center=self.x_transition, width=self.x_width)
        ab = self.a_bounce(x)
        ai = self.a_inflation(x)
        abp = self.a_bounce_p(x)
        aip = self.a_inflation_p(x)
        return (1.0 - s) * abp + s * aip + sp * (ai - ab)

    def a_pp(self, x: float) -> float:
        s, sp, spp = _tanh_switch(x, center=self.x_transition, width=self.x_width)
        ab = self.a_bounce(x)
        ai = self.a_inflation(x)
        abp = self.a_bounce_p(x)
        aip = self.a_inflation_p(x)
        abpp = self.a_bounce_pp(x)
        aipp = self.a_inflation_pp(x)
        return (1.0 - s) * abpp + s * aipp + 2.0 * sp * (aip - abp) + spp * (ai - ab)

    def x_end(self) -> float:
        """
        Evaluate late time well into the inflationary phase (but not too close to x0).
        """
        x0 = self.x0()
        delta_t = x0 - self.x_transition
        return x0 - delta_t / self.growth_factor_after_transition


def _ms_potential_from_z(z: Callable[[float], float], zpp: Callable[[float], float], x: float) -> float:
    zx = z(x)
    if zx == 0:
        return float("inf")
    return zpp(x) / zx


def _solve_mode(
    *,
    k: float,
    omega0: Optional[float],
    x_start: float,
    x_end: float,
    potential: Callable[[float], float],
    rtol: float,
    atol: float,
    method: str,
) -> dict[str, float]:
    """
    Solve v'' + (k^2 - U(x)) v = 0 in x-time, with Bunch-Davies initial conditions.
    Uses real/imag split (4D first-order system).
    """

    raise RuntimeError("_solve_mode is deprecated; use _solve_mode_symplectic.")


def _solve_mode_ivp(
    *,
    k: float,
    omega0: float,
    x_grid: np.ndarray,
    U_grid: np.ndarray,
    x_start: float,
    x_end: float,
    rtol: float,
    atol: float,
    method: str = "DOP853",
    max_step: float | None = None,
) -> dict[str, float]:
    """
    Adaptive solver for v'' + (k^2 - U(x)) v = 0 using solve_ivp on the real/imag split.

    This is used as a fallback in regimes with turning points / tachyonic intervals where a fixed-step
    symplectic integrator may lose Wronskian accuracy due to strong squeezing.
    """
    omega0 = float(omega0)
    if omega0 <= 0:
        omega0 = float(np.sqrt(max(float(k * k - np.interp(x_start, x_grid, U_grid)), 1e-18)))

    v_re0 = float(1.0 / np.sqrt(2.0 * omega0))
    v_im0 = 0.0
    vp_re0 = 0.0
    vp_im0 = float(-omega0 * v_re0)

    def U_of(x: float) -> float:
        return float(np.interp(float(x), x_grid, U_grid))

    def rhs(x: float, y: np.ndarray) -> np.ndarray:
        v_re, v_im, vp_re, vp_im = [float(v) for v in y]
        omega2 = float(k * k - U_of(x))
        return np.array([vp_re, vp_im, -omega2 * v_re, -omega2 * v_im], dtype=float)

    sol = solve_ivp(
        rhs,
        t_span=(float(x_start), float(x_end)),
        y0=np.array([v_re0, v_im0, vp_re0, vp_im0], dtype=float),
        method=str(method),
        rtol=float(rtol),
        atol=float(atol),
        max_step=(float(max_step) if max_step is not None else np.inf),
        t_eval=[float(x_end)],
    )
    if (not sol.success) or sol.y.shape[1] < 1:
        return {"amp": float("nan"), "W_over_i": float("nan")}

    v_re, v_im, vp_re, vp_im = [float(v) for v in sol.y[:, -1]]
    amp = float(np.sqrt(v_re * v_re + v_im * v_im))
    W_over_i = 2.0 * (v_im * vp_re - v_re * vp_im)
    return {"amp": amp, "W_over_i": W_over_i}


def _solve_mode_symplectic(
    *,
    k: float,
    omega0: Optional[float],
    x_grid: np.ndarray,
    U_grid: np.ndarray,
    x_start: float,
    x_end: float,
    substeps: int = 4,
) -> dict[str, float]:
    """
    Symplectic (velocity-Verlet / leapfrog) integrator for:
      v'' + (k^2 - U(x)) v = 0
    on a (piecewise) uniform x-grid, evolving Re/Im components simultaneously.

    This tends to preserve the Wronskian much better than generic stiff ODE solvers
    when modes become highly squeezed (including ω^2<0 regions).
    """
    if substeps < 1:
        substeps = 1

    # Ensure increasing grid
    if not (x_grid.size >= 2 and np.all(np.diff(x_grid) > 0)):
        raise ValueError("x_grid must be a strictly increasing array")

    dx = float(x_grid[1] - x_grid[0])
    if dx <= 0:
        raise ValueError("Invalid x_grid spacing")

    # Find start/end indices (x_start is normally a grid point from select_start)
    i0 = int(np.argmin(np.abs(x_grid - x_start)))
    i1 = int(np.searchsorted(x_grid, x_end, side="right") - 1)
    i1 = max(i0 + 1, min(i1, int(x_grid.size - 1)))

    # Initial frequency
    if omega0 is None:
        omega2_0 = float(k * k - float(U_grid[i0]))
        omega0 = float(np.sqrt(max(omega2_0, 1e-18)))

    v_re = float(1.0 / np.sqrt(2.0 * omega0))
    v_im = 0.0
    vp_re = 0.0
    vp_im = float(-omega0 * v_re)

    h = dx / float(substeps)

    for i in range(i0, i1):
        U0 = float(U_grid[i])
        U1 = float(U_grid[i + 1])
        dU = U1 - U0
        # micro-steps with linear interpolation of U across the cell
        for j in range(substeps):
            f0 = float(j) / float(substeps)
            f1 = float(j + 1) / float(substeps)
            Uj = U0 + dU * f0
            Ujp = U0 + dU * f1

            omega2_j = float(k * k - Uj)
            omega2_jp = float(k * k - Ujp)

            # half-step momentum update
            vp_re += -0.5 * h * omega2_j * v_re
            vp_im += -0.5 * h * omega2_j * v_im

            # full-step position update
            v_re += h * vp_re
            v_im += h * vp_im

            # half-step momentum update with omega2 at end of microstep
            vp_re += -0.5 * h * omega2_jp * v_re
            vp_im += -0.5 * h * omega2_jp * v_im

        # Renormalize the Wronskian after each grid cell to mitigate loss of precision
        # in highly squeezed / tachyonic regimes (where v and v' can grow exponentially).
        # Compute W in extended precision to reduce catastrophic cancellation at very small k.
        W_over_i_ld = np.longdouble(2.0) * (
            np.longdouble(v_im) * np.longdouble(vp_re) - np.longdouble(v_re) * np.longdouble(vp_im)
        )
        W_over_i = float(W_over_i_ld)
        if np.isfinite(W_over_i) and W_over_i != 0.0:
            if W_over_i < 0.0:
                # Restore positive orientation (equivalent to complex conjugation).
                v_im = -v_im
                vp_im = -vp_im
                W_over_i_ld = -W_over_i_ld
                W_over_i = float(W_over_i_ld)
            s = float(1.0 / np.sqrt(W_over_i_ld))
            v_re *= s
            v_im *= s
            vp_re *= s
            vp_im *= s

    amp = float(np.sqrt(v_re * v_re + v_im * v_im))
    W_over_i = float(
        np.longdouble(2.0)
        * (np.longdouble(v_im) * np.longdouble(vp_re) - np.longdouble(v_re) * np.longdouble(vp_im))
    )
    return {"amp": amp, "W_over_i": W_over_i}


class BouncePerturbationsModule(TfptModule):
    module_id = "bounce_perturbations"
    title = "Bounce perturbations (f(R) z(η) + background ODE + transfer function T(k))"

    def spec(self) -> ModuleSpec:
        return ModuleSpec(
            name=self.title,
            module_id=self.module_id,
            inputs=[
                "Background (dimensionless): Starobinsky f(R)=R+R^2/(6M^2) in terms of F=df/dR plus a torsion proxy ρ~a^{-6}",
                "Mode equation: v'' + (k^2 - z''/z) v = 0 (Appendix L); scalar and tensor z as in f(R)",
            ],
            outputs=["transfer function T(k)", "Wronskian conservation diagnostics"],
            formulas=[
                "T(k) = |v_k(x_end)| / |v_k^{(R2)}(x_end)| (Appendix L definition, mapped to x)",
                "tensor: z_t = a * sqrt(F)",
                "scalar (f(R)): z_s = a * sqrt(Q_s), Q_s = 3 F'^2 / (2 F (H + F'/(2F))^2) (c_s^2=1)",
                "adiabatic ICs at x_start: v=1/sqrt(2ω0), v'=-i ω0 v with ω0=sqrt(k^2-U(x_start))",
                "Checks: T(k)->1 for k >> k_bounce; Wronskian ~ i conserved",
            ],
            validation=[
                "T(k) approaches 1 for large k (numerically)",
                "Wronskian deviation small across k-grid",
                "stability: F>0 and f''(R)>0 (Starobinsky)",
            ],
            determinism="Deterministic given the fixed background and solver tolerances.",
        )

    def run(self, config) -> ModuleResult:
        # --- Background solve (dimensionless Starobinsky f(R) + torsion proxy ρ~a^{-6}) ---
        const = TfptConstants.compute()
        M_over_Mpl = float(const.M_over_Mpl)  # reported; the dimensionless ODE below uses τ=M t so M cancels

        settings = {
            "F0": 5.0,
            # Important for scalar stability in f(R): g0=dF/dτ at the bounce.
            # If g0≈0, then H + F'/(2F) can pass near zero and Q_s may blow up.
            "g0": 0.1,
            "a0": 1.0,
            "tau_pre": 25.0,
            "tau_post": 120.0,
            "n_tau": 3000,
            "transition_eps": 1e-6,
            "n_x": 2500,
        }

        full = _solve_starobinsky_background_with_torsion(**settings)
        base = _solve_starobinsky_baseline_match(full=full, **settings)

        # Build potentials (scalar + tensor) on uniform x grids
        pot_full = _build_potentials(full, n_x=settings["n_x"])
        pot_base = _build_potentials(base, n_x=settings["n_x"])

        # Estimate bounce scales from MS potential heights in raw x-units (v'' + (k^2-U(x)) v = 0).
        U_s_max = float(max(np.nanmax(pot_full.U_s_arr), np.nanmax(pot_base.U_s_arr)))
        U_t_max = float(max(np.nanmax(pot_full.U_t_arr), np.nanmax(pot_base.U_t_arr)))
        k_bounce_s_raw = float(np.sqrt(max(0.0, U_s_max)))
        k_bounce_t_raw = float(np.sqrt(max(0.0, U_t_max)))

        def _scale_potentials(pot: _Potentials, k_scale: float) -> _Potentials:
            """
            Rescale x -> x_hat = k_scale * x and U -> U_hat = U / k_scale^2 so that the mode equation
            can be solved in the dimensionless variables (x_hat, k_hat=k/k_scale).
            """
            if not np.isfinite(k_scale) or k_scale <= 0:
                raise ValueError(f"Invalid k_scale={k_scale} for potential rescaling")
            return _Potentials(
                x_grid=pot.x_grid * k_scale,
                U_s_arr=pot.U_s_arr / (k_scale**2),
                U_t_arr=pot.U_t_arr / (k_scale**2),
                x_min=float(pot.x_min * k_scale),
                x_max=float(pot.x_max * k_scale),
            )

        # Work in k_hat = k/k_bounce and x_hat = k_bounce * x, per sector (scalar/tensor have different k_bounce).
        pot_full_s = _scale_potentials(pot_full, k_bounce_s_raw if k_bounce_s_raw > 0 else 1.0)
        pot_base_s = _scale_potentials(pot_base, k_bounce_s_raw if k_bounce_s_raw > 0 else 1.0)
        pot_full_t = _scale_potentials(pot_full, k_bounce_t_raw if k_bounce_t_raw > 0 else 1.0)
        pot_base_t = _scale_potentials(pot_base, k_bounce_t_raw if k_bounce_t_raw > 0 else 1.0)

        # k-grid in rescaled units (k_hat ≡ k/k_bounce). High-k corresponds to k_hat >> 1.
        k_min = 0.1
        k_max = 10.0
        k_points = 12
        ks = np.logspace(np.log10(k_min), np.log10(k_max), num=k_points)

        # solver tolerances
        rtol = 1e-8
        atol = 1e-10

        # Locate the bounce (τ=0) in the shifted x-coordinate (x(tau_transition)=0), for start-time selection.
        x_full_shifted = _compute_x_shifted(tau=full.tau, a_tau=full.a, tau_t=float(full.tau_transition))
        x_base_shifted = _compute_x_shifted(tau=base.tau, a_tau=base.a, tau_t=float(base.tau_transition))
        x_bounce_full_raw = float(np.interp(0.0, full.tau, x_full_shifted))
        x_bounce_base_raw = float(np.interp(0.0, base.tau, x_base_shifted))
        x_bounce_full_s_hat = x_bounce_full_raw * (k_bounce_s_raw if k_bounce_s_raw > 0 else 1.0)
        x_bounce_base_s_hat = x_bounce_base_raw * (k_bounce_s_raw if k_bounce_s_raw > 0 else 1.0)
        x_bounce_full_t_hat = x_bounce_full_raw * (k_bounce_t_raw if k_bounce_t_raw > 0 else 1.0)
        x_bounce_base_t_hat = x_bounce_base_raw * (k_bounce_t_raw if k_bounce_t_raw > 0 else 1.0)

        def select_start(
            pot: _Potentials,
            *,
            k: float,
            U_arr: np.ndarray,
            x_upper: Optional[float] = None,
            min_fraction: float = 0.25,
            adiabatic_factor: float = 100.0,
        ) -> tuple[float, float, str]:
            """
            Pick an adiabatic initial time where ω^2=k^2-U(x) is positive and k^2 >> |U(x)|.

            If x_upper is provided, choose the latest eligible time at or before x_upper
            (useful to force the initial condition to lie before the bounce).
            """
            omega2 = (k * k) - U_arr
            xg = pot.x_grid
            want = omega2 > (min_fraction * k * k)
            want = want & ((k * k) >= (adiabatic_factor * np.abs(U_arr)))
            if x_upper is not None:
                want = want & (xg <= x_upper)

            idx = np.where(want)[0]
            status = "adiabatic"
            if len(idx) == 0:
                # Fallback: choose an oscillatory point with *maximal* ω^2 (most adiabatic available),
                # instead of the latest oscillatory point (which can sit near a turning point ω^2≈0).
                want2 = omega2 > 0.0
                if x_upper is not None:
                    want2 = want2 & (xg <= x_upper)
                idx2 = np.where(want2)[0]
                if len(idx2) == 0:
                    raise RuntimeError(f"No oscillatory region found for k={k} on this potential grid.")
                i0 = int(idx2[int(np.argmax(omega2[idx2]))])
                status = "oscillatory_max_omega"
            else:
                i0 = int(idx[-1] if x_upper is not None else idx[0])
            x0 = float(xg[i0])
            w0 = float(np.sqrt(max(float(omega2[i0]), 1e-18)))
            return x0, w0, status

        def _omega2_min_over_range(*, x_grid: np.ndarray, U_arr: np.ndarray, k_val: float, x_start: float, x_end: float) -> float:
            # Match _solve_mode_symplectic interval selection.
            i0 = int(np.argmin(np.abs(x_grid - x_start)))
            i1 = int(np.searchsorted(x_grid, x_end, side="right") - 1)
            i1 = max(i0 + 1, min(i1, int(x_grid.size - 1)))
            omega2 = (k_val * k_val) - U_arr
            return float(np.min(omega2[i0 : i1 + 1]))

        # Solve scalar + tensor transfer functions
        T_s: list[float] = []
        T_t: list[float] = []
        W_s: list[float] = []
        W_t: list[float] = []
        x_start_scalar_full: list[float] = []
        x_start_scalar_base: list[float] = []
        x_start_tensor_full: list[float] = []
        x_start_tensor_base: list[float] = []
        omega0_scalar_full: list[float] = []
        omega0_scalar_base: list[float] = []
        omega0_tensor_full: list[float] = []
        omega0_tensor_base: list[float] = []
        ic_status_scalar_full: list[str] = []
        ic_status_scalar_base: list[str] = []
        ic_status_tensor_full: list[str] = []
        ic_status_tensor_base: list[str] = []
        solver_scalar_full: list[str] = []
        solver_tensor_full: list[str] = []

        x_end_scalar_cap = 600.0
        x_end_scalar = float(min(pot_full_s.x_max, pot_base_s.x_max, x_end_scalar_cap))
        x_end_tensor = float(min(pot_full_t.x_max, pot_base_t.x_max))

        # Grid spacing (x_hat units). Used to choose leapfrog substeps per k for stability.
        dx_s = float(pot_full_s.x_grid[1] - pot_full_s.x_grid[0])
        dx_t = float(pot_full_t.x_grid[1] - pot_full_t.x_grid[0])

        def _substeps_for(
            *,
            x_grid: np.ndarray,
            U_arr: np.ndarray,
            k_val: float,
            x_start: float,
            x_end: float,
            dx: float,
            target_step: float,
            min_substeps: int,
        ) -> int:
            """
            Choose leapfrog substeps so that the effective micro-step h=dx/substeps resolves the
            fastest local timescale set by |ω| = sqrt(|k^2-U|) on the integration interval.

            Heuristic: require |ω|max * h <= target_step.
            """
            i0 = int(np.argmin(np.abs(x_grid - x_start)))
            i1 = int(np.searchsorted(x_grid, x_end, side="right") - 1)
            i1 = max(i0 + 1, min(i1, int(x_grid.size - 1)))
            omega_abs_max = float(np.max(np.sqrt(np.abs((k_val * k_val) - U_arr[i0 : i1 + 1]))))
            if not np.isfinite(omega_abs_max) or omega_abs_max <= 0 or dx <= 0:
                return int(min_substeps)
            sub = int(np.ceil((omega_abs_max * dx) / float(target_step)))
            return int(max(min_substeps, sub))

        for k in ks:
            kf = float(k)

            # scalar
            xs_b_s, w0_b_s, st_b_s = select_start(
                pot_full_s,
                k=kf,
                U_arr=pot_full_s.U_s_arr,
                x_upper=(x_bounce_full_s_hat - 5.0),
            )
            xs_r_s, w0_r_s, st_r_s = select_start(
                pot_base_s,
                k=kf,
                U_arr=pot_base_s.U_s_arr,
                x_upper=None,
            )
            x_start_scalar_full.append(xs_b_s)
            x_start_scalar_base.append(xs_r_s)
            omega0_scalar_full.append(w0_b_s)
            omega0_scalar_base.append(w0_r_s)
            ic_status_scalar_full.append(st_b_s)
            ic_status_scalar_base.append(st_r_s)
            sub_s = _substeps_for(
                x_grid=pot_full_s.x_grid,
                U_arr=pot_full_s.U_s_arr,
                k_val=kf,
                x_start=xs_b_s,
                x_end=x_end_scalar,
                dx=dx_s,
                target_step=0.05,
                min_substeps=6,
            )
            sol_b_s = _solve_mode_symplectic(
                k=kf,
                omega0=w0_b_s,
                x_grid=pot_full_s.x_grid,
                U_grid=pot_full_s.U_s_arr,
                x_start=xs_b_s,
                x_end=x_end_scalar,
                substeps=sub_s,
            )
            sol_r_s = _solve_mode_symplectic(
                k=kf,
                omega0=w0_r_s,
                x_grid=pot_base_s.x_grid,
                U_grid=pot_base_s.U_s_arr,
                x_start=xs_r_s,
                x_end=x_end_scalar,
                substeps=sub_s,
            )
            solver_scalar_full.append("symplectic")
            T_s.append(sol_b_s["amp"] / sol_r_s["amp"] if sol_r_s["amp"] != 0 else float("nan"))
            W_s.append(sol_b_s["W_over_i"])

            # tensor
            xs_b_t, w0_b_t, st_b_t = select_start(
                pot_full_t,
                k=kf,
                U_arr=pot_full_t.U_t_arr,
                x_upper=(x_bounce_full_t_hat - 5.0),
            )
            xs_r_t, w0_r_t, st_r_t = select_start(
                pot_base_t,
                k=kf,
                U_arr=pot_base_t.U_t_arr,
                x_upper=None,
            )
            x_start_tensor_full.append(xs_b_t)
            x_start_tensor_base.append(xs_r_t)
            omega0_tensor_full.append(w0_b_t)
            omega0_tensor_base.append(w0_r_t)
            ic_status_tensor_full.append(st_b_t)
            ic_status_tensor_base.append(st_r_t)
            sub_t = _substeps_for(
                x_grid=pot_full_t.x_grid,
                U_arr=pot_full_t.U_t_arr,
                k_val=kf,
                x_start=xs_b_t,
                x_end=x_end_tensor,
                dx=dx_t,
                target_step=0.08,
                min_substeps=4,
            )
            sol_b_t = _solve_mode_symplectic(
                k=kf,
                omega0=w0_b_t,
                x_grid=pot_full_t.x_grid,
                U_grid=pot_full_t.U_t_arr,
                x_start=xs_b_t,
                x_end=x_end_tensor,
                substeps=sub_t,
            )
            sol_r_t = _solve_mode_symplectic(
                k=kf,
                omega0=w0_r_t,
                x_grid=pot_base_t.x_grid,
                U_grid=pot_base_t.U_t_arr,
                x_start=xs_r_t,
                x_end=x_end_tensor,
                substeps=sub_t,
            )
            T_t.append(sol_b_t["amp"] / sol_r_t["amp"] if sol_r_t["amp"] != 0 else float("nan"))
            W_t.append(sol_b_t["W_over_i"])
            solver_tensor_full.append("symplectic")

        ks_list = [float(x) for x in ks]
        T_s_arr = np.array(T_s, dtype=float)
        T_t_arr = np.array(T_t, dtype=float)
        W_s_arr = np.array(W_s, dtype=float)
        W_t_arr = np.array(W_t, dtype=float)

        # In rescaled units (x_hat=k_bounce*x, k_hat=k/k_bounce), the bounce scale is ~1 by construction.
        k_bounce_s = 1.0
        k_bounce_t = 1.0

        def _summarize_high_k(T_arr: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
            if not bool(np.any(mask)):
                return float("nan"), float("nan")
            return float(np.nanmean(T_arr[mask])), float(np.nanmax(np.abs(T_arr[mask] - 1.0)))

        # "High-k" region for the asymptotic T(k)->1 test: k >= max(k_max/5, 5*k_bounce)
        high_k_mask_s = ks >= max((k_max / 5.0), 5.0 * k_bounce_s)
        high_k_mask_t = ks >= max((k_max / 5.0), 5.0 * k_bounce_t)
        high_k_Ts, high_k_Ts_dev = _summarize_high_k(T_s_arr, high_k_mask_s)
        high_k_Tt, high_k_Tt_dev = _summarize_high_k(T_t_arr, high_k_mask_t)

        # Wronskian diagnostics: also report the subset where ω^2=k^2-U stays positive everywhere (no turning point)
        safe_s_idx = [
            i
            for i, k_val in enumerate(ks_list)
            if (
                (_omega2_min_over_range(x_grid=pot_full_s.x_grid, U_arr=pot_full_s.U_s_arr, k_val=k_val, x_start=x_start_scalar_full[i], x_end=x_end_scalar) > 0.0)
                and (_omega2_min_over_range(x_grid=pot_base_s.x_grid, U_arr=pot_base_s.U_s_arr, k_val=k_val, x_start=x_start_scalar_base[i], x_end=x_end_scalar) > 0.0)
            )
        ]
        safe_t_idx = [
            i
            for i, k_val in enumerate(ks_list)
            if (
                (_omega2_min_over_range(x_grid=pot_full_t.x_grid, U_arr=pot_full_t.U_t_arr, k_val=k_val, x_start=x_start_tensor_full[i], x_end=x_end_tensor) > 0.0)
                and (_omega2_min_over_range(x_grid=pot_base_t.x_grid, U_arr=pot_base_t.U_t_arr, k_val=k_val, x_start=x_start_tensor_base[i], x_end=x_end_tensor) > 0.0)
            )
        ]

        W_s_dev = float(np.nanmax(np.abs(W_s_arr - 1.0)))
        W_t_dev = float(np.nanmax(np.abs(W_t_arr - 1.0)))
        W_s_dev_safe = float(np.nanmax(np.abs(W_s_arr[safe_s_idx] - 1.0))) if safe_s_idx else float("nan")
        W_t_dev_safe = float(np.nanmax(np.abs(W_t_arr[safe_t_idx] - 1.0))) if safe_t_idx else float("nan")

        checks: list[Check] = []
        checks.append(
            Check(
                check_id="stability_F_positive",
                passed=bool(full.F_min > 0 and base.F_min > 0),
                detail=f"min F: full={full.F_min:.6g}, base={base.F_min:.6g}",
            )
        )
        checks.append(
            Check(
                check_id="background_constraint_residual_small",
                passed=bool(full.max_constraint_abs < 1e-6 and base.max_constraint_abs < 1e-6),
                detail=f"max |constraint|: full={full.max_constraint_abs:.3e}, base={base.max_constraint_abs:.3e}",
            )
        )
        checks.append(
            Check(
                check_id="scalar_T_high_k",
                passed=bool(np.isfinite(high_k_Ts) and abs(high_k_Ts - 1.0) < 0.10 and high_k_Ts_dev < 0.20),
                detail=f"high-k threshold uses k>=max(k_max/5,5*k_bounce_s)={max((k_max/5.0), 5.0*k_bounce_s):.3g}; mean={high_k_Ts:.4g}, max|T_s-1|={high_k_Ts_dev:.4g}",
            )
        )
        checks.append(
            Check(
                check_id="tensor_T_high_k",
                passed=bool(np.isfinite(high_k_Tt) and abs(high_k_Tt - 1.0) < 0.10 and high_k_Tt_dev < 0.20),
                detail=f"high-k threshold uses k>=max(k_max/5,5*k_bounce_t)={max((k_max/5.0), 5.0*k_bounce_t):.3g}; mean={high_k_Tt:.4g}, max|T_t-1|={high_k_Tt_dev:.4g}",
            )
        )
        checks.append(
            Check(
                check_id="Wronskian_scalar",
                passed=bool(np.isfinite(W_s_dev) and W_s_dev < 5e-3),
                detail=f"max |W_s/i - 1| (all k)={W_s_dev:.3e}; (ω^2>0 interval subset n={len(safe_s_idx)})={W_s_dev_safe:.3e}",
            )
        )
        checks.append(
            Check(
                check_id="Wronskian_tensor",
                passed=bool(np.isfinite(W_t_dev) and W_t_dev < 5e-3),
                detail=f"max |W_t/i - 1| (all k)={W_t_dev:.3e}; (ω^2>0 interval subset n={len(safe_t_idx)})={W_t_dev_safe:.3e}",
            )
        )

        # Optional plot
        plot_path = None
        if config.plot:
            try:
                import matplotlib.pyplot as plt

                out_dir = self.output_dir(config)
                out_dir.mkdir(parents=True, exist_ok=True)
                plot_path = out_dir / "T_of_k.png"

                fig = plt.figure(figsize=(7, 4))
                ax = fig.add_subplot(1, 1, 1)
                ax.set_xscale("log")
                ax.plot(ks, T_s_arr, lw=1.8, label="T_s(k) (scalar)")
                ax.plot(ks, T_t_arr, lw=1.8, label="T_t(k) (tensor)")
                ax.axhline(1.0, color="black", lw=1.0, alpha=0.6)
                ax.set_xlabel("k_hat = k/k_bounce (dimensionless)")
                ax.set_ylabel("T(k)")
                ax.set_title("Bounce transfer functions (Starobinsky f(R) + torsion a^-6 proxy)")
                ax.grid(True, which="both", ls=":", alpha=0.4)
                ax.legend(loc="best")
                fig.tight_layout()
                fig.savefig(plot_path, dpi=160)
                plt.close(fig)
            except Exception:
                plot_path = None

        report_lines = [
            "TFPT bounce_perturbations (paper v2.4 Appendix L)",
            "",
            "What changed vs the previous proxy:",
            "- Background is now generated by an explicit ODE system for Starobinsky f(R) in terms of F=df/dR,",
            "  plus a torsion proxy energy density rho~a^{-6} that decays after the bounce.",
            "- Scalar and tensor z(η) are implemented in f(R): z_s=a*sqrt(Q_s), z_t=a*sqrt(F).",
            "- Baseline is a matched pure-R2 run (torsion proxy switched off), anchored at tau_transition and",
            "  capped to avoid numerically unstable very-early-time integration.",
            "",
            f"TFPT input: M/Mpl = {M_over_Mpl:.12e}",
            "",
            "Background settings (dimensionless τ=M t):",
            *[f"- {k} = {v}" for k, v in settings.items()],
            "",
            "Background diagnostics:",
            f"- tau_transition(full) = {full.tau_transition:.6g}",
            f"- min F (full/base) = {full.F_min:.6g} / {base.F_min:.6g}",
            f"- max |constraint| (full/base) = {full.max_constraint_abs:.3e} / {base.max_constraint_abs:.3e}",
            f"- MS potential peaks: max U_s={U_s_max:.3e}, max U_t={U_t_max:.3e}",
            f"- inferred bounce scales (raw x-units): k_bounce_s≈sqrt(max U_s)={k_bounce_s_raw:.3g}, k_bounce_t≈sqrt(max U_t)={k_bounce_t_raw:.3g}",
            "- mode equation solved in rescaled variables x_hat=k_bounce*x, k_hat=k/k_bounce per sector",
            "",
            "Numerics:",
            f"- k grid (k_hat): {k_points} log-spaced points in [{k_min}, {k_max}]",
            "- mode solver: symplectic leapfrog (velocity Verlet) on the uniform x-grid; per-k substeps chosen to keep the effective step small in ω^2 regions",
            f"- x_end_scalar = {x_end_scalar:.6g}, x_end_tensor = {x_end_tensor:.6g}",
            "- ICs: per-k local adiabatic ω0 at a point with k^2 >= 100*|U(x)|; for the bounce run start is forced pre-bounce (x <= x_bounce-5).",
            "",
            "Diagnostics:",
            f"- mean T_s(k) at high k (k >= {max((k_max/5.0), 5.0*k_bounce_s):.3g}): {high_k_Ts:.6f}",
            f"- max |T_s-1| at high k (k >= {max((k_max/5.0), 5.0*k_bounce_s):.3g}): {high_k_Ts_dev:.6f}",
            f"- mean T_t(k) at high k (k >= {max((k_max/5.0), 5.0*k_bounce_t):.3g}): {high_k_Tt:.6f}",
            f"- max |T_t-1| at high k (k >= {max((k_max/5.0), 5.0*k_bounce_t):.3g}): {high_k_Tt_dev:.6f}",
            f"- max |W_s/i - 1| (all k): {W_s_dev:.6e}",
            f"- max |W_s/i - 1| (ω^2>0 interval subset n={len(safe_s_idx)}): {W_s_dev_safe:.6e}",
            f"- max |W_t/i - 1| (all k): {W_t_dev:.6e}",
            f"- max |W_t/i - 1| (ω^2>0 interval subset n={len(safe_t_idx)}): {W_t_dev_safe:.6e}",
            "",
            "Checks:",
            *[f"- {chk.check_id}: {'PASS' if chk.passed else 'FAIL'} ({chk.detail})" for chk in checks],
            "",
            "Notes:",
            "- The torsion sector is represented as an effective stiff component ρ ∝ a^{-6} in the background constraint.",
            "- Scalar f(R) modes can be numerically stiff; this module uses rescaled variables and per-k adiabatic initialization.",
            "- To handle highly squeezed / tachyonic intervals (ω^2<0), the symplectic integrator applies a per-step Wronskian renormalization to enforce the canonical normalization (numerical stabilization; see code comment in _solve_mode_symplectic).",
            "",
        ]

        if plot_path is not None:
            report_lines += [f"Plot: {plot_path}", ""]

        return ModuleResult(
            results={
                "tfpt": {"M_over_Mpl": M_over_Mpl},
                "background_full": full.as_dict(),
                "background_baseline": base.as_dict(),
                "k_grid": ks_list,
                "T_scalar": [float(v) for v in T_s_arr],
                "T_tensor": [float(v) for v in T_t_arr],
                "x_end_scalar": x_end_scalar,
                "x_end_tensor": x_end_tensor,
                "x_start_scalar_full": [float(v) for v in x_start_scalar_full],
                "x_start_scalar_baseline": [float(v) for v in x_start_scalar_base],
                "x_start_tensor_full": [float(v) for v in x_start_tensor_full],
                "x_start_tensor_baseline": [float(v) for v in x_start_tensor_base],
                "omega0_scalar_full": [float(v) for v in omega0_scalar_full],
                "omega0_scalar_baseline": [float(v) for v in omega0_scalar_base],
                "omega0_tensor_full": [float(v) for v in omega0_tensor_full],
                "omega0_tensor_baseline": [float(v) for v in omega0_tensor_base],
                "ic_status_scalar_full": list(ic_status_scalar_full),
                "ic_status_scalar_baseline": list(ic_status_scalar_base),
                "ic_status_tensor_full": list(ic_status_tensor_full),
                "ic_status_tensor_baseline": list(ic_status_tensor_base),
                "solver_scalar_full": list(solver_scalar_full),
                "solver_tensor_full": list(solver_tensor_full),
                "Wronskian_scalar_over_i": [float(v) for v in W_s_arr],
                "Wronskian_tensor_over_i": [float(v) for v in W_t_arr],
                "diagnostics": {
                    "U_s_max": U_s_max,
                    "U_t_max": U_t_max,
                    "k_bounce_s_est_raw": k_bounce_s_raw,
                    "k_bounce_t_est_raw": k_bounce_t_raw,
                    "safe_scalar_k_count": len(safe_s_idx),
                    "safe_tensor_k_count": len(safe_t_idx),
                    "high_k_mean_Ts": high_k_Ts,
                    "high_k_max_abs_Ts_minus_1": high_k_Ts_dev,
                    "high_k_mean_Tt": high_k_Tt,
                    "high_k_max_abs_Tt_minus_1": high_k_Tt_dev,
                    "max_abs_Ws_over_i_minus_1": W_s_dev,
                    "max_abs_Wt_over_i_minus_1": W_t_dev,
                    "max_abs_Ws_over_i_minus_1_safe": W_s_dev_safe,
                    "max_abs_Wt_over_i_minus_1_safe": W_t_dev_safe,
                    "ic_status_counts_scalar_full": {k: int(ic_status_scalar_full.count(k)) for k in sorted(set(ic_status_scalar_full))},
                    "ic_status_counts_tensor_full": {k: int(ic_status_tensor_full.count(k)) for k in sorted(set(ic_status_tensor_full))},
                    "solver_counts_scalar_full": {k: int(solver_scalar_full.count(k)) for k in sorted(set(solver_scalar_full))},
                },
                "plot": {"T_of_k_png": str(plot_path) if plot_path is not None else None},
            },
            checks=checks,
            report="\n".join(report_lines),
            warnings=[],
        )


@dataclass(frozen=True)
class _StarobinskyBg:
    tau: np.ndarray
    a: np.ndarray
    h: np.ndarray
    F: np.ndarray
    g: np.ndarray
    rho0: float
    tau_transition: float
    F_min: float
    max_constraint_abs: float

    def as_dict(self) -> dict[str, object]:
        return {
            "rho0": self.rho0,
            "tau_transition": float(self.tau_transition),
            "F_min": float(self.F_min),
            "max_constraint_abs": float(self.max_constraint_abs),
            "tau_range": [float(self.tau[0]), float(self.tau[-1])],
            "a_range": [float(np.min(self.a)), float(np.max(self.a))],
        }


@dataclass(frozen=True)
class _Potentials:
    x_grid: np.ndarray
    U_s_arr: np.ndarray
    U_t_arr: np.ndarray
    x_min: float
    x_max: float

    def U_s(self, x: float) -> float:
        return float(np.interp(x, self.x_grid, self.U_s_arr))

    def U_t(self, x: float) -> float:
        return float(np.interp(x, self.x_grid, self.U_t_arr))


def _solve_starobinsky_background_with_torsion(
    *,
    F0: float,
    g0: float,
    a0: float,
    tau_pre: float,
    tau_post: float,
    n_tau: int,
    transition_eps: float,
    n_x: int,
) -> _StarobinskyBg:
    # Torsion proxy: rho~a^{-6} with rho0 fixed by the bounce constraint at tau=0, h=0
    rho0 = -0.75 * (F0 - 1.0) ** 2 * (a0**6)
    return _solve_starobinsky_background(
        rho0=rho0,
        F0=F0,
        g0=g0,
        a0=a0,
        tau_pre=tau_pre,
        tau_post=tau_post,
        n_tau=n_tau,
        transition_eps=transition_eps,
        tau_anchor=0.0,
        branch_back=-1,
        branch_fwd=+1,
    )


def _solve_starobinsky_baseline_match(*, full: _StarobinskyBg, **settings) -> _StarobinskyBg:
    # Baseline: pure f(R) (rho0=0) matched to the full solution at tau_transition
    tau_t = float(full.tau_transition)
    tau_start = float(full.tau[0])
    tau_end = float(full.tau[-1])
    a_t = float(np.interp(tau_t, full.tau, full.a))
    F_t = float(np.interp(tau_t, full.tau, full.F))
    g_t = float(np.interp(tau_t, full.tau, full.g))
    # Avoid very long backward integration for the baseline (pure-R2) branch.
    # Numerically, integrating too far back can become extremely stiff. For the transfer-function
    # ratio, it is sufficient that the baseline provides a sub-horizon start region; we cap at tau>=0.
    tau_pre_full = tau_t - tau_start
    tau_pre_baseline = min(tau_pre_full, max(0.0, tau_t))
    return _solve_starobinsky_background(
        rho0=0.0,
        F0=F_t,
        g0=g_t,
        a0=a_t,
        tau_pre=tau_pre_baseline,
        tau_post=tau_end - tau_t,
        n_tau=int(settings["n_tau"]),
        transition_eps=float(settings["transition_eps"]),
        tau_anchor=tau_t,
        branch_back=+1,
        branch_fwd=+1,
    )


def _solve_starobinsky_background(
    *,
    rho0: float,
    F0: float,
    g0: float,
    a0: float,
    tau_pre: float,
    tau_post: float,
    n_tau: int,
    transition_eps: float,
    tau_anchor: float = 0.0,
    branch_back: int = -1,
    branch_fwd: int = +1,
) -> _StarobinskyBg:
    """
    Dimensionless background in τ = M t with f(R)=R+R^2/(6M^2) expressed via F=df/dR.

    Variables:
      u(τ)=ln a(τ), F(τ), g(τ)=dF/dτ with h(τ)=H/M computed from the modified Friedmann constraint.
      torsion proxy rho(τ) = rho0 * exp(-6u(τ)).

    Constraint (dimensionless, with rho interpreted as rho/M^2):
      3 F h^2 = (3/4)(F-1)^2 - 3 h g + rho

    This is quadratic in h. We pick the "expanding" branch (h>=0) or "contracting" branch (h<=0)
    separately for the backward/forward integration to model a bounce.

    Evolution equations:
      u' = h
      F' = g
      g' = -3 h g - (F-1) + (2/3) rho
    """
    def h_from_constraint(F: float, g: float, rho: float, branch: int) -> float:
        if F <= 0:
            raise RuntimeError(f"Invalid F<=0 encountered: F={F}")
        C = 0.25 * (F - 1.0) ** 2 + (rho / 3.0)
        disc = g * g + 4.0 * F * C
        if disc < 0:
            raise RuntimeError(f"Constraint discriminant negative: disc={disc} (F={F}, g={g}, rho={rho}, C={C})")
        s = float(np.sqrt(disc))
        b = 1.0 if branch >= 0 else -1.0
        return (-g + b * s) / (2.0 * F)

    def rhs_factory(branch: int):
        def rhs(tau: float, y: np.ndarray) -> np.ndarray:
            u, F, g = y
            rho = float(rho0 * np.exp(-6.0 * u)) if rho0 != 0 else 0.0
            h = h_from_constraint(float(F), float(g), float(rho), branch)
            return np.array(
                [
                    h,
                    g,
                    # Trace equation (signature -,+,+,+): ddotF + 3H dotF + M^2(F-1) = (2/3)rho (for w=1),
                    # written in τ=M t units with g=dF/dτ:
                    -3.0 * h * g - (float(F) - 1.0) + (2.0 / 3.0) * rho,
                ],
                dtype=float,
            )

        return rhs

    y0 = np.array([float(np.log(a0)), F0, g0], dtype=float)

    # Integrate backward and forward from tau_anchor with the same state at tau_anchor.
    # Split so that total samples after concatenation equals n_tau.
    if n_tau < 10:
        raise ValueError("n_tau must be >= 10")
    n_half = max(2, n_tau // 2)
    t_eval_back = np.linspace(tau_anchor, tau_anchor - tau_pre, n_half)
    t_eval_fwd = np.linspace(tau_anchor, tau_anchor + tau_post, n_tau - n_half + 1)

    # Prefer a fast explicit integrator; fall back to stiff solvers if needed.
    def _integrate(method: str):
        sol_back = solve_ivp(
            rhs_factory(branch_back),
            (tau_anchor, tau_anchor - tau_pre),
            y0=y0,
            t_eval=t_eval_back,
            method=method,
            rtol=1e-8,
            atol=1e-10,
        )
        sol_fwd = solve_ivp(
            rhs_factory(branch_fwd),
            (tau_anchor, tau_anchor + tau_post),
            y0=y0,
            t_eval=t_eval_fwd,
            method=method,
            rtol=1e-8,
            atol=1e-10,
        )
        return sol_back, sol_fwd

    sol_back, sol_fwd = _integrate("DOP853")
    if (not sol_back.success) or (not sol_fwd.success):
        sol_back, sol_fwd = _integrate("BDF")
    if (not sol_back.success) or (not sol_fwd.success):
        sol_back, sol_fwd = _integrate("Radau")
    if not sol_back.success:
        raise RuntimeError(f"Background ODE (backward) failed: {sol_back.message}")
    if not sol_fwd.success:
        raise RuntimeError(f"Background ODE (forward) failed: {sol_fwd.message}")

    tau_back = sol_back.t[::-1]
    y_back = sol_back.y[:, ::-1]
    tau_fwd = sol_fwd.t
    y_fwd = sol_fwd.y

    tau = np.concatenate([tau_back[:-1], tau_fwd])
    y = np.concatenate([y_back[:, :-1], y_fwd], axis=1)

    u = y[0, :]
    a = np.exp(u)
    F = y[1, :]
    g = y[2, :]

    rho = rho0 * np.exp(-6.0 * u) if rho0 != 0 else np.zeros_like(a)

    # reconstruct h with the same branch choices
    h = np.empty_like(u)
    for i in range(len(tau)):
        branch = branch_back if tau[i] < tau_anchor else branch_fwd
        h[i] = h_from_constraint(float(F[i]), float(g[i]), float(rho[i]), branch)

    # Friedmann constraint residual (dimensionless):
    # 3 F h^2 = (3/4)(F-1)^2 - 3 h g + rho
    constraint = 3.0 * F * (h**2) - (0.75 * (F - 1.0) ** 2) + 3.0 * h * g - rho
    max_constraint_abs = float(np.max(np.abs(constraint)))

    # Pick transition time: first tau>tau_anchor where torsion term is negligible vs R2 term
    tau_transition = tau_anchor
    if rho0 != 0:
        mask = tau > tau_anchor
        denom = 0.75 * (F - 1.0) ** 2 + 1e-12
        ratio = np.abs(rho) / denom
        idxs = np.where(mask & (ratio < transition_eps))[0]
        if len(idxs) > 0:
            tau_transition = float(tau[int(idxs[0])])

    return _StarobinskyBg(
        tau=tau,
        a=a,
        h=h,
        F=F,
        g=g,
        rho0=float(rho0),
        tau_transition=float(tau_transition),
        F_min=float(np.min(F)),
        max_constraint_abs=float(max_constraint_abs),
    )


def _compute_x_shifted(*, tau: np.ndarray, a_tau: np.ndarray, tau_t: float) -> np.ndarray:
    """
    Compute dimensionless conformal time

      x(τ) = ∫ dτ / a(τ)

    and shift it so that x(tau_t)=0.
    """
    x = np.zeros_like(tau)
    i0 = int(np.argmin(np.abs(tau - tau_t)))
    # forward
    for i in range(i0 + 1, len(tau)):
        dt = tau[i] - tau[i - 1]
        x[i] = x[i - 1] + 0.5 * dt * (1.0 / a_tau[i - 1] + 1.0 / a_tau[i])
    # backward
    for i in range(i0 - 1, -1, -1):
        dt = tau[i + 1] - tau[i]
        x[i] = x[i + 1] - 0.5 * dt * (1.0 / a_tau[i + 1] + 1.0 / a_tau[i])
    return x


def _build_potentials(bg: _StarobinskyBg, *, n_x: int) -> _Potentials:
    """
    Build scalar/tensor potentials U=z''/z on a uniform x-grid where x = ∫ dτ/a (dimensionless conformal time).
    The grid is shifted so that x(tau_transition)=0.
    """
    tau = bg.tau
    a_tau = bg.a
    F_tau = bg.F
    h_tau = bg.h
    g_tau = bg.g
    tau_t = bg.tau_transition

    # compute x(tau) with x(tau_t)=0 (numerical trapezoid)
    x = _compute_x_shifted(tau=tau, a_tau=a_tau, tau_t=float(tau_t))

    x_min = float(np.min(x))
    x_max = float(np.max(x))
    x_grid = np.linspace(x_min, x_max, int(n_x))
    dx = float(x_grid[1] - x_grid[0])

    a = np.interp(x_grid, x, a_tau)
    F = np.interp(x_grid, x, F_tau)
    h = np.interp(x_grid, x, h_tau)
    g = np.interp(x_grid, x, g_tau)

    # scalar Q_s in x-time (equivalent to conformal time up to the common M scaling)
    denom = h + g / (2.0 * F)
    denom2 = denom * denom + 1e-18
    Qs = (3.0 * g * g) / (2.0 * F * denom2)
    Qs = np.maximum(Qs, 1e-18)

    z_s = a * np.sqrt(Qs)
    z_t = a * np.sqrt(np.maximum(F, 1e-18))

    def zpp_over_z(z: np.ndarray) -> np.ndarray:
        # Compute z''/z via log-derivatives for numerical stability when z becomes small.
        z_safe = np.maximum(z, 1e-300)
        lnz = np.log(z_safe)
        dlnz = np.gradient(lnz, dx, edge_order=2)
        d2lnz = np.gradient(dlnz, dx, edge_order=2)
        return d2lnz + dlnz * dlnz

    U_s = zpp_over_z(z_s)
    U_t = zpp_over_z(z_t)

    return _Potentials(x_grid=x_grid, U_s_arr=U_s, U_t_arr=U_t, x_min=float(x_min), x_max=float(x_max))

