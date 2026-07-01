"""
exOptimalTrajectoryCasadi.py

Minimum-time trajectory generation using CasADi.

Kinematic model:
    x_n_dot          = U_g cos(chi)
    y_n_dot          = U_g sin(chi)
    chi_dot          = omega_chi
    U_g_dot          = a_g
    a_g_dot          = u_1
    omega_chi_dot    = u_2

The cost function minimizes final time while penalizing large jerk u_1 and
course angular acceleration u_2. The script generates the plots used in
the trajectory-optimization example.

Requirements:
    pip install casadi numpy matplotlib

Author: Thor I. Fossen
Date:   2026-06-29
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import casadi as ca
except ImportError as exc:
    raise ImportError("CasADi is required. Install it with: pip install casadi") from exc


def main() -> None:
    # ------------------------------------------------------------------
    # Discretization and weights
    # ------------------------------------------------------------------
    N = 200     # Number of control intervals (shooting intervals) 
    w1 = 1.0    # Final-time weight
    w2 = 0.05   # Jerk weight
    w3 = 0.05   # Angular-acceleration weight

    L = 2.0     # USV length [m]

    # ------------------------------------------------------------------
    # Boundary conditions
    # State vector: X = [x_n, y_n, chi, U_g, a_g, omega_chi]^T
    # ------------------------------------------------------------------
    x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([30.0, 20.0, np.deg2rad(10.0), 1.0, 0.0, 0.0])

    tf_min, tf_max = 0.1, 1000.0
    Ug_max = 2.0                    # Maximum ground speed [m/s]
    ag_min, ag_max = -0.4, 0.4
    omega_min, omega_max = np.deg2rad(-3.0), np.deg2rad(3.0)
    u1_min, u1_max = -0.5, 0.5
    u2_min, u2_max = -0.5, 0.5

    # Curvature constraint for a rudder-controlled USV:
    # omega_chi = kappa_chi U_g, |kappa_chi| <= kappa_max.
    R_min = 3.0 * L                 # Minimum turning radius [m]
    kappa_max = 1.0 / R_min         # [1/m]

    # ------------------------------------------------------------------
    # CasADi optimization problem
    # ------------------------------------------------------------------
    opti = ca.Opti()

    X = opti.variable(6, N + 1)      # x_n, y_n, chi, U_g, a_g, omega_chi
    U = opti.variable(2, N)          # u_1: jerk, u_2: angular acceleration
    tf = opti.variable()             # Final time
    dt = tf / N

    def f(x, u):
        chi = x[2]
        Ug = x[3]
        ag = x[4]
        omega_chi = x[5]
        u1 = u[0]
        u2 = u[1]
        return ca.vertcat(
            Ug * ca.cos(chi),
            Ug * ca.sin(chi),
            omega_chi,
            ag,
            u1,
            u2,
        )

    # Dynamics: direct multiple shooting with RK4 integration
    for k in range(N):
        xk = X[:, k]
        uk = U[:, k]

        k1 = f(xk, uk)
        k2 = f(xk + 0.5 * dt * k1, uk)
        k3 = f(xk + 0.5 * dt * k2, uk)
        k4 = f(xk + dt * k3, uk)

        x_next = xk + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)
        opti.subject_to(X[:, k + 1] == x_next)

    # Boundary conditions
    opti.subject_to(X[:, 0] == x0)
    opti.subject_to(X[:, -1] == xf)

    opti.subject_to(opti.bounded(tf_min, tf, tf_max))

    # Zero-speed rotation is prevented by the curvature constraint below.
    opti.subject_to(opti.bounded(0.0, X[3, :], Ug_max))

    opti.subject_to(opti.bounded(ag_min, X[4, :], ag_max))

    # Absolute course-rate bound.
    opti.subject_to(opti.bounded(omega_min, X[5, :], omega_max))

    # Speed-dependent course-rate bound.
    opti.subject_to(X[5, :] <=  kappa_max * X[3, :])
    opti.subject_to(X[5, :] >= -kappa_max * X[3, :])

    opti.subject_to(opti.bounded(u1_min, U[0, :], u1_max))
    opti.subject_to(opti.bounded(u2_min, U[1, :], u2_max))

    # Objective: minimum time with smoothness penalties
    J = w1 * tf + w2 * dt * ca.sumsqr(U[0, :]) + w3 * dt * ca.sumsqr(U[1, :])
    opti.minimize(J)

    # Initial guess
    s = np.linspace(0.0, 1.0, N + 1)
    X_guess = np.outer(x0, 1 - s) + np.outer(xf, s)

    # Feasible speed guess: accelerate from rest to final speed.
    X_guess[3, :] = np.linspace(x0[3], xf[3], N + 1)
    X_guess[3, :] = np.maximum(X_guess[3, :], 0.1)
    X_guess[3, 0] = x0[3]

    # Initial guess
    opti.set_initial(X, X_guess)
    opti.set_initial(U, 0.0)
    opti.set_initial(tf, 25.0)

    # Solver settings
    p_opts = {"expand": True}
    s_opts = {
        "max_iter": 1000,
        "tol": 1e-8,
        "print_level": 0,
        "sb": "yes",
    }
    opti.solver("ipopt", p_opts, s_opts)

    sol = opti.solve()

    X_opt = sol.value(X)
    U_opt = sol.value(U)
    tf_opt = float(sol.value(tf))

    t = np.linspace(0.0, tf_opt, N + 1)
    tu = np.linspace(0.0, tf_opt, N)

    x_n = X_opt[0, :]
    y_n = X_opt[1, :]
    chi = X_opt[2, :]
    Ug = X_opt[3, :]
    ag = X_opt[4, :]
    omega_chi = X_opt[5, :]
    u1 = U_opt[0, :]
    u2 = U_opt[1, :]

    print(f"USV length: {L:.2f} m")
    print(f"Minimum turning radius: {R_min:.2f} m")
    print(f"Maximum curvature: {kappa_max:.3f} 1/m")
    print(f"Optimal final time: {tf_opt:.3f} s")
    print(f"Cost: {float(sol.value(J)):.6f}")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    plt.rcParams.update({
        "lines.linewidth": 1.0,
        "axes.grid": True,
        "grid.linestyle": ":",
        "grid.linewidth": 0.6,
        "font.size": 10,
    })

    fig = plt.figure(figsize=(10, 8))

    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(y_n, x_n, "k-", linewidth=2, label="Trajectory")
    ax1.plot(y_n[0], x_n[0], "ko", markersize=6, label="Start")
    ax1.plot(y_n[-1], x_n[-1], "kx", markersize=7, label="Final")
    ax1.set_xlabel(r"$y^n$ [m]")
    ax1.set_ylabel(r"$x^n$ [m]")
    ax1.set_title("Ground track")
    ax1.axis("equal")
    ax1.legend(loc="best")

    ax2 = plt.subplot(3, 2, 2)
    ax2.plot(t, np.rad2deg(chi), "k-", linewidth=2)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel(r"$\chi$ [deg]")
    ax2.set_title("Course")

    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(t, Ug, "k-", linewidth=2, label=r"$U_g$")
    ax3.axhline(Ug_max, color="k", linestyle="--", linewidth=1,
                label=r"$U_g^{\max}$")
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel(r"$U_g$ [m/s]")
    ax3.set_title("Speed over ground")
    ax3.legend(loc="best")

    ax4 = plt.subplot(3, 2, 4)
    ax4.plot(t, ag, "k-", linewidth=2, label=r"$a_g$")
    ax4.axhline(ag_max, color="k", linestyle="--", linewidth=1,
                label=r"$a_g^{\max}$")
    ax4.axhline(ag_min, color="k", linestyle=":", linewidth=1,
                label=r"$a_g^{\min}$")
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel(r"$a_g$ [m/s$^2$]")
    ax4.set_title("Acceleration over ground")
    ax4.legend(loc="best")

    ax5 = plt.subplot(3, 2, 5)
    ax5.plot(t, np.rad2deg(omega_chi), "k-", linewidth=2,
             label=r"$\omega_\chi$")
    ax5.axhline(np.rad2deg(omega_max), color="k", linestyle="--", linewidth=1,
                label=r"$\omega_\chi^{\max}$")
    ax5.axhline(np.rad2deg(omega_min), color="k", linestyle=":", linewidth=1,
                label=r"$\omega_\chi^{\min}$")
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel(r"$\omega_\chi$ [deg/s]")
    ax5.set_title("Course rate")
    ax5.legend(loc="best")

    ax6 = plt.subplot(3, 2, 6)
    ax6.step(tu, u1, where="post", color="k", linestyle="-",
             linewidth=2, label=r"$u_1$")
    ax6.step(tu, u2, where="post", color="k", linestyle=":",
             linewidth=2, label=r"$u_2$")
    ax6.axhline(u1_max, color="k", linestyle=":", linewidth=1)
    ax6.axhline(u1_min, color="k", linestyle=":", linewidth=1)
    ax6.set_xlabel("Time [s]")
    ax6.set_ylabel("Input")
    ax6.set_title("Virtual inputs")
    ax6.legend(loc="best")

    fig.tight_layout()

    script_dir = Path(__file__).resolve().parent
    eps_file = script_dir / "optimalTrajectoryCasadi.eps"
    pdf_file = script_dir / "optimalTrajectoryCasadi.pdf"

    fig.savefig(eps_file, format="eps", bbox_inches="tight")
    fig.savefig(pdf_file, format="pdf", bbox_inches="tight")

    print(f"Saved figure: {eps_file}")
    print(f"Saved figure: {pdf_file}")

    plt.show()


if __name__ == "__main__":
    main()