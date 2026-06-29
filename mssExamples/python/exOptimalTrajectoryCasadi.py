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
import os

try:
    import casadi as ca
except ImportError as exc:
    raise ImportError(
        "CasADi is required. Install it with: pip install casadi"
    ) from exc


def main() -> None:
    # ------------------------------------------------------------------
    # Discretization and weights
    # ------------------------------------------------------------------
    N = 100                         # number of control intervals
    w1 = 1.0                        # final-time weight
    w2 = 0.05                       # jerk weight
    w3 = 0.05                       # angular-acceleration weight

    # ------------------------------------------------------------------
    # Boundary conditions
    # State vector: X = [x_n, y_n, chi, U_g, a_g, omega_chi]^T
    # ------------------------------------------------------------------
    x0 = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    xf = np.array([20.0, 10.0, np.deg2rad(0.0), 1.0, 0.0, 0.0])

    # ------------------------------------------------------------------
    # Bounds
    # ------------------------------------------------------------------
    tf_min, tf_max = 0.1, 100.0
    Ug_min, Ug_max = 0.0, 3.0
    ag_min, ag_max = -0.4, 0.4
    omega_min, omega_max = np.deg2rad(-3.0), np.deg2rad(3.0)
    u1_min, u1_max = -0.5, 0.5
    u2_min, u2_max = -0.5, 0.5
    # ------------------------------------------------------------------
    # CasADi optimization problem
    # ------------------------------------------------------------------
    opti = ca.Opti()

    X = opti.variable(6, N + 1)      # x_n, y_n, chi, U_g, a_g, omega_chi
    U = opti.variable(2, N)          # u_1: jerk, u_2: angular acceleration
    tf = opti.variable()             # final time
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

    # Bounds on time, states, and virtual inputs
    opti.subject_to(opti.bounded(tf_min, tf, tf_max))
    opti.subject_to(opti.bounded(Ug_min, X[3, :], Ug_max))
    opti.subject_to(opti.bounded(ag_min, X[4, :], ag_max))
    opti.subject_to(opti.bounded(omega_min, X[5, :], omega_max))
    opti.subject_to(opti.bounded(u1_min, U[0, :], u1_max))
    opti.subject_to(opti.bounded(u2_min, U[1, :], u2_max))

    # Objective: minimum time with smoothness penalties
    J = w1 * tf + w2 * dt * ca.sumsqr(U[0, :]) + w3 * dt * ca.sumsqr(U[1, :])
    opti.minimize(J)

    # Initial guess
    s = np.linspace(0.0, 1.0, N + 1)
    X_guess = np.outer(x0, 1 - s) + np.outer(xf, s)
    X_guess[3, :] = 1.0             # initial guess for speed
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

    print(f"Optimal final time: {tf_opt:.3f} s")
    print(f"Cost: {float(sol.value(J)):.6f}")

    # ------------------------------------------------------------------
    # Plots
    # ------------------------------------------------------------------
    fig = plt.figure(figsize=(10, 8))

    ax1 = plt.subplot(3, 2, 1)
    ax1.plot(y_n, x_n, linewidth=2)
    ax1.plot(y_n[0], x_n[0], "o", label="Start")
    ax1.plot(y_n[-1], x_n[-1], "x", label="Final")
    ax1.set_xlabel(r"$y^n$ [m]")
    ax1.set_ylabel(r"$x^n$ [m]")
    ax1.set_title("Ground track")
    ax1.grid(True)
    ax1.axis("equal")
    ax1.legend(loc="best")

    ax2 = plt.subplot(3, 2, 2)
    ax2.plot(t, np.rad2deg(chi), linewidth=2)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel(r"$\chi$ [deg]")
    ax2.set_title("Course angle")
    ax2.grid(True)

    ax3 = plt.subplot(3, 2, 3)
    ax3.plot(t, Ug, linewidth=2)
    ax3.axhline(Ug_max, linestyle="--", linewidth=1)
    ax3.set_xlabel("Time [s]")
    ax3.set_ylabel(r"$U_g$ [m/s]")
    ax3.set_title("Speed over ground")
    ax3.grid(True)

    ax4 = plt.subplot(3, 2, 4)
    ax4.plot(t, ag, linewidth=2)
    ax4.axhline(ag_max, linestyle="--", linewidth=1)
    ax4.axhline(ag_min, linestyle="--", linewidth=1)
    ax4.set_xlabel("Time [s]")
    ax4.set_ylabel(r"$a_g$ [m/s$^2$]")
    ax4.set_title("Acceleration")
    ax4.grid(True)

    ax5 = plt.subplot(3, 2, 5)
    ax5.plot(t, np.rad2deg(omega_chi), linewidth=2)
    ax5.axhline(np.rad2deg(omega_max), linestyle="--", linewidth=1)
    ax5.axhline(np.rad2deg(omega_min), linestyle="--", linewidth=1)
    ax5.set_xlabel("Time [s]")
    ax5.set_ylabel(r"$\omega_\chi$ [deg/s]")
    ax5.set_title("Course rate")
    ax5.grid(True)

    ax6 = plt.subplot(3, 2, 6)
    ax6.step(tu, u1, where="post", linewidth=2, label=r"$u_1$")
    ax6.step(tu, u2, where="post", linewidth=2, label=r"$u_2$")
    ax6.set_xlabel("Time [s]")
    ax6.set_ylabel("Input")
    ax6.set_title("Virtual inputs")
    ax6.grid(True)
    ax6.legend(loc="best")

    fig.suptitle("Optimal Trajectory Generation Using CasADi")
    fig.tight_layout()

    # Save the figure
    script_dir = Path(__file__).resolve().parent

    for ax in fig.axes:
        ax.patch.set_alpha(1.0)
        for line in ax.lines:
            line.set_alpha(1.0)
        for coll in ax.collections:
            coll.set_alpha(1.0)

    pdf_file = script_dir / "optimalTrajectoryCasadi.pdf"
    fig.savefig(pdf_file, format="pdf", bbox_inches="tight")
    print(f"Saved figure: {pdf_file}")

    plt.show()

if __name__ == "__main__":
    main()
