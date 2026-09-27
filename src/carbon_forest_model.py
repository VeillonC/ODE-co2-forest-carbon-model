"""
Carbon exchange model between the atmosphere, trees and soil in a forest,
written as a system of three coupled ODEs:

  CT: carbon stored in tree biomass
  CS: carbon stored in soil
  CA: carbon in the atmosphere (the compartment the forest exchanges with)

The CT equation is stiff and nonlinear (logistic-like growth capped at a
carrying capacity K), so it is advanced with an implicit RK2 (Crank-Nicolson)
step solved by Newton's method at each time step. CS and CA are then updated
with an explicit RK2 step using the newly computed CT.

The second half of the script performs a one-at-a-time sensitivity analysis:
each parameter (alpha, beta, delta, gamma, K) is swept over a range of
values while the others are kept fixed, to see how it shapes the dynamics.
"""
import numpy as np
import matplotlib.pyplot as plt

### Model functions
# The ODEs
def f_CT(t, CT):
    return - (alpha / K) * CT**2 + (alpha - beta - delta - gamma) * CT

def f_CS(t, CT, CS):
    return CT * (gamma + delta) - delta * CS

def f_CA(t, CT, CS):
    return (alpha / K) * CT**2 + CT * (-alpha + beta) + delta * CS

# F(z) for the implicit RK2 step (CT only)
def F_CT(z, yn, tn, tn1, h):
    return z - yn - (h/2) * (f_CT(tn, yn) + f_CT(tn1, z))

# Derivative of F_CT with respect to z
def F_prime_CT(z, h):
    return 1 + (h/2) * (2 * (alpha / K) * z + (alpha - beta - delta - gamma))

def newton(f, fprime, x0, eps, n_max):
    x = x0
    k = 0
    fx = f(x)
    dfx = fprime(x)
    while abs(fx) > eps and k < n_max:
        if dfx == 0:
            raise ZeroDivisionError("The derivative is zero, Newton's method cannot continue.")
        x = x - fx / dfx
        k += 1
        fx = f(x)
        dfx = fprime(x)
    return x

# Implicit RK2 for the coupled system
def solve_system_rk2(CT0, CS0, CA0, t0, tf, N):
    h = (tf - t0) / N
    t = np.linspace(t0, tf, N+1)
    CT = np.zeros(N+1)
    CS = np.zeros(N+1)
    CA = np.zeros(N+1)

    CT[0] = CT0
    CS[0] = CS0
    CA[0] = CA0

    for n in range(N):
        tn = t[n]
        tn1 = t[n+1]
        CT_n = CT[n]
        CS_n = CS[n]
        CA_n = CA[n]

        # Step 1: implicit solve for CT
        def F_z(z):
            return F_CT(z, CT_n, tn, tn1, h)

        def Fp_z(z):
            return F_prime_CT(z, h)

        CT_sol = newton(F_z, Fp_z, CT_n, 1e-10, 100)
        CT[n+1] = CT_sol

        # Step 2: explicit RK2 update for CS
        k1_CS = f_CS(tn, CT_sol, CS_n)
        k2_CS = f_CS(tn+h, CT_sol, CS_n + h * k1_CS)
        CS[n+1] = CS_n + h/2 * (k1_CS + k2_CS)

        # Step 3: explicit RK2 update for CA
        k1_CA = f_CA(tn, CT_sol, CS_n)
        k2_CA = f_CA(tn+h, CT_sol, CS_n + h * k1_CS)
        CA[n+1] = CA_n + h/2 * (k1_CA + k2_CA)

    return t, CT, CS, CA

### Project code
# Model coefficients (per year)
alpha = 0.08   # tree growth (photosynthesis)
beta = 0.02    # tree respiration -> atmosphere
delta = 0.01   # transfer from trees -> soil
gamma = 0.005  # decomposition, soil -> atmosphere
K = 2000       # maximum carbon capacity of tree biomass

# Time bounds
t0 = 0
tf = 100
# number of steps
N = 100

# Initial conditions
CA0 = 880   # atmosphere
CT0 = 500   # trees
CS0 = 1500  # soil

t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
total = CT + CS + CA

# Plot the four curves
plt.plot(t, CT, label="C_T(t)")
plt.plot(t, CS, label="C_S(t)")
plt.plot(t, CA, label="C_A(t)")
plt.plot(t, total, label="C_T + C_S + C_A")
plt.xlabel("Time")
plt.ylabel("Concentrations")
plt.grid(True)
plt.legend()
plt.title("Evolution of C_T, C_S and C_A")
plt.show()

alphas = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15]
plt.figure(figsize=(12, 7))

for a in alphas:
    alpha = a  # override alpha for this run
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CT, label=f"CT(t) with alpha={a}")

plt.title("Influence of alpha on CT(t)")
plt.xlabel("Time")
plt.ylabel("Tree carbon concentration (CT)")
plt.legend()
plt.grid(True)
plt.show()

alpha = 0.08   # reset to the default value


betas = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
plt.figure(figsize=(12, 7))

for b in betas:
    beta = b  # override beta for this run
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CT, label=f"CT(t) with beta={b}")

plt.title("Influence of beta on CT(t)")
plt.xlabel("Time")
plt.ylabel("Tree carbon concentration (CT)")
plt.legend()
plt.grid(True)
plt.show()


beta = 0.02   # reset to the default value

deltas = [0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05]
plt.figure(figsize=(12, 7))

for d in deltas:
    delta = d  # override delta for this run
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CT, label=f"CT(t) with delta={d}")

plt.title("Influence of delta on CT(t)")
plt.xlabel("Time")
plt.ylabel("Tree carbon concentration (CT)")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(12, 7))

for d in deltas:
    delta = d
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CS, label=f"CS(t) with delta={d}")

plt.title("Influence of delta on CS(t) (soil carbon)")
plt.xlabel("Time")
plt.ylabel("Soil carbon concentration (CS)")
plt.legend()
plt.grid(True)
plt.show()

delta = 0.01  # reset to the default value


gammas = [0.001, 0.0025, 0.005, 0.0075, 0.01, 0.0125, 0.015, 0.0175, 0.02]
plt.figure(figsize=(12, 7))

for g in gammas:
    gamma = g  # override gamma for this run
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CT, label=f"CT(t) with gamma={g}")

plt.title("Influence of gamma on CT(t)")
plt.xlabel("Time")
plt.ylabel("Tree carbon concentration (CT)")
plt.legend()
plt.grid(True)
plt.show()


gamma = 0.005  # reset to the default value


Ks = [1000, 1500, 2000, 2500, 3000, 3500, 4000, 20000]
plt.figure(figsize=(12, 7))

for k in Ks:
    K = k  # override K for this run
    t, CT, CS, CA = solve_system_rk2(CT0, CS0, CA0, t0, tf, N)
    plt.plot(t, CT, label=f"CT(t) with K={k}")

plt.title("Influence of K on CT(t)")
plt.xlabel("Time")
plt.ylabel("Tree carbon concentration (CT)")
plt.legend()
plt.grid(True)
plt.show()
