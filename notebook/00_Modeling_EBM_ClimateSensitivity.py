"""
Extracted Python code from:
00_Modeling_EBM_ClimateSensitivity.ipynb

Only code cells are included. Cell boundaries are marked below.
"""

##########
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams["figure.figsize"] = (7, 4)

##########
F = 3.7       # radiative forcing [W/m2]
lam = 1.2     # climate feedback parameter [W/m2/K]
C = 10.0      # effective heat capacity [W yr/m2/K]

T = 0.0       # initial temperature anomaly [K]

print("F         =", F)
print("lambda    =", lam)
print("C         =", C)
print("Initial T =", T)

##########
dt = 0.1  # time step [yr]

dTdt = (F - lam * T) / C
T_new = T + dTdt * dt

print("Old T =", T)
print("dTdt  =", dTdt)
print("New T =", T_new)

##########
for dt_test in [0.1, 1.0, 5.0]:
    dTdt = (F - lam * T) / C
    T_new = T + dTdt * dt_test
    print(f"dt = {dt_test:4.1f} yr  T_new = {T_new:.3f} K")

##########
T = 0.0
dt = 0.1

for n in range(5):
    dTdt = (F - lam * T) / C
    T_new = T + dTdt * dt

    print(f"step {n}: old T = {T:.4f}, new T = {T_new:.4f}")

    # update
    T = T_new

##########
years = 150
dt = 0.1

time = np.arange(0, years + dt, dt)
T_series = np.zeros_like(time)

T = 0.0

for i in range(len(time)):
    T_series[i] = T
    dTdt = (F - lam * T) / C
    T_new = T + dTdt * dt
    T = T_new

plt.plot(time, T_series)
plt.xlabel("Time [yr]")
plt.ylabel("Temperature anomaly [K]")
plt.title("Zero-dimensional EBM")
plt.grid(True)
plt.show()

print("Final T =", T_series[-1])
print("Equilibrium value F/lambda =", F / lam)

##########
def run_ebm(F=3.7, lam=1.2, C=10.0, years=150, dt=0.1, T0=0.0):
    time = np.arange(0, years + dt, dt)
    T = np.zeros_like(time)
    T[0] = T0

    for n in range(len(time) - 1):
        dTdt = (F - lam * T[n]) / C
        T[n + 1] = T[n] + dTdt * dt

    return time, T

##########
time, T = run_ebm(F=3.7, lam=1.2, C=10.0, years=150)

plt.figure()
plt.plot(time, T)
plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("0D EBM: response to CO2 doubling forcing")
plt.grid(True)
plt.show()

##########
F2x = 3.7
lam = 1.2

ECS = F2x / lam

print(f"ECS = {ECS:.2f} K")

##########
F = 3.7
C = 10.0
lambda_list = [0.7, 1.0, 1.2, 1.5, 2.0]

plt.figure()

for lam in lambda_list:
    time, T = run_ebm(F=F, lam=lam, C=C, years=150)
    ecs = F / lam
    plt.plot(time, T, label=f"lambda={lam}, ECS={ecs:.1f} K")

plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("Sensitivity to climate feedback parameter")
plt.legend()
plt.grid(True)
plt.show()

##########
F = 3.7
lam = 1.2
C_list = [5, 10, 30, 100]

plt.figure()

for C in C_list:
    time, T = run_ebm(F=F, lam=lam, C=C, years=300)
    plt.plot(time, T, label=f"C={C}")

plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("Sensitivity to heat capacity")
plt.legend()
plt.grid(True)
plt.show()

# %% [Code cell 12; notebook cell 28]
def forcing_1pctCO2(time, F2x=3.7):
    return F2x * time / 70.0


def run_ebm_timevarying(F_func, lam=1.2, C=10.0, years=150, dt=0.1, T0=0.0):
    time = np.arange(0, years + dt, dt)
    T = np.zeros_like(time)
    F = F_func(time)
    T[0] = T0

    for n in range(len(time) - 1):
        dTdt = (F[n] - lam * T[n]) / C
        T[n + 1] = T[n] + dTdt * dt

    return time, T, F

# %% [Code cell 13; notebook cell 29]
time, T, F = run_ebm_timevarying(forcing_1pctCO2, lam=1.2, C=10.0, years=150)

plt.figure()
plt.plot(time, F)
plt.xlabel("Time [years]")
plt.ylabel("Radiative forcing [W m$^{-2}$]")
plt.title("Idealized 1%/yr CO2 forcing")
plt.grid(True)
plt.show()

plt.figure()
plt.plot(time, T)
plt.axvline(70, linestyle="--", label="CO2 doubling")
plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("EBM response to 1%/yr CO2 increase")
plt.legend()
plt.grid(True)
plt.show()

# %% [Code cell 14; notebook cell 31]
F2x = 3.7
lam = 1.2

ECS = F2x / lam
TCR = T[np.argmin(np.abs(time - 70))]

print(f"ECS = {ECS:.2f} K")
print(f"TCR = {TCR:.2f} K")
print(f"TCR / ECS = {TCR / ECS:.2f}")

# %% [Code cell 15; notebook cell 33]
models = pd.DataFrame({
    "model": ["Model A", "Model B", "Model C", "Model D", "Model E"],
    "lambda": [0.8, 1.0, 1.2, 1.5, 1.8],
    "C": [8, 15, 10, 30, 50]
})

models

# %% [Code cell 16; notebook cell 34]
results = []

plt.figure()

for _, row in models.iterrows():
    model = row["model"]
    lam = row["lambda"]
    C = row["C"]

    time, T, F = run_ebm_timevarying(forcing_1pctCO2, lam=lam, C=C, years=150)
    ecs = 3.7 / lam
    tcr = T[np.argmin(np.abs(time - 70))]

    results.append({
        "model": model,
        "lambda": lam,
        "C": C,
        "ECS [K]": ecs,
        "TCR [K]": tcr,
        "TCR/ECS": tcr / ecs,
    })

    plt.plot(time, T, label=model)

plt.axvline(70, linestyle="--", label="CO2 doubling")
plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("Mini CMIP: different model responses")
plt.legend()
plt.grid(True)
plt.show()

df_results = pd.DataFrame(results)
df_results

# %% [Code cell 17; notebook cell 37]
my_lambda = 1.0
my_C = 20.0
my_years = 150

time, T, F = run_ebm_timevarying(
    forcing_1pctCO2,
    lam=my_lambda,
    C=my_C,
    years=my_years
)

my_ECS = 3.7 / my_lambda
my_TCR = T[np.argmin(np.abs(time - 70))]

plt.figure()
plt.plot(time, T)
plt.axvline(70, linestyle="--", label="CO2 doubling")
plt.xlabel("Time [years]")
plt.ylabel("Global mean temperature anomaly [K]")
plt.title("Your own simple climate model")
plt.legend()
plt.grid(True)
plt.show()

print(f"lambda = {my_lambda}")
print(f"C = {my_C}")
print(f"ECS = {my_ECS:.2f} K")
print(f"TCR = {my_TCR:.2f} K")
print(f"TCR/ECS = {my_TCR / my_ECS:.2f}")
