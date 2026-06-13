"""
=============================================================================
 Combustion analysis of hydrogen-enriched methane (HCNG) blends
 in a constant-volume research combustion chamber
 (spark-ignition piston engine context).

 Course: Metody komputerowe w spalaniu (MKWS) 2026
 Tool:   Cantera + Python
=============================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# -----------------------------------------------------------------------------
# 1. CONFIGURATION
# -----------------------------------------------------------------------------
MECH = "gri30.yaml"            # GRI-Mech 3.0: CH4, H2, O2, N2 chemistry

# Default hydrogen fractions (mole fraction of H2 in the fuel). These are ALWAYS
# included so the report figures stay reproducible; the user can add more at run
# time via the interactive prompt. 0.0 (pure methane) is always kept as baseline.
DEFAULT_H2_FRACTIONS = [0.0, 0.1, 0.2, 0.3]

# Working list actually used by the studies. It is filled at run time in main()
# from DEFAULT_H2_FRACTIONS + whatever the user types at the prompt.
H2_FRACTIONS = list(DEFAULT_H2_FRACTIONS)

# --- Research combustion chamber (constant-volume) ---------------------------
CHAMBER_VOLUME = 0.5e-3        # m^3  (= 0.5 litre)

# --- Baseline operating point ------------------------------------------------
T0_REF   = 1000.0             # K
P0_REF   = 20.0 * ct.one_atm  # Pa  (~20 bar)
PHI_REF  = 1.0                # -

# --- Idealized engine parameters (power / efficiency estimate) ---------------
COMPRESSION_RATIO = 10.0
ENGINE_RPM        = 3000.0
STROKES           = 4

# --- Integration settings ----------------------------------------------------
T_END  = 10.0
DT_MAX = 1.0e-5


# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
def fuel_composition(x_h2):
    return {"H2": x_h2, "CH4": 1.0 - x_h2}


def make_mixture(x_h2, T0, P0, phi):
    gas = ct.Solution(MECH)
    gas.set_equivalence_ratio(phi, fuel=fuel_composition(x_h2),
                              oxidizer={"O2": 1.0, "N2": 3.76})
    gas.TP = T0, P0
    return gas


# -----------------------------------------------------------------------------
# 3. IGNITION DELAY
# -----------------------------------------------------------------------------
def ignition_delay(x_h2, T0, P0, phi):
    gas = make_mixture(x_h2, T0, P0, phi)
    reactor = ct.IdealGasReactor(gas, clone=False)
    sim = ct.ReactorNet([reactor])
    sim.max_time_step = DT_MAX
    times, temps = [], []
    while sim.time < T_END:
        sim.step()
        times.append(sim.time)
        temps.append(reactor.T)
        if reactor.T > T0 + 400.0:
            break
    times, temps = np.array(times), np.array(temps)
    if len(times) < 3:
        return np.nan
    return times[np.argmax(np.gradient(temps, times))]


# -----------------------------------------------------------------------------
# 4. PEAK STATE (adiabatic constant-volume combustion)
# -----------------------------------------------------------------------------
def peak_state(x_h2, T0, P0, phi):
    gas = make_mixture(x_h2, T0, P0, phi)
    try:
        gas.equilibrate("UV", solver="vcs")   # vcs handles T > 3000 K
    except ct.CanteraError:
        gas.equilibrate("UV")                 # fallback
    return gas.T, gas.P


# -----------------------------------------------------------------------------
# 5. OTTO-CYCLE POWER & EFFICIENCY
# -----------------------------------------------------------------------------
def cycle_performance(x_h2, T0, P0, phi, volume):
    gas = make_mixture(x_h2, T0, P0, phi)
    mass = gas.density * volume
    u_react = gas.int_energy_mass
    gas.equilibrate("TV")
    u_prod = gas.int_energy_mass
    q_v = (u_react - u_prod)
    Q_in = q_v * mass

    fresh = make_mixture(x_h2, T0, P0, phi)
    gamma = fresh.cp_mass / fresh.cv_mass
    eta = 1.0 - (1.0 / COMPRESSION_RATIO) ** (gamma - 1.0)

    W_cycle = Q_in * eta
    power_events_per_s = ENGINE_RPM / 60.0 / (STROKES / 2.0)
    power = W_cycle * power_events_per_s
    return W_cycle, eta, power, Q_in


# -----------------------------------------------------------------------------
# 5b. EMISSIONS (equilibrium composition of the burned gas)
# -----------------------------------------------------------------------------
def emissions(x_h2, T0, P0, phi):
    """
    Pollutant content of the burned mixture (adiabatic, constant volume).
    Returns a dict with:
      - concentrations in the exhaust:  NO, NO2, NOx, CO  [ppm], CO2 [vol %]
      - emission indices EI            [g of pollutant per kg of fuel]
    """
    gas = make_mixture(x_h2, T0, P0, phi)
    # fuel mass fraction in the fresh charge (H2 + CH4)
    Y_fuel = float(gas["H2"].Y[0] + gas["CH4"].Y[0])

    try:
        gas.equilibrate("UV", solver="vcs")
    except ct.CanteraError:
        gas.equilibrate("UV")

    def ppm(sp):
        return float(gas[sp].X[0]) * 1e6

    def EI(sp):                       # g pollutant / kg fuel
        return float(gas[sp].Y[0]) / Y_fuel * 1000.0

    NO_ppm, NO2_ppm, CO_ppm = ppm("NO"), ppm("NO2"), ppm("CO")
    return {
        "NO_ppm":  NO_ppm,
        "NO2_ppm": NO2_ppm,
        "NOx_ppm": NO_ppm + NO2_ppm,
        "CO_ppm":  CO_ppm,
        "CO2_pct": float(gas["CO2"].X[0]) * 100.0,
        "EI_NOx":  EI("NO") + EI("NO2"),
        "EI_CO":   EI("CO"),
        "EI_CO2":  EI("CO2"),
    }


# -----------------------------------------------------------------------------
# 6. CHAMBER REPORT
# -----------------------------------------------------------------------------
def print_chamber_info():
    g = make_mixture(0.0, T0_REF, P0_REF, PHI_REF)
    n_moles = g.density_mole * CHAMBER_VOLUME
    mass = g.density * CHAMBER_VOLUME
    print("=" * 64)
    print(" RESEARCH COMBUSTION CHAMBER")
    print("=" * 64)
    print(" Type            : constant-volume, adiabatic, 0-D reactor")
    print("                   (constant-volume bomb / RCM analogue)")
    print(f" Chamber volume  : {CHAMBER_VOLUME*1e3:.3f} L  ({CHAMBER_VOLUME:.2e} m^3)")
    print(f" Initial temp.   : {T0_REF:.1f} K")
    print(f" Initial press.  : {P0_REF/1e5:.2f} bar ({P0_REF/ct.one_atm:.1f} atm)")
    print(f" Equiv. ratio phi: {PHI_REF:.2f}")
    print(f" Oxidizer        : air (O2 : N2 = 1 : 3.76)")
    print(f" Charge mass     : {mass*1e3:.3f} g  (ref. point, pure CH4)")
    print(f" Charge moles    : {n_moles*1e3:.3f} mmol")
    print()
    print(" Idealized engine model for power/efficiency:")
    print(f"   compression ratio = {COMPRESSION_RATIO:.1f}")
    print(f"   engine speed      = {ENGINE_RPM:.0f} rpm, {STROKES}-stroke")
    print("=" * 64)
    print()


# -----------------------------------------------------------------------------
# 7. STUDIES
# -----------------------------------------------------------------------------
def study_delay_vs_temperature():
    print("Fig 1: ignition delay vs initial temperature ...")
    T_range = np.linspace(900, 1400, 18)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        taus = [ignition_delay(x, T, P0_REF, PHI_REF) for T in T_range]
        plt.semilogy(1000.0 / T_range, np.array(taus) * 1e3, marker="o",
                     label=f"{int(x*100)}% H2")
    plt.xlabel("1000 / T  [1/K]"); plt.ylabel("Ignition delay  [ms]")
    plt.title(f"Ignition delay vs temperature\n(P = {P0_REF/ct.one_atm:.0f} atm, phi = {PHI_REF})")
    plt.grid(True, which="both", ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig1_delay_vs_temperature.png", dpi=150); plt.close()


def study_delay_vs_pressure():
    print("Fig 2: ignition delay vs initial pressure ...")
    P_atm = np.linspace(10, 40, 13)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        taus = [ignition_delay(x, T0_REF, P * ct.one_atm, PHI_REF) for P in P_atm]
        plt.plot(P_atm, np.array(taus) * 1e3, marker="s", label=f"{int(x*100)}% H2")
    plt.xlabel("Initial pressure  [atm]"); plt.ylabel("Ignition delay  [ms]")
    plt.title(f"Ignition delay vs pressure\n(T = {T0_REF:.0f} K, phi = {PHI_REF})")
    plt.grid(True, ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig2_delay_vs_pressure.png", dpi=150); plt.close()


def study_delay_vs_phi():
    print("Fig 3: ignition delay vs equivalence ratio ...")
    phi_range = np.linspace(0.6, 1.6, 11)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        taus = [ignition_delay(x, T0_REF, P0_REF, phi) for phi in phi_range]
        plt.plot(phi_range, np.array(taus) * 1e3, marker="^", label=f"{int(x*100)}% H2")
    plt.xlabel("Equivalence ratio  phi  [-]"); plt.ylabel("Ignition delay  [ms]")
    plt.title(f"Ignition delay vs equivalence ratio\n(T = {T0_REF:.0f} K, P = {P0_REF/ct.one_atm:.0f} atm)")
    plt.grid(True, ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig3_delay_vs_phi.png", dpi=150); plt.close()


def study_Tmax_vs_phi():
    print("Fig 4: maximum temperature vs equivalence ratio ...")
    phi_range = np.linspace(0.6, 1.6, 21)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        Tmax = [peak_state(x, T0_REF, P0_REF, phi)[0] for phi in phi_range]
        plt.plot(phi_range, Tmax, marker="o", label=f"{int(x*100)}% H2")
    plt.xlabel("Equivalence ratio  phi  [-]")
    plt.ylabel("Maximum temperature  T_max  [K]")
    plt.title(f"Peak (adiabatic) temperature vs phi\n(const. volume, T0 = {T0_REF:.0f} K, P0 = {P0_REF/ct.one_atm:.0f} atm)")
    plt.grid(True, ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig4_Tmax_vs_phi.png", dpi=150); plt.close()


def study_Pmax_vs_phi():
    print("Fig 5: maximum pressure vs equivalence ratio ...")
    phi_range = np.linspace(0.6, 1.6, 21)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        Pmax = [peak_state(x, T0_REF, P0_REF, phi)[1] / 1e5 for phi in phi_range]
        plt.plot(phi_range, Pmax, marker="s", label=f"{int(x*100)}% H2")
    plt.xlabel("Equivalence ratio  phi  [-]")
    plt.ylabel("Maximum pressure  P_max  [bar]")
    plt.title(f"Peak pressure vs phi\n(const. volume, T0 = {T0_REF:.0f} K, P0 = {P0_REF/ct.one_atm:.0f} atm)")
    plt.grid(True, ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig5_Pmax_vs_phi.png", dpi=150); plt.close()


def study_emissions_table():
    """Print and save (CSV) the emissions comparison at the reference point."""
    print("Emissions table at reference point ...")
    rows = []
    header = ("H2%", "NOx[ppm]", "CO[ppm]", "CO2[%vol]",
              "EI_NOx[g/kg]", "EI_CO[g/kg]", "EI_CO2[g/kg]")
    base = emissions(0.0, T0_REF, P0_REF, PHI_REF)
    print("\n  Exhaust emissions (equilibrium, baseline = 0% H2 = 'regular fuel'):")
    print("  {:>4} {:>9} {:>9} {:>10} {:>13} {:>12} {:>13}".format(*header))
    for x in H2_FRACTIONS:
        e = emissions(x, T0_REF, P0_REF, PHI_REF)
        rows.append((int(x*100), e))
        print("  {:>4} {:>9.1f} {:>9.1f} {:>10.2f} {:>13.2f} {:>12.2f} {:>13.1f}".format(
            int(x*100), e["NOx_ppm"], e["CO_ppm"], e["CO2_pct"],
            e["EI_NOx"], e["EI_CO"], e["EI_CO2"]))
    # percentage change vs baseline
    print("\n  Change vs regular fuel (0% H2):")
    print("  {:>4} {:>12} {:>12} {:>12}".format("H2%", "dNOx[%]", "dCO[%]", "dCO2[%]"))
    for x in H2_FRACTIONS:
        e = emissions(x, T0_REF, P0_REF, PHI_REF)
        dnox = 100.0 * (e["EI_NOx"] - base["EI_NOx"]) / base["EI_NOx"]
        dco  = 100.0 * (e["EI_CO"]  - base["EI_CO"])  / base["EI_CO"]
        dco2 = 100.0 * (e["EI_CO2"] - base["EI_CO2"]) / base["EI_CO2"]
        print("  {:>4} {:>12.1f} {:>12.1f} {:>12.1f}".format(int(x*100), dnox, dco, dco2))
    print()
    # save CSV
    with open("emissions_table.csv", "w") as f:
        f.write("H2_percent,NOx_ppm,CO_ppm,CO2_vol_percent,"
                "EI_NOx_g_per_kg,EI_CO_g_per_kg,EI_CO2_g_per_kg\n")
        for x in H2_FRACTIONS:
            e = emissions(x, T0_REF, P0_REF, PHI_REF)
            f.write(f"{int(x*100)},{e['NOx_ppm']:.2f},{e['CO_ppm']:.2f},"
                    f"{e['CO2_pct']:.3f},{e['EI_NOx']:.3f},"
                    f"{e['EI_CO']:.3f},{e['EI_CO2']:.2f}\n")


def study_emissions_change():
    """Fig 7: % change of NOx / CO / CO2 emission index vs regular fuel."""
    print("Fig 7: emissions change vs regular fuel ...")
    base = emissions(0.0, T0_REF, P0_REF, PHI_REF)
    fracs = [int(x*100) for x in H2_FRACTIONS]
    dNOx, dCO, dCO2 = [], [], []
    for x in H2_FRACTIONS:
        e = emissions(x, T0_REF, P0_REF, PHI_REF)
        dNOx.append(100.0 * (e["EI_NOx"] - base["EI_NOx"]) / base["EI_NOx"])
        dCO.append(100.0 * (e["EI_CO"] - base["EI_CO"]) / base["EI_CO"])
        dCO2.append(100.0 * (e["EI_CO2"] - base["EI_CO2"]) / base["EI_CO2"])

    xpos = np.arange(len(fracs)); width = 0.26
    fig, ax = plt.subplots(figsize=(8, 6))
    bars = [
        ax.bar(xpos - width, dNOx, width, label="NOx", color="#d62728"),
        ax.bar(xpos,         dCO,  width, label="CO",  color="#ff7f0e"),
        ax.bar(xpos + width, dCO2, width, label="CO2", color="#2ca02c"),
    ]
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(xpos); ax.set_xticklabels([f"{f}%" for f in fracs])
    ax.set_xlabel("Hydrogen fraction in fuel")
    ax.set_ylabel("Emission index change vs regular fuel  [%]")
    ax.set_title("Emissions vs regular fuel (0% H2)\n"
                 "(equilibrium products, per kg of fuel)")
    ax.grid(True, axis="y", ls=":"); ax.legend(title="Pollutant")
    for group in bars:
        for bar in group:
            h = bar.get_height()
            ax.annotate(f"{h:+.0f}", (bar.get_x()+bar.get_width()/2, h),
                        ha="center", va="bottom" if h >= 0 else "top", fontsize=8)
    plt.tight_layout(); plt.savefig("fig7_emissions_change.png", dpi=150); plt.close()


def study_NOx_vs_phi():
    """Fig 8: NOx concentration vs equivalence ratio for each blend."""
    print("Fig 8: NOx vs equivalence ratio ...")
    phi_range = np.linspace(0.6, 1.4, 17)
    plt.figure(figsize=(8, 6))
    for x in H2_FRACTIONS:
        nox = [emissions(x, T0_REF, P0_REF, phi)["NOx_ppm"] for phi in phi_range]
        plt.plot(phi_range, nox, marker="o", label=f"{int(x*100)}% H2")
    plt.xlabel("Equivalence ratio  phi  [-]")
    plt.ylabel("NOx in exhaust  [ppm]")
    plt.title(f"NOx emissions vs equivalence ratio\n"
              f"(equilibrium, T0 = {T0_REF:.0f} K, P0 = {P0_REF/ct.one_atm:.0f} atm)")
    plt.grid(True, ls=":"); plt.legend(title="H2 in fuel")
    plt.tight_layout(); plt.savefig("fig8_NOx_vs_phi.png", dpi=150); plt.close()


def study_power_efficiency():
    print("Fig 6: approximate power & efficiency gain/loss vs H2 fraction ...")
    W0, eta0, P0w, Q0 = cycle_performance(0.0, T0_REF, P0_REF, PHI_REF, CHAMBER_VOLUME)
    fractions, dPower, dEff = [], [], []
    print("\n  Per-cycle performance (idealized Otto cycle):")
    print("  H2%   W_cycle[J]   eta[%]   Power[kW]   dPower[%]  dEff[%]")
    for x in H2_FRACTIONS:
        W, eta, Pw, Q = cycle_performance(x, T0_REF, P0_REF, PHI_REF, CHAMBER_VOLUME)
        dp = 100.0 * (Pw - P0w) / P0w
        de = 100.0 * (eta - eta0) / eta0
        fractions.append(int(x * 100)); dPower.append(dp); dEff.append(de)
        print(f"  {int(x*100):>3}   {W:>9.2f}   {eta*100:>5.2f}   {Pw/1e3:>7.3f}   {dp:>8.2f}   {de:>6.2f}")
    print()
    xpos = np.arange(len(fractions)); width = 0.38
    fig, ax = plt.subplots(figsize=(8, 6))
    b1 = ax.bar(xpos - width/2, dPower, width, label="Power change", color="#1f77b4")
    b2 = ax.bar(xpos + width/2, dEff, width, label="Efficiency change", color="#d62728")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(xpos); ax.set_xticklabels([f"{f}%" for f in fractions])
    ax.set_xlabel("Hydrogen fraction in fuel")
    ax.set_ylabel("Change relative to pure CH4  [%]")
    ax.set_title("Approximate power & efficiency gain/loss\n(idealized Otto cycle, baseline = 0% H2)")
    ax.grid(True, axis="y", ls=":"); ax.legend()
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(f"{h:+.1f}", (bar.get_x() + bar.get_width()/2, h),
                        ha="center", va="bottom" if h >= 0 else "top", fontsize=9)
    plt.tight_layout(); plt.savefig("fig6_power_efficiency.png", dpi=150); plt.close()


# -----------------------------------------------------------------------------
# 8. USER INPUT (interactive prompt)
# -----------------------------------------------------------------------------
def get_user_fractions():
    """
    Ask the user for additional hydrogen percentage(s).

    The default set (0, 10, 20, 30 %) is ALWAYS included, and 0 % is always
    kept as the baseline. The user may add any number of extra values
    (comma- or space-separated). Pressing Enter keeps only the defaults, which
    reproduces exactly the figures used in the report.
    """
    defaults_pct = ", ".join(f"{int(x*100)}" for x in DEFAULT_H2_FRACTIONS)
    print("=" * 64)
    print(" HYDROGEN FRACTION INPUT")
    print("=" * 64)
    print(f" Default H2 fractions (always included): {defaults_pct} %")
    raw = input(" Add extra H2 percentage(s), e.g. 15, 25, 40  "
                "(Enter to skip): ").strip()

    extra = []
    if raw:
        for tok in raw.replace(";", ",").replace(" ", ",").split(","):
            tok = tok.strip().rstrip("%").strip()
            if not tok:
                continue
            try:
                val = float(tok)
            except ValueError:
                print(f"   - ignoring '{tok}' (not a number)")
                continue
            if not (0.0 <= val <= 100.0):
                print(f"   - ignoring {val:g} (must be between 0 and 100)")
                continue
            extra.append(round(val / 100.0, 4))

    # Combine defaults + user values + mandatory 0 %, remove duplicates, sort
    fractions = sorted(set(DEFAULT_H2_FRACTIONS) | set(extra) | {0.0})
    used = ", ".join(f"{int(round(x*100))}" for x in fractions)
    print(f" Fractions used in this run: {used} %")
    print("=" * 64)
    print()
    return fractions


# -----------------------------------------------------------------------------
# 9. MAIN
# -----------------------------------------------------------------------------
def main():
    global H2_FRACTIONS
    H2_FRACTIONS = get_user_fractions()

    print_chamber_info()
    tau = ignition_delay(0.0, T0_REF, P0_REF, PHI_REF)
    Tmx, Pmx = peak_state(0.0, T0_REF, P0_REF, PHI_REF)
    print("Sample (pure CH4 @ ref point):")
    print(f"  ignition delay = {tau*1e3:.3f} ms")
    print(f"  T_max          = {Tmx:.1f} K")
    print(f"  P_max          = {Pmx/1e5:.2f} bar\n")
    study_delay_vs_temperature()
    study_delay_vs_pressure()
    study_delay_vs_phi()
    study_Tmax_vs_phi()
    study_Pmax_vs_phi()
    study_power_efficiency()
    study_emissions_table()
    study_emissions_change()
    study_NOx_vs_phi()
    print("Done. Saved figures:")
    for f in ["fig1_delay_vs_temperature.png", "fig2_delay_vs_pressure.png",
              "fig3_delay_vs_phi.png", "fig4_Tmax_vs_phi.png",
              "fig5_Pmax_vs_phi.png", "fig6_power_efficiency.png",
              "fig7_emissions_change.png", "fig8_NOx_vs_phi.png"]:
        print("  " + f)
    print("Saved table: emissions_table.csv")


if __name__ == "__main__":
    main()
