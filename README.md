# 🔥 HCNG Combustion Analysis — Hydrogen-Enriched Methane in a Piston Engine

> A single, self-contained Cantera engine that quantifies what happens when you
> start mixing hydrogen into methane and burn it in a spark-ignition piston
> engine: **how fast it ignites, how hot and how hard it hits, how much power
> and efficiency you gain — and what it costs you in NOₓ.**

Course project for **Metody komputerowe w spalaniu** (Computer Methods in
Combustion) — Warsaw University of Technology, Faculty of Power and
Aeronautical Engineering.

![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Cantera](https://img.shields.io/badge/Cantera-3.x-orange)
![Mechanism](https://img.shields.io/badge/Mechanism-GRI--Mech%203.0-green)
![Status](https://img.shields.io/badge/runs%20in-~40s-success)

---

## What this does

A zero-dimensional, constant-volume, adiabatic reactor (a constant-volume
combustion bomb / RCM analogue — the textbook model for end-gas autoignition in
a piston engine) is solved with [Cantera](https://cantera.org) for fuels ranging
from **pure methane to 30 % hydrogen by mole** — and you can add any extra
hydrogen fraction you like through an interactive prompt (see
[How to run](#how-to-run)). Pure methane (0 % H₂) is the
conventional-fuel baseline that everything is compared against.

From one run you get **eight publication-quality figures**, a **CSV emissions
table**, and a printed chamber/performance report — covering four physical
stories:

| # | What's computed | How |
|---|---|---|
| ⏱️ | **Ignition delay** vs temperature, pressure, equivalence ratio | Time-integrated reactor, τ = time of max d*T*/d*t* |
| 🌡️ | **Peak temperature & pressure** vs equivalence ratio | Adiabatic constant-volume (*UV*) equilibrium |
| ⚙️ | **Power & efficiency** gain/loss vs H₂ fraction | Idealised Otto cycle, η = 1 − (1/r)^(γ−1) |
| 💨 | **Emissions** (NOₓ, CO, CO₂) vs the conventional fuel | Equilibrium composition of the burned gas |

---

## The headline result — the HCNG trade-off

Adding hydrogen is not a free lunch. It buys you cleaner carbon numbers and pays
for it in nitrogen oxides:

| At **30 % H₂** vs. pure methane | Change |
|---|---:|
| 🟢 CO₂ emissions | **−7.2 %** |
| 🟢 Ignition delay | **shorter** (faster, more reactive) |
| 🟢 Indicated efficiency | **+1.4 %** |
| 🟢 Indicated power | **+0.8 %** |
| 🔴 Thermal NOₓ | **+7.0 %** |
| 🔴 Peak temperature | **higher** (drives the NOₓ) |

> **TL;DR** — moderate hydrogen enrichment is a legitimate route to
> decarbonising an existing SI engine, but the higher flame temperature feeds
> the Zeldovich mechanism and pushes NOₓ up by roughly as much as it pulls CO₂
> down. Lean operation or after-treatment is the price of admission.

---

## Requirements

**Python 3.10+** and three packages:

```bash
pip install cantera numpy matplotlib
```

(Conda users: `conda install -c cantera cantera numpy matplotlib`.)

---

## How to run

One command. That's it.

```bash
python ignition_delay.py
```

The script is **interactive**: on start it asks how much hydrogen you want to
study.

```
 HYDROGEN FRACTION INPUT
 Default H2 fractions (always included): 0, 10, 20, 30 %
 Add extra H2 percentage(s), e.g. 15, 25, 40  (Enter to skip):
```

- The **default set (0, 10, 20, 30 %) is always included**, and 0 % (pure
  methane) is always kept as the baseline — so pressing **Enter** reproduces
  exactly the figures used in the report.
- Type any extra value(s), comma- or space-separated (`25` or `15, 25, 40`), and
  those blends are added to every figure and table. Out-of-range or non-numeric
  entries are ignored with a warning, duplicates are removed, and the list is
  sorted automatically.

In about **40 seconds** it prints the combustion-chamber definition, a
per-cycle performance table and an emissions table to the console, then writes
all figures and the CSV to the working directory.

> **Non-interactive use** — to run in a pipeline or CI, just feed an empty line:
> `echo "" | python ignition_delay.py` keeps the defaults.

### Run it in Docker (the course `canteracontainer` way)

```bash
docker build -t cantera .
docker run -it --name kinetics -v ${PWD}:/root/Simulation/ cantera
python ignition_delay.py
```

### Run it with zero install (Google Colab)

```python
!pip install cantera numpy matplotlib
!echo "" | python ignition_delay.py        # Enter = default fractions
# or feed your own:  !echo "25, 40" | python ignition_delay.py
```

---

## Outputs

| File | Figure |
|---|---|
| `fig1_delay_vs_temperature.png` | Ignition delay vs temperature (Arrhenius plot) |
| `fig2_delay_vs_pressure.png` | Ignition delay vs initial pressure |
| `fig3_delay_vs_phi.png` | Ignition delay vs equivalence ratio |
| `fig4_Tmax_vs_phi.png` | Peak temperature vs equivalence ratio |
| `fig5_Pmax_vs_phi.png` | Peak pressure vs equivalence ratio |
| `fig6_power_efficiency.png` | Power & efficiency gain/loss vs H₂ fraction |
| `fig7_emissions_change.png` | NOₓ / CO / CO₂ change vs conventional fuel |
| `fig8_NOx_vs_phi.png` | NOₓ vs equivalence ratio |
| `emissions_table.csv` | Full emissions data (ppm + emission indices) |

---

## The research chamber

All results are computed for one clearly-defined operating point, printed at the
top of every run:

| Quantity | Value |
|---|---|
| Chamber type | constant-volume, adiabatic, 0-D reactor |
| Chamber volume | 0.5 L (5 × 10⁻⁴ m³) |
| Initial temperature *T₀* | 1000 K |
| Initial pressure *p₀* | 20 atm (≈ 20.3 bar) |
| Equivalence ratio *φ* | 1.0 (reference) |
| Oxidiser | air, O₂ : N₂ = 1 : 3.76 |
| Fuel | CH₄ + H₂, *x*(H₂) = 0, 0.1, 0.2, 0.3 |
| Mechanism | GRI-Mech 3.0 (53 species, 325 reactions) |
| Idealised engine | compression ratio 10, 3000 rpm, 4-stroke |

Want a different study? Everything above lives in the `CONFIGURATION` block at
the top of `ignition_delay.py` — change the H₂ fractions, the operating point,
the chamber volume or the engine parameters and re-run.

---

## How the numbers are made

- **Ignition delay** — the reactor is marched in time and τ is taken as the
  instant of steepest temperature rise (max d*T*/d*t*), a standard and robust
  criterion.
- **Peak state** — the burned mixture is brought to chemical equilibrium at
  constant internal energy and volume (*UV*), giving the adiabatic
  constant-volume *T*ₘₐₓ and *p*ₘₐₓ.
- **Power & efficiency** — the charge energy is the constant-volume heating
  value of the mixture (the internal-energy drop when products are returned to
  *T₀*); it is run through an ideal Otto cycle, η = 1 − (1/r)^(γ−1), with γ of
  the fresh charge. Intended for **relative** fuel-to-fuel comparison.
- **Emissions** — taken from the equilibrium composition of the burned gas,
  reported both as exhaust concentrations (ppm, vol %) and as emission indices
  (g of pollutant per kg of fuel).

---

## A word on the model (read before quoting numbers)

This is an **idealised equilibrium / single-zone model**, chosen so the whole
study runs in under a minute and the physics stays transparent. Be aware that:

- Emission concentrations are **equilibrium values at the peak temperature** —
  an *upper bound*. Real exhaust NOₓ is lower because expansion freezes the
  chemistry. Treat the emission numbers as **relative trends**, not tailpipe
  figures.
- The power/efficiency estimate **ignores heat loss, friction and finite burn
  duration**, so absolute values are optimistic; the *differences* between fuels
  are the meaningful output.
- **GRI-Mech 3.0** is a natural-gas mechanism. It is an excellent fit for
  CH₄ + H₂ + air (it contains full H₂/O₂ and NOₓ sub-chemistry), which is
  exactly why this study uses methane-based blends rather than a liquid gasoline
  surrogate.

These are deliberate, defensible simplifications — and naming them is part of
the analysis.

---

## Repository

| File | What it is |
|---|---|
| `ignition_delay.py` | The complete Cantera engine — chamber report, all parameter sweeps, figures and CSV export. |
| `MKWS Report.pdf` | Full project report (Introduction, Literature Review, Model Description, Results, Conclusions). |

---

## References

- D. G. Goodwin et al., *Cantera: An object-oriented software toolkit for
  chemical kinetics, thermodynamics, and transport processes* — <https://cantera.org>
- G. P. Smith et al., *GRI-Mech 3.0* — <http://combustion.berkeley.edu/gri-mech/>
- S. R. Turns, *An Introduction to Combustion*, McGraw-Hill.
- J. B. Heywood, *Internal Combustion Engine Fundamentals*, McGraw-Hill.

---

*Built for MKWS 2026. Chemistry by Cantera, plots by matplotlib, trade-offs by
thermodynamics.*
