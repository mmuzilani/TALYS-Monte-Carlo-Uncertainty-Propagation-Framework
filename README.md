# Nuclear Data Preparation

## Introduction to Nuclear Reaction Data

In stellar environments, heavy elements are built through a sequence of neutron-capture reactions and β-decays.  
To model such nucleosynthesis, we require reliable nuclear physics input.  
The foundation of all calculations is the neutron-capture reaction:

$$
(A,Z) + n \rightarrow (A+1,Z) + \gamma ,
$$

characterized by an energy-dependent probability called the **cross section**, $\sigma(E)$.

Since physical conditions inside stars vary with temperature, neutron energies follow statistical distributions.  
Therefore, raw $\sigma(E)$ must be transformed into astrophysical quantities such as **MACS** and **reaction rates**.

We focuses on preparing these essential nuclear inputs:

- Energy-dependent cross section σ(E)  
- Maxwellian-Averaged Cross Section (MACS)  
- β-decay half-life data  

---

## Neutron-Capture Cross Section σ(E)

The cross section $\sigma(E)$ measures the probability that a nucleus will capture a neutron of energy **E**.

- **Unit:** barn (1 b = 10⁻²⁸ m²)  
- **Behavior:** varies strongly with energy and shows resonance peaks  
- **Importance:** determines how fast isotopes transform during nucleosynthesis  

For each Fe isotope in the chain (56→57→58→59→60), TALYS provides:

```
Energy(MeV), CrossSection(b)
```

This σ(E) dataset is the **most fundamental input**, because *all other astrophysical quantities are derived from it*.

---

##  Neutron Energy in Stars and Maxwell–Boltzmann Distribution

Inside stars, neutrons are in **thermal equilibrium** with the plasma.  
They do **not** have one fixed energy. Instead, they follow the **Maxwell–Boltzmann distribution**:

$$
f(E) \propto E^{1/2} \, e^{-E/kT}.
$$

Where:

- **kT** = thermal energy (in keV)  
- **T** = stellar temperature  
- **E** = neutron kinetic energy  

Because σ(E) depends on E, the *effective reaction probability* requires weighting σ(E) by this distribution.  
Thus, a simple σ(E) curve is not enough — we must average it over the Maxwellian.

---

## Maxwellian-Averaged Cross Section (MACS)

The **MACS** is the average of σ(E) over a Maxwell–Boltzmann distribution:

$$
MACS(kT)=
\frac{2}{\sqrt{\pi} (kT)^2}
\int_{0}^{\infty}
\sigma(E)\, E\, e^{-E/kT}\, dE .
$$

### Meaning of MACS

- It is the **astrophysical cross section**, relevant for stellar interiors  
- It smooths out resonance peaks of σ(E)  
- It depends on temperature  

**Unit:** barn  
MACS is essential because stars are thermal environments — σ(E) alone cannot describe reaction rates.

---

## Experimental MACS (KADoNiS / EXFOR)

TALYS gives **theoretical** σ(E). But stars follow **real nuclear physics**, so we need experimental MACS from:

- **KADoNiS** (s-process database)  
- **EXFOR** (global reaction data library)

I prepare **four experimental MACS files**, one for each neutron-capture step:

- 56Fe(n,γ)57Fe  
- 57Fe(n,γ)58Fe  
- 58Fe(n,γ)59Fe  
- 59Fe(n,γ)60Fe  

These experimental values allow **validation of TALYS output**.

---

## Radioactive Decay Data and Decay Constants

Not all isotopes in the chain are stable:

- **Stable:** 56Fe, 57Fe, 58Fe  
- **Radioactive:** 59Fe, 60Fe  

Radioactive isotopes decay according to:

$$
N(t) = N_0 \, e^{-\lambda t},
$$

where the **decay constant** λ is:

$$
\lambda_\beta = \frac{\ln 2}{T_{1/2}}.
$$

# From MACS to Reaction Rate ⟨σv⟩


In a star, particles move rapidly due to thermal motion. Thus, the true interaction strength must consider:

- cross section σ(E)  
- neutron velocity v(E)  
- Maxwell–Boltzmann energy distribution  

The **reaction rate per particle pair** is:

$$
\langle \sigma v \rangle
= \int_0^\infty \sigma(E) \, v(E) \, f(E,T) \, dE .
$$

This tells us how frequently neutron captures actually happen in a stellar plasma.

---

## Relation between MACS and ⟨σv⟩

For neutron captures, we use a special identity:

$$
\langle \sigma v \rangle = v_T \, MACS(kT),
$$

where

$$
v_T = \sqrt{\frac{2kT}{m_n}}
$
is the thermal neutron speed.

Thus:
- **MACS** captures the energy dependence (σ(E) and Maxwellian weight)
- **v_T** adds the thermal velocity scale

Together, they produce a proper stellar reaction rate.

**Units of ⟨σv⟩:**  
m³ s⁻¹

---

# From Reaction Rate to Capture Rate λ

Stars contain an enormous number of free neutrons.  
If neutron density is \( n_n \), the **capture rate** per nucleus is:

$$
\lambda_{\text{cap}} = n_n \, \langle \sigma v \rangle .
$$

This λ is the **probability per second** that a nucleus captures a neutron.

Interpretation:

- λ = 1 s⁻¹ → reaction happens once per second  
- λ = 10⁻⁶ s⁻¹ → reaction takes ≈ 10⁶ s to occur  

In my project, λ is computed for four Fe captures:

- 56Fe → 57Fe  
- 57Fe → 58Fe  
- 58Fe → 59Fe  
- 59Fe → 60Fe  

These λ values form the backbone of the reaction network.

---

# β-Decay Rates from Half-Lives

Some isotopes are unstable:

- 59Fe → β⁻ decay  
- 60Fe → β⁻ decay  

A radioactive nucleus decays according to:

$$
N(t) = N_0 e^{-\lambda_\beta t},
$$

where the decay constant is:

$$
\lambda_\beta = \frac{\ln 2}{T_{1/2}}.
$$

Meaning:

- Short \(T_{1/2}\) → fast decay (large λβ)  
- Long \(T_{1/2}\) → slow decay (small λβ)  

These λβ values are essential in the final network equations.

---

# Constructing the 5-Isotope Nuclear Reaction Network

The isotopes form a sequential chain:

$$
^{56}Fe \rightarrow {}^{57}Fe \rightarrow {}^{58}Fe \rightarrow {}^{59}Fe \rightarrow {}^{60}Fe.
$$

## General Reaction Network Equation

## ⭐ General Reaction Network Equation

The evolution of isotope abundances in a neutron-capture chain is governed by:

$$
\frac{dY_i}{dt}
= -\lambda_i\,Y_i 
+ \lambda_{i-1}\,Y_{i-1}
- \lambda_{\beta,i}\,Y_i
$$

### Meaning of the terms
- Isotope destroyed by neutron capture:  
  $$ -\lambda_i Y_i $$
- Isotope created from the previous isotope:  
  $$ +\lambda_{i-1} Y_{i-1} $$
- Isotope destroyed by β-decay:  
  $$ -\lambda_{\beta,i} Y_i $$

---

# ⭐ Explicit Equations for the 5-Isotope Fe Reaction Network

### **1. \(^{56}\mathrm{Fe}\)** — stable
$$
\frac{dY_{56}}{dt}
= -\lambda_{56}\,Y_{56}
$$

---

### **2. \(^{57}\mathrm{Fe}\)**
$$
\frac{dY_{57}}{dt}
= \lambda_{56}\,Y_{56}
- \lambda_{57}\,Y_{57}
$$

---

### **3. \(^{58}\mathrm{Fe}\)**
$$
\frac{dY_{58}}{dt}
= \lambda_{57}\,Y_{57}
- \lambda_{58}\,Y_{58}
$$

---

### **4. \(^{59}\mathrm{Fe}\)** — radioactive
$$
\frac{dY_{59}}{dt}
= \lambda_{58}\,Y_{58}
- \left( \lambda_{59} + \lambda_{\beta,59} \right)\,Y_{59}
$$

---

### **5. \(^{60}\mathrm{Fe}\)** — radioactive
$$
\frac{dY_{60}}{dt}
= \lambda_{59}\,Y_{59}
- \lambda_{\beta,60}\,Y_{60}
$$

---

# ⭐ Interpretation of the Abundance Curves

- \(Y_{56}(t)\) → decreases steadily (consumed by neutron capture)
- \(Y_{57}(t)\) → rises at first, then slowly falls  
- \(Y_{58}(t)\) → rises later, following buildup of 57Fe  
- \(Y_{59}(t)\) → rises but then declines due to β-decay  
- \(Y_{60}(t)\) → final product; gradually accumulates  

These curves show how material flows through the Fe isotopes during nucleosynthesis.


# Matrix Form of the Network (for Fast Numerical Solving)

We can rewrite the system in **matrix form**:

$$
\frac{d\mathbf{Y}}{dt} = A \mathbf{Y},
$$

where **A** is the 5×5 rate matrix containing λ values.

Example structure:

$$
A =
\begin{pmatrix}
-\lambda_0 & 0 & 0 & 0 & 0 \\
\lambda_0  & -\lambda_1 & 0 & 0 & 0 \\
0          & \lambda_1  & -\lambda_2 & 0 & 0 \\
0          & 0          & \lambda_2  & -(\lambda_3 + \lambda_{\beta,3}) & 0 \\
0          & 0          & 0          & \lambda_3                        & -\lambda_{\beta,4}
\end{pmatrix}
$$

The solution is:

$$
Y(t) = e^{At} \, Y(0).
$$

SciPy computes this using:

- Schur decomposition  
- Matrix exponential expm()  
- Smooth time sampling  

This yields the abundance curves.

---

# Meaning of the Abundance Plots

The numerical solution gives abundance curves:

- \(Y_{56}(t)\) decreases (converted into heavier isotopes)  
- \(Y_{57}(t)\) rises, then falls  
- \(Y_{58}(t)\) rises later  
- \(Y_{59}(t)\) rises but decays due to β-decay  
- \(Y_{60}(t)\) accumulates last (final product)  

These curves represent **real astrophysical nucleosynthesis** under neutron irradiation.

---

## Cross Section Visualization — σ(E) vs Energy

The neutron-capture cross section σ(E) shows the probability of capture as a function of neutron energy.  
Plotting σ(E):

- Reveals **resonances**  
- Shows energy structure of the reaction  
- Helps identify errors in TALYS data  
- Forms the base for all later calculations  

### Interpretation of the Plot

- **X-axis:** neutron energy \(E\) (MeV)  
- **Y-axis:** cross section \(\sigma(E)\) (barn)  
- Sharp peaks → resonance regions  
- Smooth/flat parts → non-resonant behavior  

Correct σ(E) structure ensures accurate MACS and reaction rates.

---

## Maxwellian Integrand Plot — Checking the MACS Calculation

To compute MACS, the weighted integrand is:

$$
I(E) = \sigma(E)\,E\, e^{-E/kT}.
$$

Plotting this function shows:

- Which neutron energies contribute most to the MACS  
- Whether your interpolation and integration are stable  
- Whether the σ(E) dataset is realistic  

### Interpretation of the Plot

- **X-axis:** neutron energy (keV)  
- **Y-axis:** integrand value \( I(E) \) (arbitrary units)  
- Peak location → “effective energy window”  
- Width → energy spread contributing to the MACS  

For s-process temperatures (kT = 30 keV), the peak usually lies around 20–60 keV.

This plot is a **scientific validation** of MACS.

---

## MACS Comparison — TALYS vs Experimental

To assess accuracy, we compare:

- **TALYS MACS** (theory)  
- **KADoNiS/EXFOR MACS** (experiment)  

These values are plotted together as:

- Points for TALYS  
- Points for experiment  
- Text labels indicating percentage difference  

## Reaction Rate vs Temperature — ⟨σv⟩(T)

Reaction rate increases with temperature due to increased neutron velocities:

$$
\langle \sigma v \rangle = v_T \times MACS(kT),
$$

so temperature strongly shapes nucleosynthesis.

### Interpretation of the Plot

- **X-axis:** temperature (10⁷–10⁹ K)  
- **Y-axis:** reaction rate ⟨σv⟩ (m³/s)  
- Smooth rising curve  
- At low T → slow neutrons → slow capture  
- At high T → fast neutrons → rapidly increasing capture rate  
- Shape determines reaction competition and branching  

This plot connects laboratory σ(E) to stellar behavior.

---

## Abundance Evolution — Solution of the 5-Isotope Network

Using the ODE system:

$$
\frac{dY_i}{dt}
=
-\lambda_i Y_i
+
\lambda_{i-1}Y_{i-1}
-
\lambda_{\beta,i}Y_i,
$$

we compute the time evolution of all five Fe isotopes.

### Interpretation of Abundance Curves

- **56Fe** → initial seed, decreases  
- **57Fe** → rises then falls  
- **58Fe** → rises at later times  
- **59Fe** → rises but eventually decays  
- **60Fe** → final product, accumulates  

These curves represent **real nucleosynthesis flow** under neutron exposure.

---

## Monte-Carlo Uncertainty Bands — (16th–84th Percentiles)

Nuclear data has uncertainties (MACS measurement error, model error, etc.).  
To quantify their effect, we perform Monte-Carlo sampling:

$$
MACS_{\text{random}} =
MACS_{\text{mean}} \cdot (1 + \epsilon),
$$

where:

- \(\epsilon\) ~ Gaussian(0, σ)  

Repeating the network simulation many times gives many abundance curves. From all Monte-Carlo solutions:

- **50th percentile** → median curve  
- **16th percentile** → lower uncertainty limit  
- **84th percentile** → upper uncertainty limit  

- **Narrow band:** prediction is robust  
- **Wide band:** uncertainties strongly affect nucleosynthesis  
- **Shifted median:** TALYS may deviate from experiment  

Uncertainty quantification is essential in modern nuclear astrophysics.

---
