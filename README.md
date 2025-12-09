# Thesis code: High-Fidelity Illumination Modelling and N-Body Trajectory Optimisation for Enceladus Surface Groundtrack Design

![MATLAB](https://img.shields.io/badge/MATLAB-R2024b-orange)

## 📌 Abstract

This repository contains the MATLAB software developed for the Master's Thesis **["High-Fidelity Illumination Modelling and N-Body Trajectory Optimisation for Enceladus Surface Groundtrack Design"](LINK_TO_PDF_HERE)**, conducted within the framework of the **European Space Agency’s (ESA) GIGANTES project** under the supervision of Prof. J.P. Sanchez-Cuartielles at ISAE SUPAERO.

This codebase contributes to the broader [GIGANTES repository](LINK_TO_MAIN_REPO_HERE). The specific algorithms and tools developed for this thesis are located within the `refinements_development` directory.

### Context
Preliminary interplanetary mission design often relies on the linked-conics model, a computationally efficient but simplified approximation. However, this model faces significant limitations when transitioning to high-fidelity analysis, particularly for complex tours of Saturn's moons.

To address these challenges, this software suite implements two critical capabilities:

1.  **High-Fidelity Illumination:** A model that accurately predicts surface lighting conditions on Enceladus, incorporating both the terminator line and the frequent, long-lasting eclipses cast by an **oblate Saturn**. The model is validated against NASA’s SPICE toolkit with eclipse timing accurate to within one second.
2.  **N-Body Trajectory Optimisation:** A differential correction algorithm designed to bridge the gap between linked-conics and high-fidelity dynamics. By identifying principal perturbers (the Sun, Titan, and Saturn's $J_2$) and utilizing **State Transition Matrices (STM)** with **B-plane targeting**, the algorithm iteratively adjusts Trajectory Correction Manoeuvres (TCM) to ensure final groundtracks meet scientific objectives with minimal $\Delta V$.

## ☀️ High-Fidelity Illumination
This module computes and visualizes the illumination of groundtracks on Enceladus. It offers three visualization modes:

**1. Global Surface Projection**
For single flyby analysis, the tool can texture the entire surface of Enceladus based on the illumination conditions at the exact epoch of periapsis.

<div align="center">
  <img src="pics/gt1.png" width="100%" alt="Porkchop Plot1" />
</div>

**2. Groundtrack Coloring**
For sequences involving multiple groundtracks, it is more effective to map the illumination status directly onto the trajectory path.

<div align="center">
  <img src="pics/gt2.png" width="100%" alt="Porkchop Plot2" />
</div>

**3. Incidence Angle Validation**
Mission requirements often constrain the solar incidence angle for mapping instruments. The tool features a validation mode that highlights invalid incidence angles (non-compliant segments) in red.

<div align="center">
  <img src="pics/gt3.png" width="100%" alt="Porkchop Plot3" />
</div>

#### 🚀 Usage
To generate these visualizations, simply run the demo script:
[`refinements_development/illumination/DEMO_ILLUMINATION.m`](LINK_TO_FILE_HERE)

## 🪐 N-Body Trajectory Optimisation
This module computes and visualizes single flyby deltaV corrections on Enceladus.

<div align="center">
  <img src="pics/single_run.png" width="100%" alt="Porkchop Plot2" />
  <p><em>STM algorithm converging to the desired B-PLANE tartgeted flyby geometry</em></p>
</div>

In particular running the (link to true_anomaly_parfor.m) given:
- a flyby geometry:
    - nodein = [pump_angle_in, crank_angle_in, v_infinity_in]
    - nodeout = [pump_angle_out, crank_angle_out, v_infinity_out]
- an epoch for the flyby pericentre

The scripts loops the STM B-PLANE flyby optimization algorithm trough different value of true anomaly. This is the true anomaly before the flyby pericentre that the trajectory correction manoeuvre(TCM) is applied to match the desired flyby geometry.


<div align="center">
  <img src="pics/sweep_run.png" width="100%" alt="Porkchop Plot2" />
  <p><em>STM B-Plane optmization for different values of backward true anomaly</em></p>
</div>


To visualize better a single value of true anomaly the script pippo.m is nice.

It provides:

<div align="center">
  <img src="pics/170deg_gtrack.png" width="100%" alt="Porkchop Plot2" />
  <p><em>STM B-Plane optmization for different values of backward true anomaly</em></p>
</div>

<div align="center">
  <img src="pics/170_deg_altitude.png" width="100%" alt="Porkchop Plot2" />
  <p><em>STM B-Plane optmization for different values of backward true anomaly</em></p>
</div>

<div align="center">
  <img src="https://github.com/user-attachments/assets/2184db14-ec5c-4b90-aec2-5921ffe67305" width="100%" alt="Porkchop Plot2" />
  <p><em>STM B-Plane optmization for different values of backward true anomaly</em></p>
</div>
