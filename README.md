# Jansen_Mechanism Theo Jansen Walking Mechanism Simulation

A complete Python simulation of the **Theo Jansen walking mechanism** developed for **IE410: Introduction to Robotics**.

This project models the kinematics of the famous Theo Jansen leg mechanism, generates realistic walking trajectories, creates animated motion, and studies how changing link lengths affects gait performance.

---

# Features

- Full forward kinematics simulation
- Foot-tip trajectory generation
- Animated walking mechanism GIF
- Gait comparison with reference walking curve
- Parameter sweep analysis for link lengths
- Dark-theme scientific visualizations
- Uses Theo Jansen’s published “holy numbers”

---

# Output Preview

## 1. Foot Tip Trajectory

Shows the closed-loop trajectory traced by the foot during one crank revolution.
<img width="1523" height="362" alt="foot_trajectory" src="https://github.com/user-attachments/assets/dcf0d062-5115-4e6c-9aff-fccb688c2b39" />


---

## 2. Simulated Gait vs Reference Curve

Compares the generated walking trajectory against an ideal gait curve.

<img width="1382" height="332" alt="gait_comparison" src="https://github.com/user-attachments/assets/3f6bb000-2a14-4aa7-aabd-751adf3d1130" />


---

## 3. Effect of Varying Link `h`

Illustrates how changing the long connector length affects the gait trajectory.

<img width="1524" height="404" alt="link_variation_h" src="https://github.com/user-attachments/assets/29d88b19-c70f-4dfd-a7e5-c5a846224fec" />
---

## 4. Effect of Varying Link `m`

Illustrates how changing crank length modifies stride and step height.

<img width="1524" height="465" alt="link_variation_m" src="https://github.com/user-attachments/assets/cb1fcc7a-aa7d-4091-85c4-36f63d1083fc" />


---

## 5. Mechanism Animation

Animated simulation of the Theo Jansen mechanism.

<img width="800" height="800" alt="mechanism_animation-2" src="https://github.com/user-attachments/assets/096a0cb1-c239-4aed-81f1-a0b6849202e9" />


---

# Mechanism Overview

The Theo Jansen mechanism is a planar linkage system that converts rotational motion into a walking gait.

The simulation computes:

- Joint positions
- Foot trajectories
- Crank rotation motion
- Linkage assembly using circle-circle intersection kinematics

---

# Project Structure

```bash
theo_jansen_complete.py
│
├── Forward kinematics engine
├── Trajectory simulation
├── Animation rendering
├── Parameter sweep analysis
└── Plot generation
```

---

# Technologies Used

- Python
- NumPy
- Matplotlib
- ImageIO

---

# Installation

Clone the repository:

```bash
git clone https://github.com/your-username/theo-jansen-simulation.git
cd theo-jansen-simulation
```

Install dependencies:

```bash
pip install numpy matplotlib imageio
```

---

# Running the Simulation

Run the main script:

```bash
python theo_jansen_complete.py
```

---

# Generated Outputs

The program automatically generates:

| File | Description |
|---|---|
| `foot_trajectory.png` | Closed foot trajectory |
| `mechanism_animation.gif` | Animated walking mechanism |
| `gait_comparison.png` | Comparison with ideal gait |
| `link_variation_m.png` | Effect of crank length |
| `link_variation_h.png` | Effect of connector length |

---

# Theo Jansen Holy Numbers

The simulation uses the original linkage ratios proposed by Theo Jansen:

```python
DEFAULT_LINKS = {
    'a': 38.0,
    'b': 41.5,
    'c': 39.3,
    'd': 40.1,
    'e': 55.8,
    'f': 39.4,
    'g': 36.7,
    'h': 65.7,
    'i': 49.0,
    'j': 50.0,
    'k': 61.9,
    'l': 7.8,
    'm': 15.0,
}
```

---

# Kinematic Method

The mechanism is solved using:

- Circle-circle intersection geometry
- Planar linkage constraints
- Forward kinematic chain computation

The foot point is computed for each crank angle over a full revolution:

:contentReference[oaicite:0]{index=0}

---

# Results

The generated trajectory produces:

- Smooth walking motion
- Large stride length
- Stable stance phase
- Realistic gait curve

The parameter sweep also demonstrates how small geometric changes significantly affect locomotion behavior.

---


# Reference

- Theo Jansen Mechanisms
- Jansen, T. (2007)
- Shin et al. (2018), Journal of Mechanical Robotics
- Single-DOF Gait Trainer studies

---


