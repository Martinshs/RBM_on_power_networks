# Transient Stability and Random Batch Approximation of the IEEE 9-Bus System

## Overview

This repository studies transient rotor-angle dynamics in the classical three-machine representation associated with the IEEE/WSCC 9-bus benchmark. The dynamics are described by coupled swing equations, and a Random Batch Method (RBM) is used to approximate the electrical interactions between synchronous machines.

The three-machine model is obtained by reducing the electrical network to the internal generator nodes through Kron reduction. This yields a lower-dimensional generator-only ODE, while the resulting effective coupling matrix can be dense even when the original transmission network is sparse.

The objectives are to simulate the transient response after a network disturbance, study recovery or loss of synchronism, and compare the exact dynamics with a random-batch approximation.

## 1. Classical three-machine model

Let $\mathcal I=\{1,2,3\}$ denote the set of synchronous machines. For each machine $i\in\mathcal I$, we consider:

- $\delta_i(t)$: rotor electrical angle relative to a synchronous reference;
- $\omega_i(t)$: angular-frequency deviation from the synchronous reference;
- $M_i>0$: inertia coefficient;
- $D_i\geq0$: damping coefficient;
- $P_i\in\mathbb R$: effective mechanical power input;
- $K_{ij}(t)=K_{ji}(t)\geq0$: effective electrical coupling between machines $i$ and $j$, with $K_{ii}(t)=0$.

The dynamics are

$$
\dot\delta_i(t)=\omega_i(t),\qquad
M_i\dot\omega_i(t)
=
P_i-D_i\omega_i(t)
-\sum_{j\neq i}K_{ij}(t)\sin\bigl(\delta_i(t)-\delta_j(t)\bigr).
$$

In vector form,

$$
\dot\delta=\omega,\qquad
M\dot\omega=P-D\omega-F(\delta,t),
$$

where $M=\operatorname{diag}(M_1,M_2,M_3)$, $D=\operatorname{diag}(D_1,D_2,D_3)$, $P=(P_1,P_2,P_3)^\top$, and

$$
F_i(\delta,t)=\sum_{j\neq i}K_{ij}(t)\sin(\delta_i-\delta_j).
$$

The coefficients $K_{ij}$ are determined by the effective electrical network obtained after Kron reduction. In the classical lossless approximation, they depend on the reduced network susceptances and the internal voltage magnitudes of the generators.

The simulation starts from a synchronous equilibrium $(\delta^\star,0)$ satisfying

$$
P_i=\sum_{j\neq i}K_{ij}^{\mathrm{pre}}\sin(\delta_i^\star-\delta_j^\star),\qquad i\in\mathcal I.
$$

Hence $\dot\delta(0)=0$ and $\dot\omega(0)=0$. Since only angle differences enter the dynamics, the angular reference is fixed by taking $\delta_1^\star=0$.

## 2. Network disturbance

Let $0<t_f<t_c<T$ denote the disturbance and clearing times. The effective coupling matrix is modeled by

$$
K(t)=
\begin{cases}
K^{\mathrm{pre}}, & 0\leq t<t_f,\\
K^{\mathrm{fault}}, & t_f\leq t<t_c,\\
K^{\mathrm{post}}, & t_c\leq t\leq T.
\end{cases}
$$

The state remains continuous at the switching times, while the vector field changes through $K(t)$. Immediately after the disturbance,

$$
M_i\dot\omega_i(t_f^+)
=
P_i-
\sum_{j\neq i}
K_{ij}^{\mathrm{fault}}
\sin(\delta_i^\star-\delta_j^\star).
$$

This quantity is generally nonzero, so the machines accelerate or decelerate relative to the synchronous reference and a transient response develops.

The disturbance duration is $\tau=t_c-t_f$. By varying $\tau$, we can investigate whether the system recovers synchronism after the disturbance.

## 3. Transient stability and synchronization

To determine whether the machines remain synchronized, we monitor two quantities.

First, the **relative rotor angles** $\delta_i(t)-\delta_j(t)$ measure how far the machines separate in phase. Stable synchronization requires these relative angles to remain bounded rather than grow continuously.

Second, the **frequency deviations** $\omega_i(t)$ describe how fast each rotor angle moves relative to the synchronous reference. Recovery of a common frequency means that $\omega_i(t)-\omega_j(t)\to0$ for every pair of machines. Recovery of the nominal synchronous frequency additionally means $\omega_i(t)\to0$.

The main outputs of the deterministic simulation are therefore the relative rotor angles, the frequency deviations, and the maximum pairwise angular separation.

## 4. Random Batch Method

The interaction field can be decomposed over the three unordered pairs $\mathcal P=\{\{1,2\},\{1,3\},\{2,3\}\}$. For $1\leq i<j\leq3$, define

$$
\Phi_{ij}(\delta,t)
=
K_{ij}(t)\sin(\delta_i-\delta_j)(e_i-e_j),
$$

where $e_i$ is the $i$-th canonical vector of $\mathbb R^3$. The complete interaction satisfies

$$
F(\delta,t)=\Phi_{12}(\delta,t)+\Phi_{13}(\delta,t)+\Phi_{23}(\delta,t).
$$

Let $h>0$ be the RBM switching time and let $t_n=nh$. At the beginning of each interval $[t_n,t_{n+1})$, one pair $E_n\in\mathcal P$ is sampled uniformly, so that $\mathbb P(E_n=\{i,j\})=1/3$.

If $E_n=\{i,j\}$, the interaction field is approximated on $[t_n,t_{n+1})$ by

$$
\widehat F_n(\delta,t)=3\Phi_{ij}(\delta,t).
$$

The factor $3$ gives the unbiasedness property

$$
\mathbb E\!\left[\widehat F_n(\delta,t)\mid\delta\right]
=
F(\delta,t).
$$

The random-batch dynamics are

$$
\dot\delta^h(t)=\omega^h(t),\qquad
M\dot\omega^h(t)
=
P-D\omega^h(t)-\widehat F_n(\delta^h(t),t),
\qquad t\in[t_n,t_{n+1}).
$$

They use the same initial condition as the deterministic system: $\delta^h(0)=\delta^\star$ and $\omega^h(0)=0$.

Thus, on every switching interval, the complete electrical interaction is replaced by one randomly selected and appropriately rescaled pair interaction. For the three-machine system this provides a transparent proof of concept; for large reduced networks, the same mechanism replaces a dense interaction field by small random batches.

## 5. Numerical experiments

The repository considers four experiments.

### Equilibrium verification

We verify that the prescribed initial state satisfies the pre-disturbance equilibrium condition by computing

$$
R_i
=
P_i-\sum_{j\neq i}K_{ij}^{\mathrm{pre}}
\sin(\delta_i^\star-\delta_j^\star).
$$

The equilibrium is numerically consistent when $|R_i|\leq\varepsilon_{\mathrm{tol}}$ for every $i$.

### Transient response

The deterministic system is integrated through the pre-disturbance, fault, and post-disturbance intervals. We record the relative rotor angles, frequency deviations, and maximum angular separation.

### Clearing-time study

The disturbance duration $\tau=t_c-t_f$ is varied while the remaining parameters are fixed. This allows us to identify the transition between trajectories that recover synchronism and trajectories that lose synchronism.

### Random-batch approximation

The deterministic solution $(\delta,\omega)$ is compared with the random-batch solution $(\delta^h,\omega^h)$. We measure

$$
e_\delta^h(t)=\|\delta^h(t)-\delta(t)\|_2,
\qquad
e_\omega^h(t)=\|\omega^h(t)-\omega(t)\|_2.
$$

Repeated independent realizations are used to estimate mean and mean-square errors and to study their dependence on the switching time $h$.

## 6. Repository structure

```text
power-grid-rbm/
├── README.md
├── src/
│   ├── model.py
│   ├── rbm.py
│   ├── simulation.py
│   └── parameters.py
├── experiments/
│   ├── transient_stability.py
│   ├── clearing_time.py
│   └── rbm_convergence.py
├── figures/
└── requirements.txt
```

## References

1. S. Jin, L. Li, and J.-G. Liu, **Random Batch Methods (RBM) for interacting particle systems**, *Journal of Computational Physics* 400 (2020), 108877. https://doi.org/10.1016/j.jcp.2019.108877

2. U. Biccari and E. Zuazua, **A Stochastic Approach to the Synchronization of Coupled Oscillators**, *Frontiers in Energy Research* 8 (2020), 115. https://doi.org/10.3389/fenrg.2020.00115

3. D. W. M. Veldman and E. Zuazua, **A framework for randomized time-splitting in linear-quadratic optimal control**, *Numerische Mathematik* 151 (2022), 495--549. https://doi.org/10.1007/s00211-022-01290-3

4. A. Álvarez-López and M. Hernández, **Convergence, design and training of continuous-time dropout as a random batch method**, arXiv:2510.13134. https://arxiv.org/abs/2510.13134

5. M. Hernández and E. Zuazua, **Random Batch Methods for Discretized PDEs on Graphs**, arXiv:2506.11809. https://arxiv.org/abs/2506.11809

6. M. Hernández, **Random domain decomposition for parabolic PDEs on graphs**, arXiv:2508.21557. https://arxiv.org/abs/2508.21557

7. R. D. Zimmerman, C. E. Murillo-Sánchez, and R. J. Thomas, **MATPOWER: Steady-State Operations, Planning, and Analysis Tools for Power Systems Research and Education**, *IEEE Transactions on Power Systems* 26(1) (2011), 12--19.

8. J. H. Chow (ed.), **Time-Scale Modeling of Dynamic Networks with Applications to Power Systems**, Lecture Notes in Control and Information Sciences 46, Springer, 1982.