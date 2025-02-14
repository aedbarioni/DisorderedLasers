# Stability Analysis of Disordered Lasers
Codes for stability analysis for delayed coupled lasers arrays with disorder introduced in the parameters and in the network.

See the references below for more details.

# Usage

- `StabilityAnalysisLK.m` : Analyses the stability of coupled lasers following the dynamics of the Lang-Kobayashi model (Eq. 1 of Ref. 1) as random disorder is continually introduced in one of the chosen parameters: linewidth enhancement factor (alpha), coupling strength (kappa), or frequency detuning (omega). We determine the stability at each level of disorder by calculating the Lyapunov exponent associated to the synchronous state where all lasers have same frequency, but can have different amplitudes and a small (constant) phase difference.

- `LaserDyn.m` : Solves the laser dynamics and calculates the following quantities: the standard deviation of the power spectrum 'stdFFT', the power of the combined beam 'Ptot', and the order coherence of the combined beam 'OrderCoh'. This function can be added to `StabilityAnalysisLK.m` to calculate other measures of stability such as standard deviation of the power spectrum peaks. The latter is the measure used to scale the system up to 1024 lasers.

- `Multistability.m` : Calculates all the solutions of the transcendental equations for homogeneous systems over a wide range of amplitudes and frequency shifts and determines it's stability by calculating the Lyapunov exponent. Given a chosen time delay, coupling strength, and network topology, the code plots the pairs of amplitude and frquency shift with the stability color-coded be red (unstable) and blue (stable).

- `StabilityAnalysisMulti.m` : Analyses the stability of coupled lasers as random disorder is continually introduced in a combination of the chosen parameters: linewidth enhancement factor (alpha), coupling strength (kappa), or frequency detuning (omega). To compare this case with the constrained heterogeneity, the transcendental equations remain invariant and do not need to be recalculated at each step, moreover since the frequency detuning do not affect the jacobians J1 and J2, we only need to calculate the Lyapunov exponent for the disordered alpha.

- `StabilityAnalysisLRE.m` : This code analyses the stability of coupled lasers following the dynamics of the laser rate model (Eq. S9 of Ref. 1) as parameter disorder is continually introduced in the frequency detuning (omega). We determine the stability at each level of disorder by calculating the Lyapunov exponent associated to the synchronous state where all lasers have same frequency, but can have different amplitudes and a small (constant) phase difference.


# Dependences

- `LyapExpAnalysis` : This folder contains transcendental equations for homogeneous systems (`sync.m`, `sync_hom_var.m`) used to calculate the solutions of the transcendental equations for identical parameters and for when one of the parameters varies but is kept homogeneous among the agents, respectively. The codes to solve the transcendental equations for heterogeneous systems (`sync_het.m`, `sync_het_multi.m`, `sync_het_net.m`) used when one parameter is heterogeneous, multiple parameters are heterogeneous, and when the network is disordered, respectively. This folder also includes the codes to calculate the Jacobian matrices in each case (`Jac_hom.m`, `Jac_het.m`, `Jac_het_multi.m`, `Jac_het_net.m`). Moreover we also include the functions needed to calculate the rightmost Lyapunov exponent, `dde_rightmost_eig.m` and it's dependence `stst_stabil_mod`. The latter calculation requires installation of the DDE-BIFTOOL toolbox. These codes were tested on version 3.0. See Ref. 2 and their website (https://sourceforge.net/projects/ddebiftool/) for more details.
  
- `MSFAnalysis` : This folder contains the routines to generate the MSF surface landscape (`MSFlandscape.m`) and the study of the dependence of the region of stability in the coupling strength and degree (`MSF_kappa_degree.m`, `MSF_kappa_degree_reduced.m`) for the case with identical indegrees and the reduced system, respectively.

- `NetworkOpt` : This folder contains optimization routines (`NetOpt.m`) used to minimize the Lyapunov exponent (calculated in `ObjAdj.m`) associated with the corresponding constraints. We also include the code to sparsify the network (`Sparsify.m`) by removing edges with wheight below a given cutoff and renormalizing the remaining edges.

- `LyapExpSurfaces` : This folder contains codes to calculate the Lyapunov exponent for all the directions of perturbation for a network of 3 nodes (`LyapExpLandscape.m`) connected through a decaying topology and (`LyapExpLandscapeVarTopoCH.m`) connected through different topologies and with disdorder introduce with the constrained heterogeneity strategy.

- `BasinOfAttraction` : This folder contains codes to calculate the average size of the basin of attraction and the convergence score (1 if it converges and 0 if it doesn't) for each perturbed initial condition. 

All codes were tested and run in MATLAB 2024b. To run the codes, download all files in this repository to a folder of your choice and run one of the `main` scripts of your choice. All codes generate/include the required data to run the simulations and optimization; simulations can take a few minutes on a standard laptop.

# License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

The full text of the GNU General Public License can be found in the file "LICENSE.txt".

# References
1.  AED Barioni, AN Montanari, AE Motter. Interpretable mechanisms for disorder-promoted synchronization in coupled laser networks. (2025)
2.  K Engelborghs, T Luzyanina, D Roose. Numerical bifurcation analysis of delay differential equations using DDE-BIFTOOL. *ACM Transactions on Mathematical Software*, **28**:1-21 (2002).

