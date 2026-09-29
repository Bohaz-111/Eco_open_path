# Eco open paths

Fortran codes and accompanying notes for calculating the equilibrium momentum distribution of a one-dimensional anharmonic oscillator with Hamiltonian

$$
H(p,q)=\frac{p^2}{2}+\frac{q^2}{2}+\frac{q^4}{4}, \qquad m=\hbar=k_B=1.
$$

This project explores an open-path extension of the frequency-fitting approach in Zeng and Manolopoulos, [*Economised path integrals* (2026)](https://arxiv.org/abs/2607.06414), using the open-path formulation of Kapil, Cuzzocrea and Ceriotti, [*The Anisotropy of the Proton Momentum Distribution in Water* (2018)](https://arxiv.org/abs/1805.01193). Readers are assumed to be familiar with both papers and have basic experience with molecular dynamics.

## Codes

| File                   | Purpose                                                      |
| ---------------------- | ------------------------------------------------------------ |
| `open_quartic.f90`     | Standard Trotter open-path molecular dynamics with a Bussi thermostat; calculates the momentum distribution from the end-to-end estimator. |
| `Eco_open_quartic.f90` | Adds fitted internal-mode frequencies using Levenberg–Marquardt optimisation. The joint fit gives equal weight to fractional errors in harmonic mean-square end-to-end separation, including the endpoint kernel, and mean-square open-path radius of gyration. Set `mode = 0` for Trotter or `mode = 1` for Eco. |
| `Exact_quartic.f90`    | Quantum benchmark obtained by diagonalising the Hamiltonian in a harmonic-oscillator basis and constructing the thermal momentum distribution directly in momentum space. Requires LAPACK/BLAS and convergence with respect to basis size. |

## Notes

| File                                      | Contents                                                     |
| ----------------------------------------- | ------------------------------------------------------------ |
| `P1_Open_path.pdf`                        | Derives the open-path momentum estimator, including the Gaussian endpoint kernel and normalisation. |
| `P2_Cholesky_factorisation.pdf`           | Explains Cholesky factorisation in Dirac notation.           |
| `P3_Cholesky_Solve.pdf`                   | Derives the factorisation algorithm and forward/back substitution used to solve the fitting equations. |
| `P4_Levenberg_Marquardt_least_square.pdf` | Develops nonlinear least squares, Gauss–Newton and Levenberg–Marquardt, explaining the damping used in the frequency fit. |
| `P5_Quantum_benchmark.pdf`                | Derives the harmonic-basis Hamiltonian, momentum-space recurrence and thermal momentum distribution implemented in `Exact_quartic.f90`. |

Simulation and basis parameters are set near the beginning of each source file; recompile after changing them.
