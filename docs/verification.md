# Executed verification and demonstrations

Version 0.2.0, 10 November 2026.

Runtime: GNU Octave 7.3.0 on Linux x86_64, SUNDIALS/IDA, sparse Jacobian. The production source is MATLAB code. No licensed MATLAB run is claimed. The plotting helper has not been rendered in MATLAB; publication figures were generated from saved data using Python.

18 regression groups passed, plus three exact-solution benchmark families and physical-regime output assertions. Maximum scaled RHS fixture error: 1.11e-16.

## Exact solutions

Diffusion orders: 2.000738, 2.000150, 1.999893, 1.999393.
Advection L1 orders: 0.551694, 0.767813, 0.866945, 0.929463.
Normalized first-order reaction error at RelTol=1e-9: 5.20472661e-09.

All exact problems use the production solver. See benchmark_definitions.md for the closed-form solutions, cell-average convention, norms and tolerance settings.

## Coupled mesh/cost results

Reference: 160 cells with tighter tolerances and InitialStep=1e-6 s. Differences are sampled at 25 times; they are not exact coupled-model errors. Timings are medians of three sequential complete simulate calls after warm-up, excluding export/plotting.

Case | N | Maximum T difference [C] | Maximum X difference [kg/kg] | Median [s] | Min [s] | Max [s]
--- | --- | --- | --- | --- | --- | ---
heating | 10 | 2.74783 | 2.77556e-17 | 3.159 | 2.645 | 4.021
heating | 20 | 1.49741 | 0 | 3.125 | 2.904 | 3.200
heating | 40 | 0.684039 | 0 | 4.612 | 4.307 | 6.501
heating | 80 | 0.235846 | 0 | 4.467 | 3.510 | 4.485
phase_cycle | 10 | 1.23205 | 0.0173015 | 5.964 | 4.966 | 6.242
phase_cycle | 20 | 0.605875 | 0.0102393 | 6.481 | 5.733 | 7.226
phase_cycle | 40 | 0.270705 | 0.00480621 | 7.484 | 7.177 | 8.253
phase_cycle | 80 | 0.0868809 | 0.00157342 | 7.465 | 6.695 | 8.032
oxygen_limited | 10 | 1.46409 | 0.000509117 | 2.688 | 2.643 | 2.989
oxygen_limited | 20 | 0.827544 | 0.000295323 | 3.752 | 3.707 | 3.941
oxygen_limited | 40 | 0.390798 | 0.000141715 | 2.988 | 2.879 | 3.960
oxygen_limited | 80 | 0.136675 | 5.01012e-05 | 3.896 | 3.407 | 4.401

## Demonstration outcomes

* heating: N=80, final mean matrix temperature 39.999991 C, end energy residual 2.58707e-05% (global end balance), minimum oxygen ratio 0.23299986 kg/kg.
* phase_cycle: N=80, final mean matrix temperature 39.709747 C, end energy residual 4.03356e-06% (global end balance), minimum oxygen ratio 0.2220298 kg/kg.
* oxygen_limited: N=80, final mean matrix temperature 40.746891 C, end energy residual 4.62638e-07% (global end balance), minimum oxygen ratio 0.00020194976 kg/kg.
* oxygen_supplied: N=80, final mean matrix temperature 45.323610 C, end energy residual 9.10661e-08% (global end balance), minimum oxygen ratio 0.024468177 kg/kg.

Phase cycle: dry-stage water loss 3.873207 kg; humid-stage gain 2.058369 kg.

Coefficients are prescribed illustrative inputs, not experimental material data. Numerical balance and exact limiting-case verification do not validate the physical closure for a new material. The original independent Python reference remains specific to the compost benchmark.

## Reproduction

`test_porousltne(true); run_benchmarks; run_mesh_cost_study(3); run_demonstrations(80);`

The complete reports, timings and environment are in results/*.json and results/*.mat. Figure generators and CSV summaries are supplied. Staged reference runs use a smaller initial step because tight tolerances make the abrupt boundary restart more demanding.
