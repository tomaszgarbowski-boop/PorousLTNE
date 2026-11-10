# PorousLTNE 0.2.0

A small MATLAB package for coupled heat, moisture and oxygen transport in
ventilated wet porous beds, with separate gas and wet-matrix temperatures.
Authors of the accompanying manuscript: Tomasz Garbowski and Dawid Karasiewicz.
Correspondence: tomasz.garbowski@up.poznan.pl.

## Start

Unzip, open the `PorousLTNE` folder in MATLAB, then run:

```matlab
startup_porousltne
[p,s,o] = demo_case('heating',40);
out = plt.simulate(p,s,o);
plt.plotResult(out);
plt.export(out,fullfile('results','my_heating_run'));
```

The production code uses base MATLAB and `ode15s`, without external toolboxes,
MEX compilation, Python or network services. MATLAB R2020b or later is the
intended target. The recorded execution used GNU Octave 7.3.0 on Linux with
SUNDIALS/IDA and sparse Jacobian support; actual MATLAB execution remains
required before submission. See `docs/verification.md` for measured results.

## Three distinct demonstrations

* `demo_case('heating',N)`: warm-air heating of a nonreactive moist granulate;
  phase transfer is disabled to isolate sensible heating.
* `demo_case('phase_cycle',N)`: a dry-air stage followed by humid warm air;
  evaporation and condensation use the same conservative phase-transfer law.
* `demo_case('oxygen_limited',N)`: reactive self-heating with transported oxygen
  and low supply. `oxygen_supplied` increases flow for comparison.

`run_demonstrations(80)` runs and exports all four configurations. The presets
are explicit numerical scenarios, not measured material databases. The
original compost application remains available as `example_compost`.

## Verification and reproduction

```matlab
test_porousltne(false)    % contracts, invariants, restart, Jacobian and export
test_porousltne(true)     % also independent Python compost fixtures
run_benchmarks            % exact diffusion, advection and reaction limits
run_mesh_cost_study(3)    % repeated timings and mesh comparisons
example_user_material    % local custom kinetics via p.user
example_user_exchange    % custom water activity and phase exchange
```

`tests/fixtures` and `verification/reference_python.py` preserve an independent
compost implementation. Python is only needed to regenerate those fixtures
or publication figures. Install `verification/requirements.txt`, then run
`python verification/make_figures.py` to regenerate figures from saved results.
The scripts write to `results` and overwrite outputs in the same named folders.

## Interfaces and documentation

The package exposes `plt.defaults`, `plt.material`, `plt.schedule`,
`plt.simulate`, `plt.plotResult`, `plt.export` and `plt.sweep`. Local
reaction and exchange callbacks alter constitutive laws independently of
transport. A 15-color sparse finite-difference Jacobian supports the stiff
integration; it assumes callbacks are cell-local.

* `docs/extension_guide.md`: material parameters, units, initial states,
  custom kinetics, custom exchange, verification and export.
* `docs/model.md`: conservation equations, enthalpy convention and scope.
* `docs/benchmark_definitions.md`: exact solutions, error norms and timing protocol.
* `docs/verification.md`: actual runtime, results and limitations.
* `docs/release.md`: author and public-release fields remaining before submission.

Times are seconds, temperatures Celsius, humidity a fraction, and other
parameters SI. Outputs retain all five fields, seven cumulative diagnostics,
parameters, schedule, runtime and event state. Callback handles are saved as
text in portable MAT files; restore named handles explicitly before reusing
loaded parameters. The exporter never evaluates saved text.

## Limits and release status

The model assumes fixed total porosity, a uniform one-dimensional grid,
nonnegative prescribed flow and built-in properties over 0–60 Celsius.
The default phase law has water activity one. It does not represent bound
water, complete drying, shrinkage, carbon dioxide or substrate depletion.
Dry-carrier mass is prescribed and not reduced by oxygen uptake. Reduced
model balance checks do not establish experimental validity. The conditioning
energy indicator is idealized, not measured plant electricity.

Version 0.2.0 is prepared for author review. No public repository or permanent
archive identifier has been assigned. `Licence.txt` contains proposed MIT
wording requiring approval and identification of the copyright holder.
