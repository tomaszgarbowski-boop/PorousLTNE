# Model and conventions

The bed has length L, cross-sectional area A, total porosity eps and dry-solid
bulk density rho_bd. The gas-filled volume fraction is
`eg = eps - rho_bd*X/rho_l`. The cell state is
`[Tg, Ts, qv, X, qo]`, with `qv=eg*rv` and `qo=eg*rd*Y`.
Here Y is kg oxygen per kg dry carrier gas, not humid-air mass fraction.
Cell temperatures are Celsius; ideal-gas expressions use T+273.15 kelvin.

## Transport and boundaries

Dry-carrier flux is prescribed from inlet conditions: `gda=rd_in*u_in`.
Local superficial speed is `u=gda/rd`. The pressure drop is the bed integral
of `mu*u/K + rg*CF*u^2/sqrt(K)`; it does not determine the imposed flow.
The initial permeability K0 is corrected for local gas-filled porosity using
a normalized Kozeny–Carman relation. Hydrostatic pressure is excluded.

Each finite-volume balance uses `-(F_right-F_left)/dz`. Advection uses the
upstream value; nonnegative flow is required. Diffusion uses harmonic face
coefficients. Inlet gas temperature, vapor density and oxygen ratio are
Dirichlet data with a half-cell diffusive distance. Their outlet diffusive
flux is zero. Wet-solid conduction and liquid-water diffusion have zero
flux at both ends. Piecewise constant controls restart integration without
resetting states or cumulative quantities.

Water balances are `dqv/dt = -div(Fv)-Sp` and
`rho_bd*dX/dt = -div(Jw)+Sp`, with
`Fv = gda*(rv/rd)_upwind - eg*Dv*eg^nD*grad(rv)` and
`Jw=-rho_bd*Dw*grad(X)`. Oxygen transport uses
`Foxygen = gda*Y_upwind - eg*rd*Do*eg^nO*grad(Y)` with sink R.
The `oxygen=false` diagnostic mode fixes oxygen storage and evaluates
reaction using Yin; it is not a conservative oxygen model in that mode.

## Enthalpy convention

Dry gas, dry solid and liquid water have zero sensible enthalpy at 0 Celsius;
water vapor additionally stores Lv0. Per bed volume:

`Eg = (eg*rd*cda + qv*cv)*Tg + qv*Lv0`

`Es = rho_bd*(cd + cl*X)*Ts`

Gas energy includes dry-air and vapor advective enthalpy, vapor diffusive
enthalpy, and conduction. Wet-solid energy includes conduction and the
enthalpy transported by liquid redistribution. Interphase sensible transfer
is `hvol*(Tg-Ts)`. Phase mass transfer carries `Sp*(Lv0+cv*Ts)` from gas to
wet solid, with opposite signs in the two balances. Chain-rule terms from
changing vapor storage, gas porosity and liquid content are retained when
the enthalpy equations are converted to temperature rates. The phase-transfer
terms cancel exactly in the total balance. The latent contribution reported
separately is `max(Sp,0)*(Lv0+(cv-cl)*Ts)`.

Gas density follows the ideal-gas mixture at fixed total pressure.
Dry-carrier flow is quasi-steady: its accumulation, oxygen-driven changes
in carrier mass and the effects of CO2 evolution are omitted. Energy and
water conservation in this reduced formulation therefore do not establish
the validity of an unrestricted multicomponent reacting-gas model.

The default reaction combines temperature and moisture Gaussian factors with
`Y/(Ko+Y)` oxygen limitation. Reaction heat is `Ho*R`. Exchange coefficients
use standard packed-bed Nusselt/Sherwood forms with effective morphology
corrections. Water activity is one. The saturation pressure is the 1996 Buck
form reproduced by NCAR, DOI 10.26023/S2SX-K278. Default coefficients are
documented reference assumptions, not fitted material properties.

## Events and balance diagnostics

The activation criterion uses sorted cell-center Ts values with explicit
linear quantile interpolation: position `1+(N-1)*0.05`. The first rising
crossing of 23 Celsius is found by the ODE event routine. This is a continuous
surrogate for activation of 95% of the bed, not an exact volume measurement.
Alternative thresholds and quantiles are simulation options.

Seven extra states integrate external fluxes and sources using the same
time solver. The total energy residual is final minus initial stored enthalpy
minus net inlet/outlet, reaction and wall contributions. It is normalized by
the sum of the absolute integrated contributions (floor 1e-8 MJ).
The water residual is final minus initial water inventory minus net vapor
input; its normalization uses absolute net input (floor 1e-8 kg). Absolute
residuals are also returned, since percentage residuals are uninformative
near zero net exchange. These are global checks and do not replace spatial
or temporal convergence studies.

## Input contract

Reaction and exchange callbacks must be local and deterministic and return
finite N-by-1 vectors. Negative uptake or heat-transfer coefficients are
rejected. Negative airflow is rejected. The code monitors gas porosity,
water and species storage, density, and the 0–60 Celsius property range.
Crossing an invalid state stops the run with an error rather than presenting
an extrapolated trajectory as a result. Internal nonlinear trial states are
not regularized: extreme inputs may cause the stiff solver to fail before
the event is located. Such a failure requires revising the model or inputs.

## Sparse differentiation (version 0.2.0)

The production solver supplies `plt.jacobian` to `ode15s`. Five field colors
combined with three spatial colors separate the nearest-neighbour stencil.
One base and at most 15 perturbed RHS evaluations recover all field columns.
The seven diagnostic integrals are decomposed into cell-local contributions;
this recovers their global Jacobian rows without assigning one color per
cell. Diagnostic columns vanish because accumulated integrals do not feed
back into transport. The method is standard structural finite differencing,
not an analytical Jacobian or a new coloring algorithm. An independent
column-by-column finite-difference test checks the implementation.

`o.useJacobian=false` delegates differentiation to the runtime (while still
providing its sparsity pattern). The default is true. Runtime implementations
of internal differentiation can differ; recorded Octave times must not be
presented as MATLAB performance. Callbacks must remain strictly cell-local.

A 1e-30 placeholder retains a fixed stored sparse Jacobian pattern at exact zero derivatives for sparse-solver backend compatibility. It does not modify the physical RHS; round-off-scale diagnostic values can therefore appear around zero.
