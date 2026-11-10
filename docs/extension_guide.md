# Adding a material, kinetics or exchange law

## 1. Reproduce a supplied case first

From the package root in MATLAB:

```matlab
startup_porousltne
[p,s,o] = demo_case('phase_cycle',40);
out = plt.simulate(p,s,o);
plt.plotResult(out);
```

`demo_case` also supports `heating`, `oxygen_limited` and `oxygen_supplied`.
These are computed numerical demonstrations with explicitly prescribed
coefficients. They are not experimental datasets or validated material cards.
The `compost` preset retains the companion application's parameterization.

## 2. Define the material independently of operating conditions

Use `p = plt.material('wet_granulate')` or `p = plt.defaults()` as a complete
parameter structure, then replace known fields. Do not copy only a few fields
into an otherwise empty structure. The validator rejects missing and unknown
names to catch misspellings.

| Quantity | Fields | Units and interpretation |
|---|---|---|
| Geometry | `L`, `A`, `N` | Length m, area m2, integer number of uniform cells |
| Solid and pores | `rho_bd`, `eps` | Dry-solid mass per bed volume kg/m3; total porosity |
| Initial moisture | `X0`, `rho_l` | kg water/kg dry solid; liquid density kg/m3 |
| Thermal properties | `cd`, `cl`, `ks`, `kg` | Specific heats J/(kg K); effective conductivities W/(m K) |
| Pore transport | `Dv`, `Do`, `Dw` | Vapor, oxygen and liquid diffusivities m2/s |
| Flow resistance | `K0`, `CF` | Initial permeability m2; Forchheimer coefficient |
| Exchange geometry | `deq`, `chia`, `chih`, `chim` | Effective particle length m and dimensionless corrections |
| Reaction | `Rbio`, `Ho`, `Ko` | Maximum oxygen uptake kg/(m3 s), heat J/kg oxygen, oxygen ratio |
| Ambient loss | `Uw`, `AwV`, `Tamb` | W/(m2 K), wall area per bed volume 1/m, degrees Celsius |

The initial gas porosity is `eps-rho_bd*X0/rho_l` and must exceed 0.02.
Changing `X0` changes stored water, heat capacity and available gas volume;
these must not be calibrated independently as if they were unrelated inputs.
Temperature inputs use Celsius. The built-in property range is 0–60 Celsius.
Coefficients may be zero for analytical limiting cases, including zero
conductivity, species diffusivity or reaction heat. Negative coefficients
are rejected. The finite-volume equations do not permit negative airflow.

`ks` is modified by `1+alphak*(X-X0)`; use `alphak=0` for constant conductivity.
The application must keep this factor positive across the moisture range.
Preset values are illustrative; obtain measured or defensibly estimated
properties for a new application and document their provenance.

## 3. Supply the boundary program and initial state

```matlab
s = plt.schedule([0 4]*3600, [30 .10 .03; 40 .98 .03]);
o = struct('tEnd',8*3600,'target',35,'maxStep',120);
out = plt.simulate(p,s,o);
```

Each row is `[inlet_temperature_C, relative_humidity, superficial_velocity_m_s]`.
Times are seconds, begin at zero and increase strictly. For much tighter tolerances and abrupt inlet changes, set
`o.initialStep=1e-6` to resolve the initial transient of each new stage.
The full state is continuous across switches. The solver restarts at every change. It does
not interpolate controls across a step.

For nonuniform initial conditions, construct `y0 = plt.initial(p)` and edit
its five N-cell blocks, then set `o.initialState=y0`. Their order is `Tg`,
`Ts`, `qv`, `X`, `qo`. The gas storage variables are **not** relative humidity
and oxygen concentration directly: `qv=eg*rv` and `qo=eg*rd*Y`. Recompute gas
storage if changing initial temperature or water content. The last seven
entries are integrated diagnostics and should start at zero for a new run.
For continuation, retain the saved fields but reset the diagnostic entries
before starting a new independently balanced run.

## 4. Add a local reaction law

`examples/example_user_material.m` is a complete runnable example. It adds
constants inside `p.user` and assigns a named function:

```matlab
p.user.kmax=6e-5;
p.user.theta=1.06;
p.user.Koxygen=.015;
p.user.Xhalf=.3;
p.reactionFcn=@user_uptake;
```

The function has the contract `R = user_uptake(T,X,Y,p)`, where each input
field is N-by-1 and the output is nonnegative oxygen uptake in kg/(m3 bed s):

```matlab
fT=p.user.theta.^(T-25);
fX=X./(p.user.Xhalf+X);
fO=Y./(p.user.Koxygen+Y);
R=p.user.kmax*fT.*fX.*fO;
```

Use elementwise operators. Each output entry must depend only on the same
cell's input values and constants. Do not use the mean, minimum or maximum
over the whole bed inside the callback. The sparse Jacobian assumes local
constitutive coupling. Nonlocal laws need a revised sparsity/differentiation
implementation. Reaction heat is always computed as `Ho*R`, preserving the
relation between oxygen consumption and heat release. `bio=false` overrides
the callback. `oxygen=false` is a diagnostic oxygen-clamped calculation;
it is not a conservative oxygen model.

`p.user` accepts finite real numeric arrays. It is reserved for constitutive
constants and does not alter the state-vector definition. Adding a dynamic
substrate population, shrinkage or another species requires changing the
model, not placing an evolving variable in `p.user`.

## 5. Replace interphase exchange or water activity

`examples/example_user_exchange.m` demonstrates the contract
`[hvol,Sp]=f(state,p)`. Outputs are N-by-1. `hvol` is nonnegative W/(m3 K);
`Sp` is signed kg/(m3 bed s), positive for condensation and negative for
evaporation. A simple alternative is:

```matlab
hvol=p.user.hvol*ones(size(state.Ts));
rveq=p.user.aw*plt.psat(state.Ts)./(p.Rv*(state.Ts+273.15));
Sp=p.user.kmass*(state.rv-rveq);
```

The example's constant water activity is illustrative. A sorption law can
replace it using local `state.X` and `state.Ts`, with parameters supported
by measurements. The built-in saturation closure assumes freely available
water and is unsuitable for predicting complete drying or bound-water
transport. The solver stops at invalid water storage; it does not invent
a drying-rate cap. `phase=false` overrides the callback mass transfer while
retaining its sensible heat-transfer coefficient.

## 6. Verify before interpreting predictions

Run `test_porousltne(false)` after editing a law. It covers input contracts,
local balances, restart continuity, callbacks, export and the sparse
Jacobian. `run_benchmarks` exercises exact diffusion, advection and reaction
limits. `run_mesh_cost_study(3)` reproduces the reported mesh/cost comparison
and is intentionally longer. `test_porousltne(true)` additionally checks
compatibility with the retained independent compost Python fixtures.

For a new material, repeat a mesh and tolerance study, inspect local water
and oxygen values, and check energy/water residuals. A small residual means
the chosen equations are balanced; it does not validate their constitutive
assumptions. Compare physically relevant outputs with independent measurements
before using the tool for process design.

## 7. Save, reload and report a run

`plt.export(out,folder)` writes the full configuration and trajectories to
MAT plus a compact CSV. It overwrites files in the chosen output folder.
Callback identities are saved as text for portable MAT files. Restore the
original named function handle explicitly before reusing loaded parameters;
saved text is never evaluated automatically. Retain the callback source,
software version, runtime, tolerance settings and actual parameter values
with published results. No Python installation is needed for MATLAB runs.
