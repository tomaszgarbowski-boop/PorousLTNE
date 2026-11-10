# Exact verification problems

Every benchmark calls the production five-field RHS and `plt.simulate`.
Parameters suppress mechanisms not part of the exact limiting problem.
There is no separate simplified numerical solver used to obtain the reported
numerical result. These are verification problems, not empirical validation.

## Closed-end heat conduction

With no phase transfer, reaction, interphase exchange, wall loss or liquid
redistribution, the wet-matrix equation is `Ts_t=alpha*Ts_zz`, where
`alpha=ks/[rho_bd*(cd+cl*X0)]`. Both ends have zero heat flux. Initially
`Ts=20+5*cos(pi*z/L)` Celsius. The exact solution multiplies the cosine by
`exp[-alpha*(pi/L)^2*t]`. The end time is `0.2*L^2/alpha`.

Numerical unknowns and exact results are compared as cell averages. The
cosine amplitude is multiplied by `sin(pi/(2N))/(pi/(2N))`. The RMS and maximum
errors are reported. Spatial orders use `log(error_N/error_2N)/log(2)`.
The expected smooth-solution order of the centered diffusion stencil is two.
Tolerances are RelTol=1e-10 and AbsTol=1e-12 to isolate spatial error.

## Passive oxygen advection

Gas temperature and humidity are constant; reaction, oxygen diffusion and
all interphase terms vanish. The equation is `Y_t+v*Y_z=0`, with
`v=u/eg` and fixed inlet ratio 0.1. A compact initial pulse adds
`0.05*sin(pi*(z-0.10)/0.60)^4` on 0.10<z<0.70 m in a 1 m bed.
It translates by 0.15 m before reaching the outlet. Exact cell averages
come from the analytic sine-fourth primitive. L1 and maximum errors are
reported. The upwind spatial stencil approaches first order under refinement.
Time tolerances are again 1e-10 and 1e-12.

## Uniform first-order oxygen uptake

With zero transport, zero reaction heat and uniform initial conditions,
`R=k*eg*rd*Y` gives `qo(t)=qo(0)*exp(-k*t)`. Here k=2e-4 1/s and the end
is 3/k. The error is the largest absolute deviation at the adaptive output
times divided by initial storage. RelTol is varied from 1e-3 to 1e-9 with
AbsTol=RelTol*1e-3. No fixed temporal convergence order is claimed for the
adaptive variable-order stiff method.

## Coupled demonstration mesh and cost study

The heating, phase-cycle and oxygen-limited cases run on 10, 20, 40 and
80 cells. The reference uses 160 cells, RelTol=1e-8 and AbsTol=1e-11;
other runs use 2e-6 and 1e-9. Reference integrations use an initial
step of 1e-6 s to resolve the sharp inlet restart under the tighter tolerance;
coarse runs retain the default 1e-3 s. Reference fields are volume-averaged onto the
coarse cells. Both adaptive trajectories are linearly interpolated to
25 common times including the endpoints. The reported temperature metrics
are maximum and RMS field differences; moisture uses the maximum field
difference. These are differences from a finer numerical mesh, not exact
physical errors or proof of an asymptotic convergence order.

The timing measurement wraps the entire `plt.simulate` call, including input
validation, integration and diagnostics, and excludes file export and
plotting. Three sequential repetitions follow a short warm-up. Median,
minimum and maximum wall-clock times are retained. Timing studies must run
without simultaneous solver jobs; machine and runtime information accompany
the report. The 160-cell reference is not included in repeated timing points.
