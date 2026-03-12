# 2D Fluid Solver with Obstacles and Particles

**Author:** Ilya Tsivilskiy

A web browser-based, interactive 2D incompressible Navier–Stokes solver with immersed obstacles and Lagrangian tracer particles. The fluid is simulated on a Cartesian grid with free-slip boundary conditions on solid boundaries; optional tracer particles follow the flow with drag and gravity and collide with obstacles.

---

## 1. What it does

- **Fluid:** Solves the incompressible 2D Navier–Stokes equations in a square domain with constant kinematic viscosity. Velocity and pressure are defined at the same grid nodes (collocated layout).
- **Obstacles:** Solid regions are represented by a scalar mask. On the fluid–solid interface, free-slip is enforced: the normal velocity component is zero and the tangential component is preserved (velocity is projected onto the local tangent).
- **Boundaries:** Each side of the domain (top, bottom, left, right) can be set to a prescribed velocity (Dirichlet) or zero gradient (Neumann) for each velocity component.
- **Particles:** Optional point particles are advected by the flow with a drag law and gravity; they reflect off solid surfaces with optional velocity damping.
- **Interaction:** You can draw or erase obstacles, inject velocity with a brush, add particle clusters, and reset the simulation or only the obstacles. The velocity colormap scale is updated periodically from the current flow field.

### Possible applications

The app is primarily a **toy model for educational purposes**, but it can illustrate and cover such effects as:

- **Lid-driven cavity** — set a moving lid (e.g. top or bottom) via boundary conditions and observe the recirculation and secondary vortices.
- **Vortex shedding** — place an obstacle (e.g. a blob or cylinder-like shape) in the path of a uniform or shear flow (set via BCs or the velocity brush) and watch the von Kármán street.
- **Flow in a channel** — use obstacles to form walls and set inlet/outlet-like conditions on the sides to mimic channel flow.
- **Air wind tunnel** — prescribe inflow on one side and zero gradient or outflow on the opposite side; add obstacles to represent models and use particles to visualize streamlines or dust/droplets with a physically consistent drag law.
- **Particles with correct drag law** — tracer particles use the Morsi–Alexander *C*<sub>*d*</sub>(Re) correlation, so they respond realistically to the local flow and Reynolds number (Stokes and beyond), suitable for qualitative studies of particle transport, settling, or dispersion in 2D flows.

---

## 2. How to use

### Opening / launching

- **From disk:** Open `index.html` in a modern browser (e.g. Chrome, Firefox, Edge) via `file://` or by double-clicking the file.
- **With a local server (optional):** Serve the project root with any static server (e.g. `npx serve .`, `python -m http.server`) and open the given URL. No build step is required; the app is plain HTML, CSS, and JavaScript.

### Toolbar

| Button | Action |
|--------|--------|
| **Draw obstacles** | Click or drag on the canvas to add solid obstacles (circles). |
| **Erase** | Click or drag to remove obstacles. |
| **Inject velocity** | Click to add downward velocity; drag to add velocity in the drag direction. Affects fluid cells in a circular brush. |
| **Inject particles** | Click or drag to add a cluster of tracer particles at the pointer location. |
| **Reset** | Clears fluid velocities and pressure, removes all particles, clears all obstacles, and resets boundary condition inputs to zero. |
| **Reset obstacles** | Clears only the obstacle mask; fluid and particles are unchanged. |

After drawing or erasing obstacles, the solid mask is smoothed (one pass of a 3×3 filter) and binarized so the interface is well-defined for normals and particle collision.

### Boundary conditions

The four panels around the canvas (Top, Bottom, Left, Right) set boundary conditions:

- **Vx, Vy:** Prescribed velocity components in m/s. Values are clamped to a CFL-safe maximum and displayed as two decimal places (e.g. `0.00`, `-0.05`).
- **Zero ∇:** If checked, that side uses zero gradient (Neumann) for that component instead of the prescribed value; the inputs are disabled.

The solver applies these each step when updating velocity boundary values.

### Display

- **Background:** Dark blue-gray.
- **Fluid cells:** Colored by velocity magnitude (blue = low, red = high). The scale is recomputed every 10 fluid time steps from the current maximum velocity in the domain.
- **Solid cells:** Dark gray.
- **Velocity arrows:** Subsampled arrows show flow direction; arrowheads indicate direction.
- **Particles:** White dots.

---

## 3. How it works

### Governing equations

**Fluid (incompressible 2D Navier–Stokes):**

Momentum:

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\nabla p + \nu \nabla^2 \mathbf{u}$$

Continuity (incompressibility):

$$\nabla \cdot \mathbf{u} = 0$$

Here **u** = (*u*, *v*) is the velocity, *p* is pressure (divided by constant density), and *ν* is the kinematic viscosity.

**Solid boundaries:** Free-slip: zero normal velocity and no shear stress. So on the solid surface, ∂*p*/∂*n* = 0 and the velocity is projected onto the tangent direction (no penetration, no tangential friction).

**Tracer particles:** Each particle obeys

$$\frac{\mathrm{d}\mathbf{u}_p}{\mathrm{d}t} = \frac{\mathbf{u} - \mathbf{u}_p}{\tau_p} + \mathbf{g}\,\frac{\rho_p - \rho}{\rho_p} + \mathbf{a},$$

where **u**<sub>*p*</sub> is particle velocity, **u** is fluid velocity at the particle position, τ<sub>*p*</sub> is the particle relaxation time, **g** is gravity, ρ<sub>*p*</sub> and ρ are particle and fluid density, and **a** is any extra acceleration (here zero). The relaxation time is

$$\tau_p = \frac{24}{18}\,\frac{\rho_p d_p^2}{\mu C_d \mathrm{Re}},$$

with particle diameter *d*<sub>*p*</sub>, dynamic viscosity *μ*, drag coefficient *C*<sub>*d*</sub>(Re), and particle Reynolds number Re = |**u** − **u**<sub>*p*</sub>| *d*<sub>*p*</sub> / *ν*. *C*<sub>*d*</sub>(Re) is given by the Morsi–Alexander (1972) correlation for a sphere over a range of Re (Stokes regime and beyond). When a particle hits the solid, its velocity is reflected about the surface normal and multiplied by (1 − damp) to model inelastic collision.

---

### What is accounted for

- **Advection:** Semi-Lagrangian (backtrace and interpolate) so the scheme is stable for large CFL.
- **Pressure:** Pressure Poisson equation with a liquid-only Laplacian (solid cells do not contribute), so that ∇·**u** = 0 is enforced and the pressure gradient is applied only in the fluid. Boundary conditions are zero gradient on top/bottom and zero (or copy) on left/right for pressure.
- **Diffusion:** Implicit (backward Euler in time, Laplacian in space), solved with Gauss–Seidel iteration.
- **Vorticity confinement:** An extra force proportional to the gradient of vorticity is added to counteract numerical dissipation of small-scale vorticity (optional strength parameter).
- **Obstacles:** Solid mask is smoothed (e.g. 3×3 average) and binarized so that surface normals (from the gradient of the mask) are well-defined. Velocity in solid is set to zero; at fluid nodes adjacent to solid, velocity is projected onto the local tangent to enforce free-slip.
- **Particles:** Fluid velocity and solid/normal at the particle position are obtained by bilinear interpolation from the grid. Velocity is advanced with the analytical solution of d*u*/d*t* = *A* *u* + *B* over a time step; position is advanced with the updated velocity. Collision: if the new position is inside the solid, velocity is reflected and damped and position is re-advanced.

---

### Numerical methods and solver structure

The scheme follows a nodal, finite-difference style implementation on a 2D Cartesian grid (collocated *u*, *v*, *p*, and solid mask at the same nodes). The time step is fixed.

**Spatial discretization:**

- Central differences for ∂/∂*x*, ∂/∂*y*, so divergence, gradient, and curl are second-order in the interior.
- Scalar fields are sampled at arbitrary physical coordinates via bilinear interpolation from the four surrounding nodes.
- Boundary conditions are applied by overwriting boundary values (copy from interior for Neumann, prescribed for Dirichlet).

**Time step (fluid), per step:**

1. **Vorticity confinement:** Compute vorticity ω = ∂*v*/∂*x* − ∂*u*/∂*y*, then its gradient; add a force along the gradient of ω (vorticity confinement). Update **u** with this force and apply free-slip (solid).
2. **Advection:** For each component, semi-Lagrangian step: for each grid point, backtrace **x** − **u**(**x**) Δ*t* and set the new value to the interpolated value of the old field at that point. Apply velocity BCs and free-slip.
3. **Pressure (first solve):** Compute divergence of **u**, set the pressure source as divergence/Δ*t* (zero in solid). Solve the Poisson equation −∇²*p* = source in the fluid (Gauss–Seidel, with liquid-only stencil). Apply pressure BCs. Then update **u** ← **u** − Δ*t* ∇*p* and apply free-slip.
4. **Diffusion:** Solve ∂*s*/∂*t* = *ν* ∇²*s* implicitly for each velocity component (Gauss–Seidel). Apply velocity BCs and free-slip.
5. **Pressure (second solve):** Same as step 3. Final **u** and *p* are used for the next step and for particle advection.

**Particles:** After the fluid step, each particle is advanced: sample **u** at the particle, compute Re and *C*<sub>*d*</sub>, then τ<sub>*p*</sub>; form *A* = −1/τ<sub>*p*</sub> and **B** from drag and gravity; integrate velocity analytically over Δ*t*; integrate position; if the new position is in the solid, reflect and damp velocity and re-advance position. Particles that end in solid or outside the domain are removed.

**Obstacles:** Drawn/erased as circles in physical space; the solid mask is set to 1 (solid) or 0 (fluid) in affected cells. After a draw/erase action, the mask is smoothed (3×3 box filter) and binarized at 0.5 so the interface is smooth and normals are stable.

The algorithms and equation ordering (vorticity → advect → pressure → diffuse → pressure, then particles) match a reference nodal implementation used for validation; the same governing equations and boundary treatments are used there.

---

### Project layout

```
index.html
css/
  index.css
js/
  index.js
  core/
    FluidApp.js
    math/
      MathUtils.js
    solvers/
      FluidSolver.js
      ParticleSolver.js
    obstacles/
      ObstacleManager.js
    graphics/
      ObstacleContour.js
      Renderer.js
  controls/
    Interaction.js
  utils/
    Utils.js
```

- **`index.html`** — Entry page: toolbar, BC panels, canvas, and script load order.
- **`css/index.css`** — Layout and styling for toolbar, BC panels, and canvas.
- **`js/index.js`** — Entry point: creates and starts the app.
- **`js/core/FluidApp.js`** — Configuration (grid, time step, viscosity, BC, particle params), state (fluid, particles), main loop (fixed-step integration, colormap limit every 10 steps), reset, BC UI, and toolbar wiring.
- **`js/core/math/MathUtils.js`** — Shared math: clamp, lerp, and small epsilon for safe division and comparisons.
- **`js/core/solvers/FluidSolver.js`** — Grid indexing, physical/grid coordinates, interpolation, derivatives (ddx, ddy, div, grad), BCs, vorticity, solid normals/tangent/projection, free-slip, advection, diffusion, pressure Poisson, and full fluid step.
- **`js/core/solvers/ParticleSolver.js`** — Drag coefficient *C*<sub>*d*</sub>(Re), reflection, linear ODE integration, and particle step (drag, gravity, collision, removal).
- **`js/core/obstacles/ObstacleManager.js`** — Draw/erase circles on the solid mask, smooth, binarize; solid normals are exposed via FluidSolver.
- **`js/core/graphics/ObstacleContour.js`** — Marching-squares contour of the obstacle mask for rendering.
- **`js/core/graphics/Renderer.js`** — Map physical domain to canvas, velocity magnitude to color (HSL), draw cells, arrows, and particles.
- **`js/controls/Interaction.js`** — Tool modes (obstacle, eraser, velocity, particles), pointer handlers, velocity brush, and particle injection.
- **`js/utils/Utils.js`** — General utilities (e.g. FPS meter).

---

## License

Copyright (c) 2026 Ilya Tsivilskiy

This repository is provided under an **all rights reserved** notice. No permission is granted to use, copy, modify, distribute, or exploit this software, including for commercial purposes, without prior written permission from the author. No patent license is granted.

See [LICENSE](./LICENSE) for the full terms.

