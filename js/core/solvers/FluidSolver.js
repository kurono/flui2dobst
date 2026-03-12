/**
 * 2D incompressible Navier–Stokes fluid solver on a Cartesian grid with obstacles.
 *
 * Equations:
 *   Momentum: dU/dt + (U·∇)U = -∇P + ν∇²U
 *   Continuity: ∇·U = 0
 * Free-slip on solid: velocity projected onto tangent; pressure uses liquid-only Laplacian.
 *
 * Provides: grid indexing, phys/grid coords, sample/lerp, derivatives (ddx, ddy, div, grad),
 * BCs (zeroGrad, velocity, pressure), vorticity confinement, solid normals/tangent/project,
 * applySolidU, semi-Lagrangian advection, implicit diffusion, pressure Poisson, stepFluid.
 *
 * @author Ilya Tsivilskiy
 */
class FluidSolver {
  /**
   * Linear index from row/col. Grid is row-major: (iy, ix) -> (iy-1)*res + (ix-1).
   * @param {number} res - Grid resolution
   * @returns {(iy: number, ix: number) => number} Index function
   */
  static idx(res) {
    return (iy, ix) => (iy - 1) * res + (ix - 1);
  }

  /**
   * Convert grid node indices to physical coordinates (center of domain at 0).
   * y = (iy-1)*h - L/2, x = (ix-1)*h - L/2.
   * @param {number} iy - Row index (1-based)
   * @param {number} ix - Column index (1-based)
   * @param {number} h - Cell size
   * @param {number} L - Domain size
   * @returns {{ y: number, x: number }}
   */
  static gridToPhys(iy, ix, h, L) {
    return { y: (iy - 1) * h - L / 2, x: (ix - 1) * h - L / 2 };
  }

  /**
   * Convert physical (x,y) to float grid indices (for interpolation).
   * fiy = (L/2+y)/h+1, fix = (L/2+x)/h+1.
   * @param {number} y - Physical y
   * @param {number} x - Physical x
   * @param {number} h - Cell size
   * @param {number} L - Domain size
   * @returns {{ fiy: number, fix: number }}
   */
  static physToGrid(y, x, h, L) {
    return { fiy: (L / 2 + y) / h + 1, fix: (L / 2 + x) / h + 1 };
  }

  /**
   * Sample scalar field F at physical (x,y) using bilinear interpolation from 4 nodes.
   * Uses physToGrid then bilinear lerp from 4 neighbors (y then x).
   * @param {Float32Array} F - Scalar field (length res*res)
   * @param {number} y - Physical y
   * @param {number} x - Physical x
   * @param {number} h - Cell size
   * @param {number} L - Domain size
   * @param {number} res - Resolution
   * @returns {number} Sampled value
   */
  static sample(F, y, x, h, L, res) {
    const { fiy, fix } = this.physToGrid(y, x, h, L);
    const ixr = MathUtils.clamp(Math.ceil(fix), 1, res);
    const ixl = MathUtils.clamp(Math.floor(fix), 1, res);
    const iyt = MathUtils.clamp(Math.ceil(fiy), 1, res);
    const iyb = MathUtils.clamp(Math.floor(fiy), 1, res);
    const I = this.idx(res);
    const ftl = F[I(iyt, ixl)];
    const ftr = F[I(iyt, ixr)];
    const fbl = F[I(iyb, ixl)];
    const fbr = F[I(iyb, ixr)];
    const fl = MathUtils.lerp(fbl, ftl, iyb, fiy, iyt);
    const fr = MathUtils.lerp(fbr, ftr, iyb, fiy, iyt);
    return MathUtils.lerp(fl, fr, ixl, fix, ixr);
  }

  /**
   * Zero-gradient (Neumann) BC: copy from interior to boundary.
   * @param {Float32Array} F - Field (modified in place)
   * @param {number} res - Resolution
   */
  static zeroGradBC(F, res) {
    const I = this.idx(res);
    for (let ix = 1; ix <= res; ix++) {
      F[I(1, ix)] = F[I(2, ix)];
      F[I(res, ix)] = F[I(res - 1, ix)];
    }
    for (let iy = 1; iy <= res; iy++) {
      F[I(iy, 1)] = F[I(iy, 2)];
      F[I(iy, res)] = F[I(iy, res - 1)];
    }
  }

  /**
   * Partial derivative dF/dx via central differences.
   * @param {Float32Array} F - Scalar field
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {Float32Array} dF/dx
   */
  static ddx(F, h, res) {
    const out = new Float32Array(res * res);
    const I = this.idx(res);
    for (let iy = 2; iy < res; iy++) {
      for (let ix = 2; ix < res; ix++) {
        out[I(iy, ix)] = (F[I(iy, ix + 1)] - F[I(iy, ix - 1)]) / (2 * h);
      }
    }
    this.zeroGradBC(out, res);
    return out;
  }

  /**
   * Partial derivative dF/dy via central differences.
   * @param {Float32Array} F - Scalar field
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {Float32Array} dF/dy
   */
  static ddy(F, h, res) {
    const out = new Float32Array(res * res);
    const I = this.idx(res);
    for (let iy = 2; iy < res; iy++) {
      for (let ix = 2; ix < res; ix++) {
        out[I(iy, ix)] = (F[I(iy + 1, ix)] - F[I(iy - 1, ix)]) / (2 * h);
      }
    }
    this.zeroGradBC(out, res);
    return out;
  }

  /**
   * Divergence of (Fx, Fy): div = dFx/dx + dFy/dy.
   * @param {Float32Array} Fx - X component
   * @param {Float32Array} Fy - Y component
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {Float32Array} Divergence field
   */
  static div(Fx, Fy, h, res) {
    const dFxDx = this.ddx(Fx, h, res);
    const dFyDy = this.ddy(Fy, h, res);
    const out = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) out[i] = dFxDx[i] + dFyDy[i];
    this.zeroGradBC(out, res);
    return out;
  }

  /**
   * Gradient of scalar F: (dF/dx, dF/dy).
   * @param {Float32Array} F - Scalar field
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {{ Gx: Float32Array, Gy: Float32Array }}
   */
  static grad(F, h, res) {
    return { Gx: this.ddx(F, h, res), Gy: this.ddy(F, h, res) };
  }

  /**
   * Apply velocity BC per side: either prescribed ux/uy or zero-gradient.
   * Prescribed or zero-gradient per edge.
   * @param {Float32Array} F - Velocity component (modified in place)
   * @param {number} res - Resolution
   * @param {Object} bc - BC config { top, bottom, left, right } with ux, uy, zeroGrad
   * @param {'x'|'y'} comp - Which component ('x' => ux, 'y' => uy)
   */
  static velocityBC(F, res, bc, comp) {
    const I = this.idx(res);
    const val = (side) => comp === 'x' ? side.ux : side.uy;
    for (let ix = 1; ix <= res; ix++) {
      F[I(1, ix)] = bc.bottom.zeroGrad ? F[I(2, ix)] : val(bc.bottom);
      F[I(res, ix)] = bc.top.zeroGrad ? F[I(res - 1, ix)] : val(bc.top);
    }
    for (let iy = 1; iy <= res; iy++) {
      F[I(iy, 1)] = bc.left.zeroGrad ? F[I(iy, 2)] : val(bc.left);
      F[I(iy, res)] = bc.right.zeroGrad ? F[I(iy, res - 1)] : val(bc.right);
    }
  }

  /**
   * Pressure BC: zero gradient on top/bottom, zero on left/right (or copy).
   * Top/bottom: copy from interior; left/right: zero or copy.
   * @param {Float32Array} F - Pressure field (modified in place)
   * @param {number} res - Resolution
   */
  static pressureBC(F, res) {
    const I = this.idx(res);
    for (let ix = 1; ix <= res; ix++) {
      F[I(1, ix)] = F[I(2, ix)];
      F[I(res, ix)] = F[I(res - 1, ix)];
    }
    for (let iy = 1; iy <= res; iy++) {
      F[I(iy, 1)] = F[I(iy, 2)];
      F[I(iy, res)] = F[I(iy, res - 1)];
    }
  }

  /**
   * Z-component of curl of 2D velocity (U,V): dV/dx - dU/dy.
   * @param {Float32Array} U - X velocity
   * @param {Float32Array} V - Y velocity
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {Float32Array} Vorticity (scalar)
   */
  static curlZ(U, V, h, res) {
    const R = new Float32Array(res * res);
    const I = this.idx(res);
    for (let iy = 2; iy < res; iy++) {
      for (let ix = 2; ix < res; ix++) {
        R[I(iy, ix)] = (V[I(iy, ix + 1)] - V[I(iy, ix - 1)]) / (2 * h) -
          (U[I(iy + 1, ix)] - U[I(iy - 1, ix)]) / (2 * h);
      }
    }
    this.zeroGradBC(R, res);
    return R;
  }

  /**
   * Vorticity confinement force: amplifies vorticity to preserve small-scale structure.
   * @param {Float32Array} Ux - X velocity
   * @param {Float32Array} Uy - Y velocity
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @param {number} [scale=15] - Strength
   * @returns {{ Wx: Float32Array, Wy: Float32Array }} Acceleration field
   */
  static vorticityAcceleration(Ux, Uy, h, res, scale = 5) {
    const vort = this.curlZ(Ux, Uy, h, res);
    const { Gx: vdx, Gy: vdy } = this.grad(vort, h, res);
    const mag = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) mag[i] = Math.hypot(vdx[i], vdy[i]) + MathUtils.EPS;
    const Wx = new Float32Array(res * res);
    const Wy = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) {
      const vx = vdx[i] / mag[i];
      const vy = vdy[i] / mag[i];
      Wx[i] = -vy * vort[i] * h * scale;
      Wy[i] = vx * vort[i] * h * scale;
    }
    return { Wx, Wy };
  }

  /**
   * Outer unit normal to solid (from gradient of solid mask, negated and normalized).
   * @param {Float32Array} solid - Solid mask (1 = solid, 0 = fluid)
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {{ snx: Float32Array, sny: Float32Array }}
   */
  static solidNormal(solid, h, res) {
    const { Gx: nx, Gy: ny } = this.grad(solid, h, res);
    const snx = new Float32Array(res * res);
    const sny = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) {
      const m = Math.hypot(nx[i], ny[i]) + MathUtils.EPS;
      snx[i] = -nx[i] / m;
      sny[i] = -ny[i] / m;
    }
    return { snx, sny };
  }

  /**
   * Tangent vector to solid surface: 90° rotation of normal (stx = -sny, sty = snx).
   * @param {Float32Array} snx - Normal x
   * @param {Float32Array} sny - Normal y
   * @param {number} res - Resolution
   * @returns {{ stx: Float32Array, sty: Float32Array }}
   */
  static tangent(snx, sny, res) {
    const stx = new Float32Array(res * res);
    const sty = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) {
      stx[i] = -sny[i];
      sty[i] = snx[i];
    }
    return { stx, sty };
  }

  /**
   * Project vector (vx,vy) onto direction (dx,dy).
   * @param {number} vx - Vector x
   * @param {number} vy - Vector y
   * @param {number} dx - Direction x
   * @param {number} dy - Direction y
   * @returns {{ px: number, py: number }} Projected vector
   */
  static project(vx, vy, dx, dy) {
    const d = vx * dx + vy * dy;
    const dd = dx * dx + dy * dy + MathUtils.EPS;
    return { px: (d / dd) * dx, py: (d / dd) * dy };
  }

  /**
   * Free-slip BC on solid: zero velocity inside solid; at liquid nodes adjacent to solid,
   * project velocity onto local tangent.
   * @param {Float32Array} Ux - X velocity (read)
   * @param {Float32Array} Uy - Y velocity (read)
   * @param {Float32Array} solid - Solid mask
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {{ Ux: Float32Array, Uy: Float32Array }} New velocity arrays
   */
  static applySolidU(Ux, Uy, solid, h, res) {
    const I = this.idx(res);
    const Uxs = Ux.slice();
    const Uys = Uy.slice();
    for (let i = 0; i < res * res; i++) {
      if (solid[i] > 0.5) { Uxs[i] = 0; Uys[i] = 0; }
    }
    const { snx, sny } = this.solidNormal(solid, h, res);
    const { stx, sty } = this.tangent(snx, sny, res);
    // For each solid node, project neighbor liquid velocities onto local tangent
    for (let iy = 2; iy < res; iy++) {
      for (let ix = 2; ix < res; ix++) {
        if (solid[I(iy, ix)] <= 0.5) continue;
        for (let oy = -1; oy <= 1; oy++) {
          for (let ox = -1; ox <= 1; ox++) {
            const iyn = iy + oy, ixn = ix + ox;
            if (solid[I(iyn, ixn)] > 0.5) continue;
            const { px, py } = this.project(Ux[I(iyn, ixn)], Uy[I(iyn, ixn)], stx[I(iy, ix)], sty[I(iy, ix)]);
            Uxs[I(iyn, ixn)] = px;
            Uys[I(iyn, ixn)] = py;
          }
        }
      }
    }
    return { Ux: Uxs, Uy: Uys };
  }

  /**
   * Semi-Lagrangian advection: S(r) = Sprev(r - U*dt).
   * @param {Float32Array} Sprev - Scalar field at previous time
   * @param {Float32Array} Ux - X velocity
   * @param {Float32Array} Uy - Y velocity
   * @param {number} dt - Time step
   * @param {number} h - Cell size
   * @param {number} L - Domain size
   * @param {number} res - Resolution
   * @param {function(Float32Array, number): void} bcFunc - BC to apply after advection
   * @returns {Float32Array} Advected field
   */
  static advect(Sprev, Ux, Uy, dt, h, L, res, bcFunc) {
    const out = new Float32Array(res * res);
    const I = this.idx(res);
    for (let iy = 2; iy < res; iy++) {
      for (let ix = 2; ix < res; ix++) {
        const { y, x } = this.gridToPhys(iy, ix, h, L);
        // Backtrace: sample Sprev at (x - U*dt, y - V*dt)
        out[I(iy, ix)] = this.sample(Sprev, y - Uy[I(iy, ix)] * dt, x - Ux[I(iy, ix)] * dt, h, L, res);
      }
    }
    bcFunc(out, res);
    return out;
  }

  /**
   * Implicit diffusion: dS/dt = nu*Laplacian(S). Gauss–Seidel iteration.
   * @param {Float32Array} Sprev - Field at previous time
   * @param {number} nu - Kinematic viscosity
   * @param {number} dt - Time step
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @param {function(Float32Array, number): void} bcFunc - BC each iteration
   * @param {number} maxit - Max iterations
   * @returns {Float32Array} Diffused field
   */
  static diffuse(Sprev, nu, dt, h, res, bcFunc, maxit) {
    const out = Sprev.slice();
    const I = this.idx(res);
    const fac = h * h + 4 * nu * dt;
    for (let it = 0; it < maxit; it++) {
      for (let iy = 2; iy < res; iy++) {
        for (let ix = 2; ix < res; ix++) {
          const sum4 = out[I(iy + 1, ix)] + out[I(iy - 1, ix)] + out[I(iy, ix + 1)] + out[I(iy, ix - 1)];
          out[I(iy, ix)] = (nu * dt * sum4 + h * h * Sprev[I(iy, ix)]) / fac;
        }
      }
      bcFunc(out, res);
    }
    return out;
  }

  /**
   * Solve continuity: compute pressure from div(U)/dt, solve Poisson (liquid-only Laplacian),
   * then apply pressure gradient to velocity.
   * @param {Float32Array} Ux - X velocity (before pressure)
   * @param {Float32Array} Uy - Y velocity (before pressure)
   * @param {Float32Array} Pprev - Pressure from previous step
   * @param {Float32Array} solid - Solid mask (Psrc zero in solid)
   * @param {number} dt - Time step
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @param {number} maxit - Gauss–Seidel iterations
   * @returns {{ Ux: Float32Array, Uy: Float32Array, P: Float32Array }}
   */
  static solvePressure(Ux, Uy, Pprev, solid, dt, h, res, maxit) {
    const divU = this.div(Ux, Uy, h, res);
    const Psrc = new Float32Array(res * res);
    const I = this.idx(res);
    for (let i = 0; i < res * res; i++) Psrc[i] = (solid[i] > 0.5) ? 0 : divU[i] / dt;
    const P = Pprev.slice();
    const liquid = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) liquid[i] = solid[i] > 0.5 ? 0 : 1;
    // Gauss–Seidel for Poisson: only liquid neighbors contribute to Laplacian
    for (let it = 0; it < maxit; it++) {
      for (let iy = 2; iy < res; iy++) {
        for (let ix = 2; ix < res; ix++) {
          const sum4 = liquid[I(iy + 1, ix)] * P[I(iy + 1, ix)] + liquid[I(iy - 1, ix)] * P[I(iy - 1, ix)] +
            liquid[I(iy, ix + 1)] * P[I(iy, ix + 1)] + liquid[I(iy, ix - 1)] * P[I(iy, ix - 1)];
          const suml = liquid[I(iy + 1, ix)] + liquid[I(iy - 1, ix)] + liquid[I(iy, ix + 1)] + liquid[I(iy, ix - 1)];
          if (suml > 0) P[I(iy, ix)] = (sum4 - Psrc[I(iy, ix)] * h * h) / suml;
        }
      }
      this.pressureBC(P, res);
    }
    const { Gx: dPdx, Gy: dPdy } = this.grad(P, h, res);
    const UxNew = new Float32Array(res * res);
    const UyNew = new Float32Array(res * res);
    for (let i = 0; i < res * res; i++) {
      UxNew[i] = Ux[i] - dPdx[i] * dt;
      UyNew[i] = Uy[i] - dPdy[i] * dt;
    }
    return { Ux: UxNew, Uy: UyNew, P };
  }

  /**
   * One time step: vorticity → applySolidU → advect → applySolidU → pressure → applySolidU
   * → diffuse → applySolidU → pressure.
   * @param {{ Ux: Float32Array, Uy: Float32Array, P: Float32Array, solid: Float32Array }} state - Current state
   * @param {number} dt - Time step
   * @param {Object} params - fluidParams (res, L, h, nu, maxit, vorticityScale, bc)
   * @returns {{ Ux: Float32Array, Uy: Float32Array, P: Float32Array, solid: Float32Array }}
   */
  static stepFluid(state, dt, params) {
    const { res, L, h, nu, maxit, vorticityScale, bc } = params;
    const bcX = (F, r) => this.velocityBC(F, r, bc, 'x');
    const bcY = (F, r) => this.velocityBC(F, r, bc, 'y');
    const { Ux, Uy, P, solid } = state;
    let UxCur = Ux.slice();
    let UyCur = Uy.slice();
    let Pcur = P.slice();
    // 1. Add vorticity confinement
    const { Wx, Wy } = this.vorticityAcceleration(Ux, Uy, h, res, vorticityScale);
    for (let i = 0; i < res * res; i++) {
      UxCur[i] = Ux[i] + dt * Wx[i];
      UyCur[i] = Uy[i] + dt * Wy[i];
    }
    let applied = this.applySolidU(UxCur, UyCur, solid, h, res);
    UxCur = applied.Ux;
    UyCur = applied.Uy;
    // 2. Advect
    UxCur = this.advect(UxCur, Ux, Uy, dt, h, L, res, bcX);
    UyCur = this.advect(UyCur, Ux, Uy, dt, h, L, res, bcY);
    applied = this.applySolidU(UxCur, UyCur, solid, h, res);
    UxCur = applied.Ux;
    UyCur = applied.Uy;
    // 3. First pressure solve
    const press1 = this.solvePressure(UxCur, UyCur, Pcur, solid, dt, h, res, maxit);
    UxCur = press1.Ux;
    UyCur = press1.Uy;
    Pcur = press1.P;
    applied = this.applySolidU(UxCur, UyCur, solid, h, res);
    UxCur = applied.Ux;
    UyCur = applied.Uy;
    // 4. Diffuse
    UxCur = this.diffuse(UxCur, UxCur, nu, dt, h, res, bcX, maxit);
    UyCur = this.diffuse(UyCur, UyCur, nu, dt, h, res, bcY, maxit);
    applied = this.applySolidU(UxCur, UyCur, solid, h, res);
    UxCur = applied.Ux;
    UyCur = applied.Uy;
    // 5. Second pressure solve
    const press2 = this.solvePressure(UxCur, UyCur, Pcur, solid, dt, h, res, maxit);
    return { Ux: press2.Ux, Uy: press2.Uy, P: press2.P, solid };
  }

  /**
   * Allocate fluid state arrays: Ux, Uy, P, solid.
   * @param {number} res - Resolution
   * @param {number} L - Domain size
   * @returns {{ Ux: Float32Array, Uy: Float32Array, P: Float32Array, solid: Float32Array }}
   */
  static createFluidState(res, L) {
    const n = res * res;
    return {
      Ux: new Float32Array(n),
      Uy: new Float32Array(n),
      P: new Float32Array(n),
      solid: new Float32Array(n)
    };
  }

  /**
   * Maximum velocity magnitude over liquid cells (solid excluded).
   * @param {Float32Array} Ux - X velocity
   * @param {Float32Array} Uy - Y velocity
   * @param {Float32Array} solid - Solid mask
   * @param {number} res - Resolution
   * @returns {number} Max |U|
   */
  static getMaxFluidVelocity(Ux, Uy, solid, res) {
    let maxMag = 0;
    for (let i = 0; i < res * res; i++) {
      if (solid[i] > 0.5) continue;
      const m = Math.hypot(Ux[i], Uy[i]);
      if (m > maxMag) maxMag = m;
    }
    return maxMag;
  }
}
