/**
 * Lagrangian particle solver: drag (Stokes + Morsi Cd), gravity, solid reflection.
 *
 * Particle equation:
 *   dUp/dt = (U - Up)/tp + g*(rhop - rho)/rhop + a
 *   tp = (24/18)*rhop*dp²/(mu*Cd*Re) — relaxation time
 *   Re = |Ur|*dp/nu, |Ur| = |U - Up|
 *
 * Uses FluidSolver.sample() for U and solid/normal at particle position.
 *
 * @author Ilya Tsivilskiy
 */
class ParticleSolver {
  /**
   * Drag coefficient Cd(Re) from Morsi & Alexander (1972).
   * @param {number} Re - Particle Reynolds number
   * @returns {number} Cd
   */
  static #dragCoeff(Re) {
    if (Re < MathUtils.EPS) Re = MathUtils.EPS;
    const Re2 = Re * Re;
    let cd = 24 / Re;
    const inRange = (a, lo, hi) => MathUtils.inRange(a, lo, hi);
    if (inRange(Re, 0.1, 1)) cd = 22.73 / Re + 0.0903 / Re2 + 3.69;
    if (inRange(Re, 1, 10)) cd = 29.1667 / Re - 3.8889 / Re2 + 1.222;
    if (inRange(Re, 10, 100)) cd = 46.5 / Re - 116.67 / Re2 + 0.6167;
    if (inRange(Re, 100, 1000)) cd = 98.33 / Re - 2778 / Re2 + 0.3644;
    if (inRange(Re, 1000, 5000)) cd = 148.62 / Re - 4.75e4 / Re2 + 0.357;
    if (inRange(Re, 5000, 10000)) cd = -490.546 / Re + 57.87e4 / Re2 + 0.46;
    if (inRange(Re, 10000, 50000)) cd = -1662.5 / Re + 5.4167e6 / Re2 + 0.5191;
    return cd;
  }

  /**
   * Reflect vector (dx,dy) about unit-normal (nx,ny): v - 2(v·n)n.
   * @param {number} dx - Vector x
   * @param {number} dy - Vector y
   * @param {number} nx - Normal x (will be normalized)
   * @param {number} ny - Normal y
   * @returns {{ rx: number, ry: number }} Reflected vector
   */
  static #reflect(dx, dy, nx, ny) {
    const mag = Math.hypot(nx, ny) + MathUtils.EPS;
    nx /= mag;
    ny /= mag;
    const dn = dx * nx + dy * ny;
    return { rx: dx - 2 * dn * nx, ry: dy - 2 * dn * ny };
  }

  /**
   * Solve du/dt = A*u + B with IC u(0)=uprev. Analytical: u = (B/A + uprev)*exp(A*dt) - B/A.
   * @param {number} uprev - Initial value
   * @param {number} A - Linear coefficient
   * @param {number} B - Forcing
   * @param {number} dt - Time step
   * @returns {number} u(dt)
   */
  static #solveTransport(uprev, A, B, dt) {
    return (B / A + uprev) * Math.exp(A * dt) - B / A;
  }

  /**
   * Create empty particle state: arrays for x, y, vx, vy.
   * @returns {{ x: number[], y: number[], vx: number[], vy: number[] }}
   */
  static createParticles() {
    return { x: [], y: [], vx: [], vy: [] };
  }

  /**
   * Maximum velocity magnitude over all particles.
   * @param {{ vx: number[], vy: number[] }} particles - Particle state
   * @returns {number} Max |v|
   */
  static getMaxParticleVelocity(particles) {
    let maxMag = 0;
    for (let i = 0; i < particles.vx.length; i++) {
      const m = Math.hypot(particles.vx[i], particles.vy[i]);
      if (m > maxMag) maxMag = m;
    }
    return maxMag;
  }

  /**
   * Add one particle at (x,y) with optional initial velocity.
   * @param {{ x: number[], y: number[], vx: number[], vy: number[] }} particles - Particle state (mutated)
   * @param {number} x - Physical x
   * @param {number} y - Physical y
   * @param {number} [vx=0] - Initial vx
   * @param {number} [vy=0] - Initial vy
   */
  static addParticle(particles, x, y, vx = 0, vy = 0) {
    particles.x.push(x);
    particles.y.push(y);
    particles.vx.push(vx);
    particles.vy.push(vy);
  }

  /**
   * Advance particles one step: sample U at position, compute drag + gravity, integrate velocity
   * (implicit transport), integrate position, resolve solid collision (reflect + damp), then
   * remove particles that are in solid or out of domain.
   * @param {{ x: number[], y: number[], vx: number[], vy: number[] }} particles - Particle state (mutated)
   * @param {{ Ux: Float32Array, Uy: Float32Array, solid: Float32Array }} fluidState - Fluid state
   * @param {Object} params - particleParams (grid, particles { nu, rho, mu, dp, rhop, gx, gy, pVelDamp }, solidNormal { snx, sny })
   * @param {number} dt - Time step
   */
  static stepParticles(particles, fluidState, params, dt) {
    const { res, L, h } = params.grid;
    const { nu, rho, mu, dp, rhop, gx, gy, pVelDamp } = params.particles;
    const { Ux, Uy, solid } = fluidState;
    const { snx, sny } = params.solidNormal;
    const rhor = rho / rhop;
    for (let i = 0; i < particles.x.length; i++) {
      let x = particles.x[i], y = particles.y[i], vx = particles.vx[i], vy = particles.vy[i];
      const ux = FluidSolver.sample(Ux, y, x, h, L, res);
      const uy = FluidSolver.sample(Uy, y, x, h, L, res);
      const Ur = Math.hypot(ux - vx, uy - vy) + MathUtils.EPS;
      const Re = (Ur * dp) / nu;
      const Cd = this.#dragCoeff(Re);
      const tp = (24 / 18) * rhop * (dp * dp) / (mu * Cd * Re + MathUtils.EPS);
      const A = -1 / tp;
      const Bx = -A * ux + gx * (1 - rhor);
      const By = -A * uy + gy * (1 - rhor);
      vx = this.#solveTransport(vx, A, Bx, dt);
      vy = this.#solveTransport(vy, A, By, dt);
      const xNew = x + dt * vx;
      const yNew = y + dt * vy;
      const s = FluidSolver.sample(solid, yNew, xNew, h, L, res);
      const nx = FluidSolver.sample(snx, yNew, xNew, h, L, res);
      const ny = FluidSolver.sample(sny, yNew, xNew, h, L, res);
      // If new position is inside solid: reflect velocity, apply damping, re-advance from old position
      if (s > 0.5) {
        const { rx, ry } = this.#reflect(vx, vy, nx, ny);
        particles.vx[i] = rx * (1 - pVelDamp);
        particles.vy[i] = ry * (1 - pVelDamp);
        particles.x[i] = x + dt * particles.vx[i];
        particles.y[i] = y + dt * particles.vy[i];
      } else {
        particles.x[i] = xNew;
        particles.y[i] = yNew;
        particles.vx[i] = vx;
        particles.vy[i] = vy;
      }
    }
    // Remove particles that ended up in solid or outside domain
    const half = L / 2;
    for (let i = particles.x.length - 1; i >= 0; i--) {
      const px = particles.x[i];
      const py = particles.y[i];
      const inSolid = FluidSolver.sample(solid, py, px, h, L, res) > 0.5;
      const outOfDomain = px < -half || px > half || py < -half || py > half;
      if (inSolid || outOfDomain) {
        particles.x.splice(i, 1);
        particles.y.splice(i, 1);
        particles.vx.splice(i, 1);
        particles.vy.splice(i, 1);
      }
    }
  }
}
