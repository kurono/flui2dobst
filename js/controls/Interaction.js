/**
 * User interaction: tool modes (obstacle, eraser, velocity, particles) and pointer event handlers.
 *
 * Obstacle/eraser: draw/erase circles on solid mask; on pointer up, smooth and binarize.
 * Velocity: add velocity in a brush (optionally along drag direction).
 * Particles: inject a cluster of particles at click/drag.
 *
 * @author Ilya Tsivilskiy
 * @see FluidSolver, ObstacleManager, ParticleSolver, Renderer
 */

/** @readonly */
const InteractionMode = Object.freeze({
  OBSTACLE: 'obstacle',
  ERASER: 'eraser',
  VELOCITY: 'velocity',
  PARTICLES: 'particles'
});

class Interaction {
  /**
   * @param {HTMLCanvasElement} canvas - Canvas to attach pointer events to
   * @param {Renderer} renderer - For canvasToPhys
   * @param {{ fluidState: Object, particles: Object }} state - Shared fluid and particle state
   * @param {{ grid: { res: number, L: number, h: number }, onObstaclesChanged?: function }} params - Grid params and optional callback when obstacles change (draw/erase stroke end)
   */
  constructor(canvas, renderer, state, params) {
    this.canvas = canvas;
    this.renderer = renderer;
    this.fluidState = state.fluidState;
    this.particles = state.particles;
    this.res = params.grid.res;
    this.L = params.grid.L;
    this.h = params.grid.h;
    this.onObstaclesChanged = params.onObstaclesChanged || (() => {});
    this.brushRadius = this.L * 0.05;
    this.particleInjectRadius = this.L * 0.08;
    this.particleInjectCount = 12;
    this.velocityStrength = 2;
    this.mode = InteractionMode.OBSTACLE;
    this.isDown = false;
    this.lastParticleTime = 0;
    this.#bindEvents();
  }

  /**
   * Map pointer event to physical (x,y) using canvas rect and renderer.canvasToPhys.
   * @param {PointerEvent} e - Pointer event
   * @returns {{ x: number, y: number }} Physical coordinates
   */
  #getPhys(e) {
    const rect = this.canvas.getBoundingClientRect();
    const scaleX = this.canvas.width / rect.width;
    const scaleY = this.canvas.height / rect.height;
    const cx = (e.clientX - rect.left) * scaleX;
    const cy = (e.clientY - rect.top) * scaleY;
    return this.renderer.canvasToPhys(cx, cy);
  }

  /**
   * Add a cluster of particles at (px, py) with random offset within particleInjectRadius.
   * @param {number} px - Center x
   * @param {number} py - Center y
   */
  #addParticleCluster(px, py) {
    for (let k = 0; k < this.particleInjectCount; k++) {
      const r = this.particleInjectRadius * Math.sqrt(Math.random());
      const theta = 2 * Math.PI * Math.random();
      const x = px + r * Math.cos(theta);
      const y = py + r * Math.sin(theta);
      ParticleSolver.addParticle(this.particles, x, y, 0, 0);
    }
  }

  /**
   * Add velocity (vx, vy) to fluid cells within brushRadius of (px, py), with falloff.
   * Skips solid cells.
   * @param {number} px - Brush center x
   * @param {number} py - Brush center y
   * @param {number} vx - Velocity x to add
   * @param {number} vy - Velocity y to add
   */
  #applyVelocityBrush(px, py, vx, vy) {
    const { fiy, fix } = FluidSolver.physToGrid(py, px, this.h, this.L);
    const I = FluidSolver.idx(this.res);
    const r = Math.ceil(this.brushRadius / this.h);
    const iy0 = Math.max(1, Math.floor(fiy - r));
    const iy1 = Math.min(this.res, Math.ceil(fiy + r));
    const ix0 = Math.max(1, Math.floor(fix - r));
    const ix1 = Math.min(this.res, Math.ceil(fix + r));
    for (let iy = iy0; iy <= iy1; iy++) {
      for (let ix = ix0; ix <= ix1; ix++) {
        if (this.fluidState.solid[I(iy, ix)] > 0.5) continue;
        const dx = (ix - 1) * this.h - this.L / 2 - px;
        const dy = (iy - 1) * this.h - this.L / 2 - py;
        if (Math.hypot(dx, dy) > this.brushRadius) continue;
        const falloff = 1 - (Math.hypot(dx, dy) / this.brushRadius) * 0.5;
        this.fluidState.Ux[I(iy, ix)] += vx * falloff * this.velocityStrength * 0.1;
        // Scale down so one stroke doesn't blow up the field
        this.fluidState.Uy[I(iy, ix)] += vy * falloff * this.velocityStrength * 0.1;
      }
    }
  }

  /** Pointer down: apply current tool at (px, py). */
  #onPointerDown = (e) => {
    e.preventDefault();
    this.isDown = true;
    const { x: px, y: py } = this.#getPhys(e);
    if (this.mode === InteractionMode.OBSTACLE) {
      ObstacleManager.drawCircle(this.fluidState.solid, px, py, this.brushRadius, this.res, this.L, this.h);
    } else if (this.mode === InteractionMode.ERASER) {
      ObstacleManager.eraseCircle(this.fluidState.solid, px, py, this.brushRadius, this.res, this.L, this.h);
    } else if (this.mode === InteractionMode.VELOCITY) {
      this.#applyVelocityBrush(px, py, 0, 1);
    } else if (this.mode === InteractionMode.PARTICLES) {
      this.#addParticleCluster(px, py);
    }
  };

  /** Pointer move while down: continue drawing/erasing/velocity/particles. */
  #onPointerMove = (e) => {
    if (!this.isDown) return;
    const { x: px, y: py } = this.#getPhys(e);
    if (this.mode === InteractionMode.OBSTACLE) {
      ObstacleManager.drawCircle(this.fluidState.solid, px, py, this.brushRadius, this.res, this.L, this.h);
    } else if (this.mode === InteractionMode.ERASER) {
      ObstacleManager.eraseCircle(this.fluidState.solid, px, py, this.brushRadius, this.res, this.L, this.h);
    } else if (this.mode === InteractionMode.VELOCITY) {
      const dx = e.movementX || 0;
      const dy = -(e.movementY || 0);
      const mag = Math.hypot(dx, dy) || 1;
      if (mag > 0) this.#applyVelocityBrush(px, py, (dx / mag) * 0.5, (dy / mag) * 0.5);
    } else if (this.mode === InteractionMode.PARTICLES && Date.now() - this.lastParticleTime > 50) {
      this.#addParticleCluster(px, py);
      this.lastParticleTime = Date.now();
    }
  };

  /** Pointer up: if obstacle/eraser, smooth and binarize solid mask and notify contour to rebuild. */
  #onPointerUp = () => {
    if (this.isDown && (this.mode === InteractionMode.OBSTACLE || this.mode === InteractionMode.ERASER)) {
      ObstacleManager.smooth(this.fluidState.solid, this.res, 1);
      ObstacleManager.binarize(this.fluidState.solid, this.res, 0.5);
      this.onObstaclesChanged();
    }
    this.isDown = false;
  };

  /** Attach pointerdown, pointermove, pointerup, pointerleave to canvas. */
  #bindEvents() {
    this.canvas.addEventListener('pointerdown', this.#onPointerDown, { passive: false });
    this.canvas.addEventListener('pointermove', this.#onPointerMove, { passive: false });
    this.canvas.addEventListener('pointerup', this.#onPointerUp);
    this.canvas.addEventListener('pointerleave', this.#onPointerUp);
  }

  /**
   * Set the active tool mode.
   * @param {string} mode - One of InteractionMode.OBSTACLE | ERASER | VELOCITY | PARTICLES
   */
  setMode(mode) {
    this.mode = mode;
  }

  /** @returns {string} Current tool mode */
  getMode() {
    return this.mode;
  }
}
