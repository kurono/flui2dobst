/**
 * Canvas 2D rendering: velocity magnitude as colored cells, velocity arrows, solid cells, particles.
 *
 * Domain is centered and scaled to fit the canvas; y is flipped so that +y is up (physical convention).
 * Draws scalar field (e.g. |U|), velocity arrows, and particles.
 *
 * @author Ilya Tsivilskiy
 */
class Renderer {
  /**
   * @param {HTMLCanvasElement} canvas - Target canvas
   * @param {number} L - Domain size [m] (used for scale and extent)
   */
  /** Upsample factor for velocity magnitude drawing (bilinear). */
  static UPSAMPLE = 2;

  constructor(canvas, L) {
    this.canvas = canvas;
    this.L = L;
    this.ctx = canvas.getContext('2d');
    /** Cached contour segments; rebuilt only when obstacles change (user draw/erase/reset). */
    this._cachedContourSegments = [];
    this._cachedContourRes = 0;
  }

  /**
   * Bilinearly sample velocity magnitude at physical (x, y) from grid. Does not modify solver state.
   * @param {number} x - Physical x
   * @param {number} y - Physical y
   * @param {Float32Array} mag - Magnitude grid (length res*res), mag[i] at node (iy, ix), i = (iy-1)*res + (ix-1)
   * @param {number} res - Grid resolution
   * @param {number} h - Cell size
   * @returns {number}
   */
  #sampleMagnitudeBilinear(x, y, mag, res, h) {
    const fix = (x + this.L / 2) / h + 1;
    const fiy = (y + this.L / 2) / h + 1;
    const ixl = Math.max(1, Math.min(res - 1, Math.floor(fix)));
    const ixr = Math.max(1, Math.min(res, ixl + 1));
    const iyb = Math.max(1, Math.min(res - 1, Math.floor(fiy)));
    const iyt = Math.max(1, Math.min(res, iyb + 1));
    const idx = (iy, ix) => (iy - 1) * res + (ix - 1);
    const tx = fix - ixl;
    const ty = fiy - iyb;
    const mbl = mag[idx(iyb, ixl)];
    const mbr = mag[idx(iyb, ixr)];
    const mtl = mag[idx(iyt, ixl)];
    const mtr = mag[idx(iyt, ixr)];
    const mb = mbl + (mbr - mbl) * tx;
    const mt = mtl + (mtr - mtl) * tx;
    return mb + (mt - mb) * ty;
  }

  /**
   * Convert physical (x,y) to canvas pixel coordinates. Origin at center; y flipped.
   * @param {number} x - Physical x
   * @param {number} y - Physical y
   * @returns {{ cx: number, cy: number }} Canvas coordinates
   */
  #physToCanvas(x, y) {
    const scale = Math.min(this.canvas.width, this.canvas.height) / (this.L * 1.1);
    return {
      cx: this.canvas.width / 2 + x * scale,
      cy: this.canvas.height / 2 - y * scale
    };
  }

  /**
   * Convert canvas pixel coordinates to physical (x,y). Inverse of #physToCanvas.
   * @param {number} cx - Canvas x
   * @param {number} cy - Canvas y
   * @returns {{ x: number, y: number }} Physical coordinates
   */
  canvasToPhys(cx, cy) {
    const scale = Math.min(this.canvas.width, this.canvas.height) / (this.L * 1.1);
    return {
      x: (cx - this.canvas.width / 2) / scale,
      y: (this.canvas.height / 2 - cy) / scale
    };
  }

  /**
   * Draw obstacle outline using marching-squares contour (iso = 0.5).
   * Uses cached segments unless obstaclesDirty; converts grid coords to canvas.
   */
  #drawObstacleContour(solid, res, h, scale, ox, oy, segments) {
    if (segments.length === 0) return;

    const cellW = (this.L * scale) / res;
    const cellH = (this.L * scale) / res;
    // Align with cell drawing: row 0 = iy=1 = bottom of first row; cell top at oy+(L-iy*h)*scale, so node at row is at top+cellH
    const gridToCanvas = (col, row) => ({
      cx: ox + col * cellW,
      cy: oy + (this.L - (row + 1) * h) * scale + cellH
    });

    this.ctx.strokeStyle = '#000';
    this.ctx.lineWidth = Math.max(3, (this.L * scale) / res * 1.2);
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    this.ctx.beginPath();

    for (const { p0, p1 } of segments) {
      const a = gridToCanvas(p0.x, p0.y);
      const b = gridToCanvas(p1.x, p1.y);
      if (![a.cx, a.cy, b.cx, b.cy].every(Number.isFinite)) continue;
      this.ctx.moveTo(a.cx, a.cy);
      this.ctx.lineTo(b.cx, b.cy);
    }
    this.ctx.stroke();
  }

  /**
   * Map velocity magnitude to HSL color (blue = low, red = high).
   * @param {number} mag - |U|
   * @param {number} maxMag - Scale reference (avoids division by zero)
   * @returns {string} CSS color
   */
  #magnitudeToColor(mag, maxMag) {
    const t = Math.min(1, mag / (maxMag + 1e-6));
    return `hsl(${(1 - t) * 240}, 80%, 50%)`;
  }

  /**
   * Draw one frame: background, cell-colored velocity magnitude, solid cells, domain border,
   * velocity arrows (subsampled), and particles as black dots.
   * @param {{ Ux: Float32Array, Uy: Float32Array, solid: Float32Array }} fluidState - Current fluid state
   * @param {{ x: number[], y: number[] }} particles - Particle positions
   * @param {number} res - Grid resolution
   * @param {number} h - Cell size
   * @param {number} [maxMag] - Velocity scale for colormap (if omitted, computed from current field)
   * @param {{ obstaclesDirty?: boolean }} [options] - If obstaclesDirty, contour is recomputed and cached
   */
  render(fluidState, particles, res, h, maxMag, options = {}) {
    const { Ux, Uy, solid } = fluidState;
    const obstaclesDirty = options.obstaclesDirty === true;
    if (obstaclesDirty || this._cachedContourRes !== res) {
      this._cachedContourSegments = ObstacleContour.build(solid, res, 0.5);
      this._cachedContourRes = res;
    }
    const contourSegments = this._cachedContourSegments;
    const width = this.canvas.width;
    const height = this.canvas.height;
    const scale = Math.min(width, height) / (this.L * 1.1);
    const ox = width / 2 - (this.L / 2) * scale;
    const oy = height / 2 - (this.L / 2) * scale;
    const cellW = (this.L * scale) / res;
    const cellH = (this.L * scale) / res;
    if (maxMag == null) {
      maxMag = 0.01;
      for (let i = 0; i < Ux.length; i++) {
        const m = Math.hypot(Ux[i], Uy[i]);
        if (m > maxMag) maxMag = m;
      }
    }
    const n = res * res;
    const mag = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      mag[i] = solid[i] > 0.5 ? 0 : Math.hypot(Ux[i], Uy[i]);
    }
    this.ctx.fillStyle = '#1a1a2e';
    this.ctx.fillRect(0, 0, width, height);
    const k = Renderer.UPSAMPLE;
    const subW = cellW / k;
    const subH = cellH / k;
    for (let sy = 0; sy < res * k; sy++) {
      for (let sx = 0; sx < res * k; sx++) {
        const x = -this.L / 2 + (sx + 0.5) * (this.L / (res * k));
        const y = -this.L / 2 + (sy + 0.5) * (this.L / (res * k));
        const sample = this.#sampleMagnitudeBilinear(x, y, mag, res, h);
        this.ctx.fillStyle = this.#magnitudeToColor(sample, maxMag);
        const cx = ox + sx * subW;
        const cy = oy + (this.L * scale) - (sy + 1) * subH;
        this.ctx.fillRect(cx, cy, subW + 1, subH + 1);
      }
    }
    // Obstacle bounds: marching-squares contour at solid/fluid interface
    this.#drawObstacleContour(solid, res, h, scale, ox, oy, contourSegments);
    // Domain border
    this.ctx.strokeStyle = '#444';
    this.ctx.lineWidth = 1;
    this.ctx.strokeRect(ox, oy, this.L * scale, this.L * scale);
    // Velocity arrows (subsampled for readability)
    const stride = Math.max(1, Math.floor(res / 12));
    this.ctx.strokeStyle = '#000';
    this.ctx.fillStyle = '#000';
    this.ctx.lineWidth = Math.max(2, cellW * 0.15);
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    const arrowShaftLen = cellW * 3.0;
    const arrowHeadLen = cellW * 1.0;
    for (let iy = 1; iy <= res; iy += stride) {
      for (let ix = 1; ix <= res; ix += stride) {
        const i = (iy - 1) * res + (ix - 1);
        if (solid[i] > 0.5) continue;
        const u = Ux[i];
        const v = Uy[i];
        const mag = Math.hypot(u, v) + 1e-9;
        const x = (ix - 1) * h - this.L / 2;
        const y = (iy - 1) * h - this.L / 2;
        const { cx, cy } = this.#physToCanvas(x, y);
        const dx = (u / mag) * arrowShaftLen;
        const dy = -(v / mag) * arrowShaftLen;
        const ex = cx + dx;
        const ey = cy + dy;
        this.ctx.beginPath();
        this.ctx.moveTo(cx, cy);
        this.ctx.lineTo(ex, ey);
        this.ctx.stroke();
        if (arrowShaftLen > arrowHeadLen * 0.5) {
          const shaftLen = Math.hypot(dx, dy) + 1e-9;
          const ux = dx / shaftLen;
          const uy = dy / shaftLen;
          const ax = ex - ux * arrowHeadLen;
          const ay = ey - uy * arrowHeadLen;
          const wing = arrowHeadLen * 0.4;
          const px = -uy * wing;
          const py = ux * wing;
          // Arrowhead as small triangle
          this.ctx.beginPath();
          this.ctx.moveTo(ex, ey);
          this.ctx.lineTo(ax + px, ay + py);
          this.ctx.lineTo(ax - px, ay - py);
          this.ctx.closePath();
          this.ctx.fill();
        }
      }
    }
    this.ctx.fillStyle = '#000';
    const particleRadius = Math.max(4, (this.L * scale) / res * 0.8);
    for (let i = 0; i < particles.x.length; i++) {
      const { cx, cy } = this.#physToCanvas(particles.x[i], particles.y[i]);
      this.ctx.beginPath();
      this.ctx.arc(cx, cy, particleRadius, 0, Math.PI * 2);
      this.ctx.fill();
    }
  }
}
