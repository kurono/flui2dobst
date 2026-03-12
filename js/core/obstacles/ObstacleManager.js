/**
 * Obstacle geometry: draw/erase circles on the solid mask, smooth, binarize, and expose solid normals.
 *
 * Solid mask: 1 = solid, 0 = fluid. Used for free-slip BC and particle collision.
 * Smoothing improves surface normals.
 *
 * @author Ilya Tsivilskiy
 * @see FluidSolver.idx, FluidSolver.gridToPhys, FluidSolver.physToGrid, FluidSolver.solidNormal
 */
class ObstacleManager {
  /**
   * Set solid mask to 1 inside the circle centered at (px, py) in physical coords.
   * @param {Float32Array} solid - Solid mask (modified in place)
   * @param {number} px - Center x [m]
   * @param {number} py - Center y [m]
   * @param {number} radius - Radius [m]
   * @param {number} res - Grid resolution
   * @param {number} L - Domain size
   * @param {number} h - Cell size
   */
  static drawCircle(solid, px, py, radius, res, L, h) {
    const I = FluidSolver.idx(res);
    const rCells = Math.max(1, Math.ceil(radius / h));
    const { fiy: cy, fix: cx } = FluidSolver.physToGrid(py, px, h, L);
    const iy0 = Math.max(1, Math.floor(cy - rCells));
    const iy1 = Math.min(res, Math.ceil(cy + rCells));
    const ix0 = Math.max(1, Math.floor(cx - rCells));
    const ix1 = Math.min(res, Math.ceil(cx + rCells));
    for (let iy = iy0; iy <= iy1; iy++) {
      for (let ix = ix0; ix <= ix1; ix++) {
        const { y, x } = FluidSolver.gridToPhys(iy, ix, h, L);
        if (Math.hypot(x - px, y - py) <= radius) solid[I(iy, ix)] = 1;
      }
    }
  }

  /**
   * Set solid mask to 0 inside the circle (erase obstacle).
   * @param {Float32Array} solid - Solid mask (modified in place)
   * @param {number} px - Center x [m]
   * @param {number} py - Center y [m]
   * @param {number} radius - Radius [m]
   * @param {number} res - Grid resolution
   * @param {number} L - Domain size
   * @param {number} h - Cell size
   */
  static eraseCircle(solid, px, py, radius, res, L, h) {
    const I = FluidSolver.idx(res);
    const rCells = Math.max(1, Math.ceil(radius / h));
    const { fiy: cy, fix: cx } = FluidSolver.physToGrid(py, px, h, L);
    const iy0 = Math.max(1, Math.floor(cy - rCells));
    const iy1 = Math.min(res, Math.ceil(cy + rCells));
    const ix0 = Math.max(1, Math.floor(cx - rCells));
    const ix1 = Math.min(res, Math.ceil(cx + rCells));
    for (let iy = iy0; iy <= iy1; iy++) {
      for (let ix = ix0; ix <= ix1; ix++) {
        const { y, x } = FluidSolver.gridToPhys(iy, ix, h, L);
        if (Math.hypot(x - px, y - py) <= radius) solid[I(iy, ix)] = 0;
      }
    }
  }

  /**
   * Smooth solid mask with a 3×3 box filter (multiple passes). Improves surface normals.
   * 9-point average for interior, keep boundaries unchanged.
   * @param {Float32Array} solid - Solid mask (modified in place)
   * @param {number} res - Grid resolution
   * @param {number} [passes=1] - Number of smoothing passes
   */
  static smooth(solid, res, passes = 1) {
    const I = FluidSolver.idx(res);
    let cur = solid;
    for (let p = 0; p < passes; p++) {
      const next = new Float32Array(res * res);
      for (let iy = 2; iy < res; iy++) {
        for (let ix = 2; ix < res; ix++) {
          let acc = 0;
          for (let oy = -1; oy <= 1; oy++) {
            for (let ox = -1; ox <= 1; ox++) acc += cur[I(iy + oy, ix + ox)];
          }
          next[I(iy, ix)] = acc / 9;
        }
      }
      for (let iy = 1; iy <= res; iy++) {
        next[I(iy, 1)] = cur[I(iy, 1)];
        next[I(iy, res)] = cur[I(iy, res)];
      }
      for (let ix = 1; ix <= res; ix++) {
        next[I(1, ix)] = cur[I(1, ix)];
        next[I(res, ix)] = cur[I(res, ix)];
      }
      cur = next;
    }
    if (passes > 0) {
      for (let i = 0; i < res * res; i++) solid[i] = cur[i];
    }
  }

  /**
   * Binarize solid mask: values > threshold become 1, else 0.
   * @param {Float32Array} solid - Solid mask (modified in place)
   * @param {number} res - Grid resolution
   * @param {number} [threshold=0.5] - Cutoff
   */
  static binarize(solid, res, threshold = 0.5) {
    for (let i = 0; i < res * res; i++) solid[i] = solid[i] > threshold ? 1 : 0;
  }

  /**
   * Get outer unit normals (and tangents) for the solid surface. Delegates to FluidSolver.solidNormal.
   * @param {Float32Array} solid - Solid mask
   * @param {number} h - Cell size
   * @param {number} res - Resolution
   * @returns {{ snx: Float32Array, sny: Float32Array }}
   */
  static getSolidNormals(solid, h, res) {
    return FluidSolver.solidNormal(solid, h, res);
  }
}
