/**
 * Obstacle contour via marching squares algorithm.
 * Cell layout: tl---tr (row r), bl---br (row r+1); columns c, c+1.
 *
 *       base = (tl) ------ (tr)  --- c
 *                  |           |
 *                  |           |
 *                  |           |
 *                 (bl) ------ (br)

 *                  |
 *                  r
 * 
 * Grid: Float32Array row-major I(iy,ix)=(iy-1)*res+(ix-1). We map (iy,ix) cell to (r,c) with r=iy-1, c=ix-1.
 * Output: segments as { p0: {x,y}, p1: {x,y} } in 0-based index space (x=col, y=row).
 * Smoothing is applied only to a copy of the mask (solver mask unchanged).
 *
 * @author Ilya Tsivilskiy
 */

const idx = (res) => (iy, ix) => (iy - 1) * res + (ix - 1);

/**
 * Inverse linear interpolation: coordinate at which the linear segment (rLow, fLow)–(rHigh, fHigh) equals f.
 * @param {number} rLow - Lower coordinate
 * @param {number} rHigh - Higher coordinate
 * @param {number} fLow - Field value at rLow
 * @param {number} fHigh - Field value at rHigh
 * @param {number} f - Target value
 * @returns {number} Interpolated coordinate where the segment attains f
 */
function linterp(rLow, rHigh, fLow, fHigh, f) {
  return rHigh - (rHigh - rLow) * (fHigh - f) / (fHigh - fLow + Number.EPSILON);
}

/**
 * One pass of 3×3 box filter on interior nodes; boundary values copied from source. Does not mutate input.
 * @param {Float32Array} a - Scalar field, row-major, length res²
 * @param {number} res - Grid resolution (rows/cols)
 * @returns {Float32Array} New array with smoothed interior and original boundary
 */
function smoothCopy(a, res) {
  const I = idx(res);
  const out = new Float32Array(a.length);
  for (let i = 0; i < a.length; i++) out[i] = a[i];
  for (let iy = 2; iy < res; iy++) {
    for (let ix = 2; ix < res; ix++) {
      let sum = 0;
      for (let dy = -1; dy <= 1; dy++)
        for (let dx = -1; dx <= 1; dx++) sum += a[I(iy + dy, ix + dx)];
      out[I(iy, ix)] = sum / 9;
    }
  }
  return out;
}

class ObstacleContour {
  /**
   * Build contour line patches at iso-value. Uses a copy of the mask; optionally smooths that copy.
   *
   * @param {Float32Array} solid - Solid mask (1=solid, 0=fluid), row-major, length res²
   * @param {number} res - Grid resolution (rows/cols)
   * @param {number} [isoVal=0.5] - Iso-value
   * @param {boolean} [smooth=true] - If true, smooth a copy for nicer contours (solver unchanged)
   * @param {number} [smoothPasses=1] - Smoothing passes when smooth is true
   * @returns {{ p0: {x: number, y: number}, p1: {x: number, y: number} }[]}
   */
  static build(solid, res, isoVal = 0.5, smooth = true, smoothPasses = 1) {
    let aC = solid;
    if (smooth && smoothPasses > 0) {
      aC = new Float32Array(solid);
      for (let p = 0; p < smoothPasses; p++) aC = smoothCopy(aC, res);
    }

    const I = idx(res);
    const linePatches = [];

    // Contour2D: for r=0..rows-1, c=0..cols-1. We use iy=r+1, ix=c+1 so r=iy-1, c=ix-1.
    for (let iy = 1; iy < res; iy++) {
      for (let ix = 1; ix < res; ix++) {
        const r = iy - 1;
        const c = ix - 1;
        const tl = { x: c, y: r, f: aC[I(iy, ix)] };
        const tr = { x: c + 1, y: r, f: aC[I(iy, ix + 1)] };
        const bl = { x: c, y: r + 1, f: aC[I(iy + 1, ix)] };
        const br = { x: c + 1, y: r + 1, f: aC[I(iy + 1, ix + 1)] };

        let p0 = null;
        let p1 = null;

        // case 1: all below
        if (tl.f < isoVal && tr.f < isoVal && bl.f < isoVal && br.f < isoVal) {
          // no segment
        }
        // case 2: bl above
        if (tl.f < isoVal && tr.f < isoVal && bl.f > isoVal && br.f < isoVal) {
          p0 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
          p1 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
        }
        // case 3: br above
        if (tl.f < isoVal && tr.f < isoVal && bl.f < isoVal && br.f > isoVal) {
          p0 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
          p1 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
        }
        // case 4: bl, br above
        if (tl.f < isoVal && tr.f < isoVal && bl.f > isoVal && br.f > isoVal) {
          p0 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
          p1 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
        }
        // case 5: tr above
        if (tl.f < isoVal && tr.f > isoVal && bl.f < isoVal && br.f < isoVal) {
          p0 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
          p1 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
        }
        // case 6: tr, bl above (saddle - skip like reference)
        if (tl.f < isoVal && tr.f > isoVal && bl.f > isoVal && br.f < isoVal) {
          // need 2 points - not implemented in reference
        }
        // case 7: tr, br above
        if (tl.f < isoVal && tr.f > isoVal && bl.f < isoVal && br.f > isoVal) {
          p0 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
          p1 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
        }
        // case 8: tr, bl, br above
        if (tl.f < isoVal && tr.f > isoVal && bl.f > isoVal && br.f > isoVal) {
          p0 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
          p1 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
        }
        // case 9: tl above
        if (tl.f > isoVal && tr.f < isoVal && bl.f < isoVal && br.f < isoVal) {
          p0 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
          p1 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
        }
        // case 10: tl, bl above
        if (tl.f > isoVal && tr.f < isoVal && bl.f > isoVal && br.f < isoVal) {
          p0 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
          p1 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
        }
        // case 11: tl, br above (saddle - skip like reference)
        if (tl.f > isoVal && tr.f < isoVal && bl.f < isoVal && br.f > isoVal) {
          // need 2 points
        }
        // case 12: tl, bl, br above
        if (tl.f > isoVal && tr.f < isoVal && bl.f > isoVal && br.f > isoVal) {
          p0 = { x: linterp(tl.x, tr.x, tl.f, tr.f, isoVal), y: tl.y };
          p1 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
        }
        // case 13: tl, tr above
        if (tl.f > isoVal && tr.f > isoVal && bl.f < isoVal && br.f < isoVal) {
          p0 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
          p1 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
        }
        // case 14: tl, tr, bl above
        if (tl.f > isoVal && tr.f > isoVal && bl.f > isoVal && br.f < isoVal) {
          p0 = { x: br.x, y: linterp(br.y, tr.y, br.f, tr.f, isoVal) };
          p1 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
        }
        // case 15: tl, tr, br above
        if (tl.f > isoVal && tr.f > isoVal && bl.f < isoVal && br.f > isoVal) {
          p0 = { x: linterp(bl.x, br.x, bl.f, br.f, isoVal), y: bl.y };
          p1 = { x: bl.x, y: linterp(bl.y, tl.y, bl.f, tl.f, isoVal) };
        }
        // case 16: all above
        if (tl.f > isoVal && tr.f > isoVal && bl.f > isoVal && br.f > isoVal) {
          // no segment
        }

        if (p0 !== null && p1 !== null) {
          linePatches.push({ p0, p1 });
        }
      }
    }

    return linePatches;
  }
}
