/**
 * Shared math utilities for the fluid solver and graphics.
 *
 * @author Ilya Tsivilskiy
 */
const MathUtils = {
  /** Small epsilon for safe division and float comparisons. */
  EPS: 1e-12,

  /**
   * Clamp value to [lo, hi].
   * @param {number} f - Value
   * @param {number} lo - Lower bound
   * @param {number} hi - Upper bound
   * @returns {number}
   */
  clamp(f, lo, hi) {
    return Math.max(lo, Math.min(hi, f));
  },

  /**
   * Linear interpolation: value at p given f1 at p1 and f2 at p2.
   * @param {number} f1 - Value at p1
   * @param {number} f2 - Value at p2
   * @param {number} p1 - Coordinate left/bottom
   * @param {number} p - Query coordinate
   * @param {number} p2 - Coordinate right/top
   * @returns {number} Interpolated value
   */
  lerp(f1, f2, p1, p, p2) {
    const len = p2 - p1;
    if (Math.abs(len) < this.EPS) return (f1 + f2) * 0.5;
    return f1 * ((p2 - p) / len) + f2 * ((p - p1) / len);
  },

  /**
   * Inverse linear interpolation: coordinate where the field equals iso
   * between (rLow, fLow) and (rHigh, fHigh). Returns midpoint if denominator is negligible; result clamped to [rLow, rHigh].
   * @param {number} rLow - Lower coordinate
   * @param {number} rHigh - Higher coordinate
   * @param {number} fLow - Field value at rLow
   * @param {number} fHigh - Field value at rHigh
   * @param {number} iso - Target value
   * @returns {number} Interpolated coordinate in [rLow, rHigh]
   */
  lerpInv(rLow, rHigh, fLow, fHigh, iso) {
    const denom = fHigh - fLow;
    if (Math.abs(denom) < this.EPS) return (rLow + rHigh) * 0.5;
    const t = (iso - fLow) / denom;
    const r = rLow + t * (rHigh - rLow);
    return Math.max(rLow, Math.min(rHigh, r));
  },

  /**
   * Check if value is in closed interval [lo, hi].
   * @param {number} a - Value
   * @param {number} lo - Lower bound
   * @param {number} hi - Upper bound
   * @returns {boolean}
   */
  inRange(a, lo, hi) {
    return lo <= a && a <= hi;
  }
};
