/**
 * Utility helpers (e.g. FPS meter) for the fluid app.
 * @author Ilya Tsivilskiy
 */

/**
 * FPS meter: rolling window of frame timestamps to compute current frames per second.
 */
class FPSMeter {
  /**
   * @param {number} [windowSize=30] - Number of recent frames to average over
   */
  constructor(windowSize = 30) {
    this.windowSize = windowSize;
    /** @type {number[]} */
    this.timestamps = [];
  }

  /**
   * Call each frame with the current high-res time (e.g. from requestAnimationFrame).
   * @param {number} t - Timestamp in milliseconds
   */
  tick(t) {
    this.timestamps.push(t);
    if (this.timestamps.length > this.windowSize) {
      this.timestamps.shift();
    }
  }

  /**
   * Current FPS over the last window, or 0 if not enough samples.
   * @returns {number}
   */
  getFPS() {
    if (this.timestamps.length < 2) return 0;
    const span = (this.timestamps[this.timestamps.length - 1] - this.timestamps[0]) / 1000;
    if (span <= 0) return 0;
    return (this.timestamps.length - 1) / span;
  }
}
