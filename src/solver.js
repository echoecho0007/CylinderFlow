/* ==========================================================================
   CylinderFlow — D2Q9 lattice Boltzmann solver (BGK collision)
   Pure numerical code: no DOM access, so the same source is unit-tested
   headlessly and inlined into index.html.
   Lattice units throughout: dx = dt = 1, reference density rho0 = 1.
   ========================================================================== */
const LBM = (function () {
  'use strict';
  // D2Q9 velocity set. Index: 0 rest, 1 E, 2 N, 3 W, 4 S, 5 NE, 6 NW, 7 SW, 8 SE
  const Q = 9;
  const CX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
  const CY = [0, 0, 1, 0, -1, 1, 1, -1, -1];
  const W = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
  const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];     // c_opp = -c
  const MIRY = [0, 1, 4, 3, 2, 8, 7, 6, 5];    // reflect c_y -> -c_y (free-slip wall)
  const CS2 = 1 / 3;                           // lattice speed of sound squared

  /** Equilibrium distribution, 2nd-order Hermite expansion. */
  function feq(i, rho, ux, uy) {
    const cu = CX[i] * ux + CY[i] * uy;
    const usq = ux * ux + uy * uy;
    return W[i] * rho * (1 + 3 * cu + 4.5 * cu * cu - 1.5 * usq);
  }

  /** nu = cs^2 (tau - 1/2) = tau/3 - 1/6 */
  const viscosityFromTau = (tau) => CS2 * (tau - 0.5);
  const tauFromViscosity = (nu) => nu / CS2 + 0.5;

  class Solver {
    /**
     * cfg: nx, ny, tau, U (inlet velocity, lattice), D (cylinder diameter, cells),
     *      cxc, cyc (cylinder centre), xBC: 'inout' | 'periodic',
     *      yBC: 'freeslip' | 'periodic', gx, gy (body acceleration, tests only),
     *      obstacle: 'cylinder' | 'none', solidFn(x,y) (optional extra solids),
     *      perturb (fractional transverse velocity seed, default 0)
     */
    constructor(cfg) {
      this.cfg = Object.assign({
        xBC: 'inout', yBC: 'freeslip', gx: 0, gy: 0, obstacle: 'cylinder', perturb: 0,
        sponge: 0, spongeTau: 1.0
      }, cfg);
      const c = this.cfg;
      this.nx = c.nx | 0; this.ny = c.ny | 0; this.N = this.nx * this.ny;
      this.tau = c.tau; this.omega = 1 / c.tau; this.nu = viscosityFromTau(c.tau);
      this.U = c.U; this.D = c.D;
      this.cxc = c.cxc !== undefined ? c.cxc : this.nx / 5;
      this.cyc = c.cyc !== undefined ? c.cyc : (this.ny - 1) / 2;
      const N = this.N;
      this.f = new Float64Array(Q * N);
      this.fn = new Float64Array(Q * N);
      this.rho = new Float64Array(N);
      this.ux = new Float64Array(N);
      this.uy = new Float64Array(N);
      this.solid = new Uint8Array(N);
      this.slow = new Uint8Array(N);   // node needs the boundary-aware gather path
      // per-column relaxation frequency: uniform, except an optional absorbing
      // ("sponge") layer next to the outlet where viscosity rises smoothly to damp
      // acoustic waves reflected by the simple outflow condition
      this.omegaX = new Float64Array(this.nx).fill(this.omega);
      if (c.xBC === 'inout' && c.sponge > 0) {
        const xs = Math.floor(this.nx * (1 - c.sponge)), tS = Math.max(c.spongeTau, c.tau);
        for (let x = xs; x < this.nx; x++) {
          const s = (x - xs) / (this.nx - 1 - xs);
          this.omegaX[x] = 1 / (c.tau + (tS - c.tau) * s * s);
        }
        this.spongeStart = xs;
      } else this.spongeStart = this.nx;
      this.t = 0;
      this.Fx = 0; this.Fy = 0;        // force ON the obstacle (momentum exchange), last step
      this.mass = 0; this.meanUx = 0; this.maxU = 0; this.minRho = 1; this.maxRho = 1;
      this.nFluid = 0;
      this.rhoIn = 1;                  // mean density on the first interior column (pressure reference)
      this.status = 'ok'; this.statusMsg = '';
      this.buildGeometry();
      this.initialise();
    }

    buildGeometry() {
      const { nx, ny, solid } = this, c = this.cfg;
      const R = this.D / 2;
      let nf = 0;
      for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
        const idx = y * nx + x;
        let s = 0;
        if (c.obstacle === 'cylinder') {
          const dx = x - this.cxc, dy = y - this.cyc;
          if (dx * dx + dy * dy <= R * R) s = 1;
        }
        if (c.solidFn && c.solidFn(x, y)) s = 1;
        solid[idx] = s;
        if (!s) nf++;
      }
      this.nFluid = nf;
      // nodes whose neighbours are solid or outside the domain use the slow path
      for (let y = 0; y < ny; y++) for (let x = 0; x < nx; x++) {
        const idx = y * nx + x;
        let s = (x === 0 || x === nx - 1 || y === 0 || y === ny - 1) ? 1 : 0;
        if (!s) for (let i = 1; i < Q; i++) if (solid[idx - CX[i] - CY[i] * nx]) { s = 1; break; }
        this.slow[idx] = s;
      }
      // boundary links (for Cp around the cylinder surface)
      this.surface = [];
      for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) {
        const idx = y * nx + x;
        if (solid[idx]) continue;
        for (let i = 1; i < Q; i++) if (solid[idx + CX[i] + CY[i] * nx]) {
          this.surface.push(idx); break;
        }
      }
    }

    initialise() {
      const { N, nx, f, solid, rho, ux, uy } = this;
      const c = this.cfg, U0 = c.xBC === 'inout' ? this.U : (c.U0 || 0);
      for (let idx = 0; idx < N; idx++) {
        const x = idx % nx;
        let u = solid[idx] ? 0 : U0, v = 0;
        // optional symmetry-breaking seed: one weak, divergence-free Gaussian
        // vortex 1.5 D behind the cylinder (identical for every Re). Peak swirl
        // velocity = perturb * U. It decays where the steady wake is stable and
        // triggers the natural instability where it is not.
        if (!solid[idx] && c.perturb) {
          const y = (idx / nx) | 0, r0 = 0.5 * this.D;
          const dx = x - (this.cxc + 1.5 * this.D), dy = y - this.cyc;
          const e = Math.exp(-(dx * dx + dy * dy) / (r0 * r0));
          const G = c.perturb * this.U * Math.sqrt(2 * Math.E) / r0;    // so that max |u'| = perturb*U
          u += -G * dy * e; v = G * dx * e;
        }
        if (c.initFn) { const r = c.initFn(x, (idx / nx) | 0); u = r[0]; v = r[1]; }
        rho[idx] = 1; ux[idx] = u; uy[idx] = v;
        for (let i = 0; i < Q; i++) f[i * N + idx] = feq(i, 1, u, v);
      }
      this.t = 0; this.status = 'ok'; this.statusMsg = '';
      this.Fx = this.Fy = 0;
      this.computeMass();
    }

    computeMass() {
      const { N, f, solid } = this; let m = 0;
      for (let idx = 0; idx < N; idx++) if (!solid[idx]) for (let i = 0; i < Q; i++) m += f[i * N + idx];
      this.mass = m; return m;
    }

    /** One time step: pull-streaming + boundary conditions + BGK collision (fused). */
    step() {
      if (this.status !== 'ok') return false;
      const { nx, ny, N, f, fn, rho: R, ux: UX, uy: UY, solid, slow } = this;
      const omX = this.omegaX, c = this.cfg;
      const gxs = this.tau * c.gx, gys = this.tau * c.gy;        // velocity-shift forcing
      const inout = c.xBC === 'inout', yper = c.yBC === 'periodic';
      const x0 = inout ? 1 : 0, x1 = inout ? nx - 2 : nx - 1;
      const N1 = N, N2 = 2 * N, N3 = 3 * N, N4 = 4 * N, N5 = 5 * N, N6 = 6 * N, N7 = 7 * N, N8 = 8 * N;
      let Fx = 0, Fy = 0, sumRho = 0, sumUx = 0, maxU2 = 0, minR = 1e9, maxR = -1e9;
      const g = [0, 0, 0, 0, 0, 0, 0, 0, 0];

      for (let y = 0; y < ny; y++) {
        for (let x = x0; x <= x1; x++) {
          const idx = y * nx + x;
          if (solid[idx]) continue;
          let f0, f1, f2, f3, f4, f5, f6, f7, f8;
          if (!slow[idx]) {
            // interior: f_i(x) <- f*_i(x - c_i)
            f0 = f[idx];
            f1 = f[N1 + idx - 1];
            f2 = f[N2 + idx - nx];
            f3 = f[N3 + idx + 1];
            f4 = f[N4 + idx + nx];
            f5 = f[N5 + idx - nx - 1];
            f6 = f[N6 + idx - nx + 1];
            f7 = f[N7 + idx + nx + 1];
            f8 = f[N8 + idx + nx - 1];
          } else {
            g[0] = f[idx];
            for (let i = 1; i < Q; i++) {
              let sx = x - CX[i], sy = y - CY[i];
              if (sx < 0 || sx >= nx) sx = (sx + nx) % nx;          // periodic in x (inout never reaches here)
              if (sy < 0 || sy >= ny) {
                if (yper) sy = (sy + ny) % ny;
                else {                                              // free-slip: specular reflection
                  const j = MIRY[i];
                  const sidx = y * nx + sx;
                  g[i] = solid[sidx] ? f[OPP[i] * N + idx] : f[j * N + sidx];
                  continue;
                }
              }
              const sidx = sy * nx + sx;
              if (solid[sidx]) {
                // half-way bounce-back; population that left x towards the wall returns reversed
                const v = f[OPP[i] * N + idx];
                g[i] = v;
                // momentum exchange: 2 c_opp f*_opp  (force on the solid)
                Fx += 2 * CX[OPP[i]] * v; Fy += 2 * CY[OPP[i]] * v;
              } else g[i] = f[i * N + sidx];
            }
            f0 = g[0]; f1 = g[1]; f2 = g[2]; f3 = g[3]; f4 = g[4]; f5 = g[5]; f6 = g[6]; f7 = g[7]; f8 = g[8];
          }
          // macroscopic moments
          const rho = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8;
          const ir = 1 / rho;
          const ux = (f1 - f3 + f5 - f6 - f7 + f8) * ir;
          const uy = (f2 - f4 + f5 + f6 - f7 - f8) * ir;
          // BGK collision towards equilibrium
          const om = omX[x];
          const ue = ux + gxs, ve = uy + gys;
          const usq = 1.5 * (ue * ue + ve * ve);
          const r9 = rho / 9, r36 = rho / 36;
          const cu1 = 3 * ue, cu2 = 3 * ve, cu5 = 3 * (ue + ve), cu6 = 3 * (ve - ue);
          fn[idx] = f0 + om * (4 * r9 * (1 - usq) - f0);
          fn[N1 + idx] = f1 + om * (r9 * (1 + cu1 + 0.5 * cu1 * cu1 - usq) - f1);
          fn[N2 + idx] = f2 + om * (r9 * (1 + cu2 + 0.5 * cu2 * cu2 - usq) - f2);
          fn[N3 + idx] = f3 + om * (r9 * (1 - cu1 + 0.5 * cu1 * cu1 - usq) - f3);
          fn[N4 + idx] = f4 + om * (r9 * (1 - cu2 + 0.5 * cu2 * cu2 - usq) - f4);
          fn[N5 + idx] = f5 + om * (r36 * (1 + cu5 + 0.5 * cu5 * cu5 - usq) - f5);
          fn[N6 + idx] = f6 + om * (r36 * (1 + cu6 + 0.5 * cu6 * cu6 - usq) - f6);
          fn[N7 + idx] = f7 + om * (r36 * (1 - cu5 + 0.5 * cu5 * cu5 - usq) - f7);
          fn[N8 + idx] = f8 + om * (r36 * (1 - cu6 + 0.5 * cu6 * cu6 - usq) - f8);
          // physical velocity with forcing = u + g/2 (zero for the cylinder runs)
          const uxp = ux + 0.5 * c.gx, uyp = uy + 0.5 * c.gy;
          R[idx] = rho; UX[idx] = uxp; UY[idx] = uyp;
          sumRho += rho; sumUx += uxp;
          const u2 = uxp * uxp + uyp * uyp;
          if (u2 > maxU2) maxU2 = u2;
          if (rho < minR) minR = rho;
          if (rho > maxR) maxR = rho;
        }
      }

      if (inout) {
        // inlet (x = 0): equilibrium with prescribed velocity, density from neighbour
        const U = this.U;
        let sr = 0, nr = 0;
        for (let y = 0; y < ny; y++) {
          const idx = y * nx;
          if (solid[idx]) continue;
          const r = solid[idx + 1] ? 1 : R[idx + 1];
          sr += r; nr++;
          for (let i = 0; i < Q; i++) fn[i * N + idx] = feq(i, r, U, 0);
          R[idx] = r; UX[idx] = U; UY[idx] = 0;
        }
        if (nr) this.rhoIn = sr / nr;
        // outlet (x = nx-1): zero-gradient (copy of the upstream column)
        for (let y = 0; y < ny; y++) {
          const idx = y * nx + nx - 1;
          for (let i = 0; i < Q; i++) fn[i * N + idx] = fn[i * N + idx - 1];
          R[idx] = R[idx - 1]; UX[idx] = UX[idx - 1]; UY[idx] = UY[idx - 1];
        }
      }

      this.f = fn; this.fn = f;
      this.t++;
      this.Fx = Fx; this.Fy = Fy;
      const nInt = this.nFluid - (inout ? 2 * ny : 0);
      this.meanRho = sumRho / nInt;
      this.meanUx = sumUx / nInt;
      this.maxU = Math.sqrt(maxU2);
      this.minRho = minR; this.maxRho = maxR;

      // ---- stability checks -------------------------------------------------
      if (!Number.isFinite(sumRho) || !Number.isFinite(Fx) || !Number.isFinite(Fy)) {
        this.status = 'unstable';
        this.statusMsg = 'Non-finite (NaN/Infinity) density or velocity detected at step ' + this.t + '.';
      } else if (minR <= 0.3 || maxR >= 3) {
        this.status = 'unstable';
        this.statusMsg = 'Density left the physically admissible range (min ' + minR.toFixed(3) +
          ', max ' + maxR.toFixed(3) + ') at step ' + this.t + '.';
      } else if (this.maxU > 0.45) {
        this.status = 'unstable';
        this.statusMsg = 'Local velocity ' + this.maxU.toFixed(3) + ' exceeded 0.45 lattice units ' +
          '(Mach ' + (this.maxU / Math.sqrt(CS2)).toFixed(2) + ') at step ' + this.t + '.';
      }
      return this.status === 'ok';
    }

    /** Drag/lift coefficients based on inlet velocity and diameter (per unit span). */
    coefficients() {
      const q = 0.5 * 1.0 * this.U * this.U * this.D;
      return { Cd: this.Fx / q, Cl: this.Fy / q };
    }

    /** Vorticity  w = dv/dx - du/dy  (central differences, y up). */
    vorticity(out) {
      const { nx, ny, ux, uy, solid } = this;
      for (let y = 1; y < ny - 1; y++) for (let x = 1; x < nx - 1; x++) {
        const idx = y * nx + x;
        out[idx] = solid[idx] ? 0 :
          0.5 * (uy[idx + 1] - uy[idx - 1]) - 0.5 * (ux[idx + nx] - ux[idx - nx]);
      }
      return out;
    }

    /** Recirculation length behind the cylinder on the centreline, in diameters (null if none). */
    recirculationLength() {
      const { nx, ux, solid } = this;
      const y = Math.round(this.cyc);
      let x = Math.floor(this.cxc + this.D / 2) + 1;
      while (x < nx && solid[y * nx + x]) x++;
      const xs = x;
      if (ux[y * nx + x] >= 0 && ux[y * nx + x + 1] >= 0) return null;
      while (x < nx - 2 && ux[y * nx + x] < 0) x++;
      if (x >= nx - 2) return null;
      const u0 = ux[y * nx + x - 1], u1 = ux[y * nx + x];
      const xr = (x - 1) + (-u0) / (u1 - u0);           // linear interpolation of u = 0
      const rear = this.cxc + this.D / 2;
      return { LrD: (xr - rear) / this.D, xReattach: xr, xStart: xs };
    }

    /** Pressure coefficient Cp = cs^2 (rho - rho_in) / (U^2/2), referenced to the mean inlet-plane pressure. */
    cp(idx) { return CS2 * (this.rho[idx] - this.rhoIn) / (0.5 * this.U * this.U); }

    /** Pressure coefficient along the cylinder surface (first fluid nodes), sorted by angle from the front stagnation point. */
    surfaceCp() {
      const { nx, rho } = this, q = 0.5 * this.U * this.U, out = [];
      for (const idx of this.surface) {
        const x = idx % nx, y = (idx / nx) | 0;
        let th = Math.atan2(y - this.cyc, x - this.cxc);          // 0 = rear
        let deg = 180 - Math.abs(th) * 180 / Math.PI;           // 0 = front stagnation
        out.push({ theta: deg, upper: th >= 0, Cp: CS2 * (rho[idx] - this.rhoIn) / q });
      }
      out.sort((a, b) => a.theta - b.theta);
      return out;
    }
  }

  /** Analyse a lift-coefficient history: amplitude, rms, and Strouhal number from zero crossings. */
  function analyseSignal(tStar, sig, fromFrac) {
    const n = sig.length;
    const i0 = Math.floor(n * (fromFrac === undefined ? 0.5 : fromFrac));
    if (n - i0 < 20) return null;
    let mean = 0, min = Infinity, max = -Infinity;
    for (let i = i0; i < n; i++) { mean += sig[i]; if (sig[i] < min) min = sig[i]; if (sig[i] > max) max = sig[i]; }
    mean /= (n - i0);
    let rms = 0;
    for (let i = i0; i < n; i++) rms += (sig[i] - mean) * (sig[i] - mean);
    rms = Math.sqrt(rms / (n - i0));
    // upward zero crossings of (sig - mean), linear interpolation
    const ups = [];
    for (let i = i0 + 1; i < n; i++) {
      const a = sig[i - 1] - mean, b = sig[i] - mean;
      if (a < 0 && b >= 0) ups.push(tStar[i - 1] + (tStar[i] - tStar[i - 1]) * (-a) / (b - a));
    }
    let period = null, periodCV = null;
    if (ups.length >= 3) {
      const ps = [];
      for (let k = 1; k < ups.length; k++) ps.push(ups[k] - ups[k - 1]);
      const pm = ps.reduce((s, v) => s + v, 0) / ps.length;
      const sd = Math.sqrt(ps.reduce((s, v) => s + (v - pm) * (v - pm), 0) / ps.length);
      period = pm; periodCV = sd / pm;
    }
    return { mean, min, max, amp: 0.5 * (max - min), rms, period, periodCV,
             St: period ? 1 / period : null, cycles: Math.max(0, ups.length - 1), window: [tStar[i0], tStar[n - 1]] };
  }

  /**
   * Summarise a run from its recorded history (t* = t U / D).
   * Statistics use the window  t* in [max(t_end/2, t_end - 40), t_end]
   * so the start-up transient is excluded. Classification thresholds are
   * documented in the Methodology panel.
   */
  function summarise(t, cd, cl) {
    const n = t.length;
    if (n < 40) return null;
    const tEnd = t[n - 1], tFrom = Math.max(0.5 * tEnd, tEnd - 40);
    let i0 = 0; while (i0 < n - 1 && t[i0] < tFrom) i0++;
    const a = analyseSignal(t, cl, i0 / n), b = analyseSignal(t, cd, i0 / n);
    if (!a || !b) return null;
    let state, stateText;
    if (tEnd < 15) { state = 'developing'; stateText = 'Developing (t* < 15): start-up transient'; }
    else if (a.amp > 0.02 && a.cycles >= 3 && a.periodCV < 0.1) { state = 'shedding'; stateText = 'Periodic vortex shedding'; }
    else if (a.amp < 0.01 && b.rms < 0.005 * Math.abs(b.mean)) { state = 'steady'; stateText = 'Steady, symmetric wake'; }
    else { state = 'transitional'; stateText = 'Unsteady, not yet periodic: keep running'; }
    return {
      window: [t[i0], tEnd], CdMean: b.mean, CdStd: b.rms, CdMin: b.min, CdMax: b.max,
      ClMean: a.mean, ClAmp: a.amp, ClRms: a.rms, St: state === 'shedding' ? a.St : null,
      cycles: a.cycles, periodCV: a.periodCV, state, stateText
    };
  }

  return { Q, CX, CY, W, OPP, MIRY, CS2, feq, viscosityFromTau, tauFromViscosity, Solver, analyseSignal, summarise };
})();
if (typeof module !== 'undefined') module.exports = LBM;
