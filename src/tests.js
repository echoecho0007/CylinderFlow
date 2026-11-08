/* ==========================================================================
   CylinderFlow — solver verification suite.
   Runs in Node during development and from the "Verification" tab in the app.
   Every check computes its own pass/fail from numbers; nothing is hard-coded.
   ========================================================================== */
function runSelfTests(LBM, quick) {
  'use strict';
  const { Q, CX, CY, W, OPP, MIRY, feq, Solver } = LBM;
  const results = [];
  const add = (group, name, pass, detail) => results.push({ group, name, pass: !!pass, detail });
  const close = (a, b, tol) => Math.abs(a - b) <= tol;
  let seed = 12345;
  const rnd = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };

  // ---------------------------------------------------------------- 1. lattice
  {
    let s = 0, sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sx4 = 0, sx2y2 = 0;
    for (let i = 0; i < Q; i++) {
      s += W[i]; sx += W[i] * CX[i]; sy += W[i] * CY[i];
      sxx += W[i] * CX[i] * CX[i]; syy += W[i] * CY[i] * CY[i]; sxy += W[i] * CX[i] * CY[i];
      sx4 += W[i] * CX[i] ** 4; sx2y2 += W[i] * CX[i] ** 2 * CY[i] ** 2;
    }
    add('Lattice', 'Weights sum to 1', close(s, 1, 1e-15), 'Σw = ' + s.toPrecision(16));
    add('Lattice', 'First moment Σ w c = 0', close(sx, 0, 1e-15) && close(sy, 0, 1e-15), `(${sx}, ${sy})`);
    add('Lattice', 'Second moment Σ w cαcβ = cs² δαβ', close(sxx, 1 / 3, 1e-15) && close(syy, 1 / 3, 1e-15) && close(sxy, 0, 1e-15),
      `xx=${sxx.toFixed(15)}, yy=${syy.toFixed(15)}, xy=${sxy}`);
    add('Lattice', 'Fourth moment isotropy (Σw cx⁴ = 3cs⁴, Σw cx²cy² = cs⁴)', close(sx4, 1 / 3, 1e-15) && close(sx2y2, 1 / 9, 1e-15),
      `${sx4.toFixed(15)}, ${sx2y2.toFixed(15)}`);
    let okOpp = true, okMir = true, okSet = true;
    const seen = new Set();
    for (let i = 0; i < Q; i++) {
      if (CX[OPP[i]] !== -CX[i] || CY[OPP[i]] !== -CY[i]) okOpp = false;
      if (CX[MIRY[i]] !== CX[i] || CY[MIRY[i]] !== -CY[i]) okMir = false;
      seen.add(CX[i] + ',' + CY[i]);
      if (Math.abs(CX[i]) > 1 || Math.abs(CY[i]) > 1) okSet = false;
    }
    add('Lattice', 'Nine distinct velocities in {-1,0,1}²', okSet && seen.size === 9, seen.size + ' distinct');
    add('Lattice', 'Opposite-direction table (bounce-back)', okOpp, 'c_opp(i) = −c_i for all i');
    add('Lattice', 'Mirror table (free-slip walls)', okMir, 'c_mir(i) = (c_x, −c_y)');
  }

  // ---------------------------------------------------------------- 2. equilibrium
  {
    let worst = 0;
    for (let k = 0; k < 200; k++) {
      const rho = 0.8 + 0.4 * rnd(), ux = 0.3 * (rnd() - 0.5), uy = 0.3 * (rnd() - 0.5);
      let m0 = 0, mx = 0, my = 0, pxx = 0, pyy = 0, pxy = 0;
      for (let i = 0; i < Q; i++) {
        const fe = feq(i, rho, ux, uy);
        m0 += fe; mx += fe * CX[i]; my += fe * CY[i];
        pxx += fe * CX[i] * CX[i]; pyy += fe * CY[i] * CY[i]; pxy += fe * CX[i] * CY[i];
      }
      worst = Math.max(worst, Math.abs(m0 - rho), Math.abs(mx - rho * ux), Math.abs(my - rho * uy),
        Math.abs(pxx - (rho * ux * ux + rho / 3)), Math.abs(pyy - (rho * uy * uy + rho / 3)), Math.abs(pxy - rho * ux * uy));
    }
    add('Equilibrium', 'Moments of f_eq: ρ, ρu, ρuu + ρcs²I (200 random states)', worst < 1e-14, 'max error ' + worst.toExponential(2));
    let pos = true;
    for (let i = 0; i < Q; i++) if (feq(i, 1, 0.1, 0.05) <= 0) pos = false;
    add('Equilibrium', 'f_eq positive at |u| ≈ 0.11', pos, 'all nine populations > 0');
  }

  // ---------------------------------------------------------------- 3. streaming
  {
    const s = new Solver({ nx: 6, ny: 6, tau: 1, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none' });
    const N = s.N;
    for (let idx = 0; idx < N; idx++) for (let i = 0; i < Q; i++) s.f[i * N + idx] = W[i];
    let ok = true;
    const marks = [];
    for (let i = 1; i < Q; i++) { const idx = 2 * 6 + 2; s.f[i * N + idx] += 0.01 * i; marks.push(i); }
    s.omegaX.fill(0);                             // pure streaming
    s.step();
    for (const i of marks) {
      const dst = (2 + CY[i]) * 6 + (2 + CX[i]);
      if (!close(s.f[i * N + dst] - W[i], 0.01 * i, 1e-15)) ok = false;
    }
    // periodic wrap: a population at x = 5 moving east arrives at x = 0
    const s2 = new Solver({ nx: 6, ny: 6, tau: 1, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none' });
    for (let idx = 0; idx < N; idx++) for (let i = 0; i < Q; i++) s2.f[i * N + idx] = W[i];
    s2.f[1 * N + 3 * 6 + 5] += 0.5; s2.omegaX.fill(0); s2.step();
    const wrap = close(s2.f[1 * N + 3 * 6 + 0] - W[1], 0.5, 1e-15);
    add('Streaming', 'Each population moves exactly one link along c_i', ok, '8 marked populations checked');
    add('Streaming', 'Periodic wrap-around', wrap, 'east-moving population x=5 → x=0');
  }

  // ---------------------------------------------------------------- 4. collision conserves ρ and ρu
  {
    const s = new Solver({ nx: 8, ny: 8, tau: 0.7, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none' });
    const N = s.N;
    for (let idx = 0; idx < N; idx++) for (let i = 0; i < Q; i++) s.f[i * N + idx] = W[i] * (0.9 + 0.2 * rnd());
    const mom = (f) => { let m = 0, px = 0, py = 0; for (let idx = 0; idx < N; idx++) for (let i = 0; i < Q; i++) { const v = f[i * N + idx]; m += v; px += v * CX[i]; py += v * CY[i]; } return [m, px, py]; };
    const a = mom(s.f);
    for (let k = 0; k < 50; k++) s.step();
    const b = mom(s.f);
    const err = Math.max(Math.abs(a[0] - b[0]) / a[0], Math.abs(a[1] - b[1]), Math.abs(a[2] - b[2]));
    add('Collision', 'BGK collision + streaming conserve mass and momentum (random non-equilibrium state, 50 steps)', err < 1e-12, 'max drift ' + err.toExponential(2));
  }

  // ---------------------------------------------------------------- 5. viscosity: decaying shear wave
  for (const tau of [0.6, 0.8, 1.2]) {
    const L = 64, A = 1e-3, k = 2 * Math.PI / L;
    const s = new Solver({ nx: 4, ny: L, tau, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none',
      initFn: (x, y) => [A * Math.sin(k * y), 0] });
    const amp = () => { let a = 0; for (let y = 0; y < L; y++) a += s.ux[y * 4] * Math.sin(k * y); return 2 * a / L; };
    // velocity arrays are filled after the first step; measure decay between steps 10 and 10+T
    for (let n = 0; n < 10; n++) s.step();
    const a0 = amp(), T = quick ? 400 : 1000;
    for (let n = 0; n < T; n++) s.step();
    const a1 = amp();
    const nuMeas = -Math.log(a1 / a0) / (k * k * T);
    const nuTh = (tau - 0.5) / 3;
    const rel = Math.abs(nuMeas - nuTh) / nuTh;
    add('Viscosity', `Shear-wave decay gives ν = cs²(τ−½) (τ = ${tau})`, rel < 0.01,
      `measured ${nuMeas.toExponential(4)}, theory ${nuTh.toExponential(4)}, error ${(100 * rel).toFixed(3)} %`);
  }

  // ---------------------------------------------------------------- 6. bounce-back: Poiseuille flow
  {
    const ny = 34, tau = 0.8, nu = (tau - 0.5) / 3, gx = 1e-6;
    const s = new Solver({ nx: 3, ny, tau, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none', gx,
      solidFn: (x, y) => y === 0 || y === ny - 1 });
    const steps = quick ? 6000 : 12000;
    for (let n = 0; n < steps; n++) s.step();
    // walls located half-way between the solid row and the first fluid row: y = 0.5 and y = ny − 1.5
    const H = ny - 2;
    let num = 0, den = 0, umax = 0;
    for (let y = 1; y < ny - 1; y++) {
      const yy = y - 0.5;
      const ua = gx / (2 * nu) * yy * (H - yy);
      const un = s.ux[y * 3 + 1];
      num += (un - ua) ** 2; den += ua * ua; umax = Math.max(umax, un);
    }
    const l2 = Math.sqrt(num / den);
    add('Boundaries', 'Half-way bounce-back reproduces Poiseuille profile (no-slip walls)', l2 < 0.01,
      `L2 error ${(100 * l2).toFixed(3)} %, u_max = ${umax.toExponential(3)} (theory ${(gx * H * H / (8 * nu)).toExponential(3)})`);
  }

  // ---------------------------------------------------------------- 7. mass conservation in a closed box
  {
    const n = 24;
    const s = new Solver({ nx: n, ny: n, tau: 0.6, U: 0, D: 0, xBC: 'periodic', yBC: 'periodic', obstacle: 'none',
      solidFn: (x, y) => x === 0 || y === 0 || x === n - 1 || y === n - 1 || ((x - 12) ** 2 + (y - 12) ** 2 < 9),
      initFn: (x, y) => [0.05 * Math.sin(x / 3), 0.05 * Math.cos(y / 4)] });
    const m0 = s.computeMass();
    for (let k = 0; k < 2000; k++) s.step();
    const m1 = s.computeMass();
    add('Boundaries', 'Closed box with bounce-back walls + obstacle conserves mass (2000 steps)', Math.abs(m1 - m0) / m0 < 1e-12,
      `relative drift ${(Math.abs(m1 - m0) / m0).toExponential(2)}`);
  }

  // ---------------------------------------------------------------- 8. momentum exchange = fluid momentum change
  {
    const s = new Solver({ nx: 60, ny: 40, tau: 0.65, U: 0.05, D: 10, cxc: 20, cyc: 19.5, xBC: 'periodic', yBC: 'periodic', U0: 0.05 });
    const N = s.N;
    const P = () => { let px = 0, py = 0; for (let idx = 0; idx < N; idx++) if (!s.solid[idx]) for (let i = 0; i < Q; i++) { px += s.f[i * N + idx] * CX[i]; py += s.f[i * N + idx] * CY[i]; } return [px, py]; };
    let worst = 0, Fmax = 0;
    for (let k = 0; k < 300; k++) {
      const a = P(); s.step(); const b = P();
      worst = Math.max(worst, Math.abs((b[0] - a[0]) + s.Fx), Math.abs((b[1] - a[1]) + s.Fy));
      Fmax = Math.max(Fmax, Math.abs(s.Fx));
    }
    add('Forces', 'Momentum-exchange force equals the fluid momentum lost each step', worst < 1e-12 * Math.max(1, Fmax * 1e3),
      `max |ΔP + F| = ${worst.toExponential(2)} (|F| up to ${Fmax.toExponential(2)})`);
  }

  // ---------------------------------------------------------------- 9. inlet/outlet consistency
  {
    const U = 0.08;
    const s = new Solver({ nx: 80, ny: 30, tau: 0.6, U, D: 0, obstacle: 'none' });
    for (let k = 0; k < 500; k++) s.step();
    let e = 0;
    for (let idx = 0; idx < s.N; idx++) e = Math.max(e, Math.abs(s.ux[idx] - U), Math.abs(s.uy[idx]), Math.abs(s.rho[idx] - 1));
    add('Boundaries', 'Empty channel: inlet/outlet/free-slip preserve uniform flow', e < 1e-12, 'max deviation ' + e.toExponential(2));
  }

  // ---------------------------------------------------------------- 10. Reynolds-number mapping
  {
    let worst = 0;
    for (const [Re, U, D] of [[20, 0.08, 20], [100, 0.08, 20], [180, 0.1, 24], [40, 0.05, 16]]) {
      const nu = U * D / Re, tau = LBM.tauFromViscosity(nu);
      const back = U * D / LBM.viscosityFromTau(tau);
      worst = Math.max(worst, Math.abs(back - Re) / Re, Math.abs(tau - (3 * nu + 0.5)));
    }
    add('Parameters', 'Re = UD/ν and τ = 3ν + ½ round-trip', worst < 1e-12, 'max relative error ' + worst.toExponential(2));
  }

  // ---------------------------------------------------------------- 11. symmetry of the cylinder flow
  {
    const s = new Solver({ nx: 160, ny: 64, tau: 0.5 + 3 * (0.08 * 12 / 20), U: 0.08, D: 12, cxc: 32, cyc: 31.5, perturb: 0 });
    let maxCl = 0;
    for (let k = 0; k < (quick ? 800 : 2000); k++) { s.step(); maxCl = Math.max(maxCl, Math.abs(s.coefficients().Cl)); }
    add('Forces', 'Symmetric set-up without seed gives zero lift (Re 20)', maxCl < 1e-9, 'max |Cl| = ' + maxCl.toExponential(2));
    add('Forces', 'Drag is positive and finite', s.coefficients().Cd > 0 && Number.isFinite(s.coefficients().Cd), 'Cd = ' + s.coefficients().Cd.toFixed(4));
    let solidU = 0;
    for (let idx = 0; idx < s.N; idx++) if (s.solid[idx]) solidU = Math.max(solidU, Math.abs(s.ux[idx]), Math.abs(s.uy[idx]));
    add('Boundaries', 'Velocity inside the cylinder stays zero', solidU === 0, 'max |u| in solid = ' + solidU);
  }

  // ---------------------------------------------------------------- 12. instability detection
  {
    const s = new Solver({ nx: 120, ny: 50, tau: 0.5003, U: 0.25, D: 12, perturb: 0.1 });
    let k = 0; while (k < 5000 && s.step()) k++;
    add('Stability', 'Divergent run is detected and stopped (τ = 0.5003, U = 0.25)', s.status === 'unstable',
      s.status === 'unstable' ? s.statusMsg : 'not detected within 5000 steps');
  }
  return results;
}
if (typeof module !== 'undefined') module.exports = runSelfTests;
