// Run the solver verification suite headlessly:  node scripts/run_tests.js
const LBM = require('../src/solver.js'), run = require('../src/tests.js');
const t0 = Date.now(), r = run(LBM, false);
for (const x of r) console.log(x.pass ? 'PASS' : 'FAIL', '|', x.group, '|', x.name, '|', x.detail);
console.log(`${r.filter(x => x.pass).length}/${r.length} passed in ${Date.now() - t0} ms`);
process.exit(r.every(x => x.pass) ? 0 : 1);
