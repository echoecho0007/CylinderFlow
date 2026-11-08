// Reproduce one row of results/reference_runs.jsonl:
//   node scripts/reference_run.js <Re> <nx> <ny> <D> <t*>   e.g.  node scripts/reference_run.js 100 400 160 20 100
const LBM=require('../src/solver.js');
const args=process.argv.slice(2).map(Number);
if(args.length<5){console.error('usage: node scripts/reference_run.js <Re> <nx> <ny> <D> <t*>');process.exit(1);}
const [Re,nx,ny,D,Tstar]=args; const U=0.08;
const s=new LBM.Solver({nx,ny,tau:LBM.tauFromViscosity(U*D/Re),U,D,cxc:Math.round(nx/5),cyc:(ny-1)/2,perturb:0.1,sponge:0.15,spongeTau:1});
const every=Math.max(1,Math.round(0.05*D/U)); const t=[],cd=[],cl=[];
const steps=Math.round(Tstar*D/U); const t0=Date.now();
for(let k=1;k<=steps;k++){ if(!s.step()){console.log('UNSTABLE');break;} if(s.t%every===0){const c=s.coefficients();t.push(s.t*U/D);cd.push(c.Cd);cl.push(c.Cl);} }
const S=LBM.summarise(t,cd,cl), lr=s.recirculationLength();
console.log(JSON.stringify({Re,grid:nx+'x'+ny,D,tau:+s.tau.toFixed(4),tEnd:+t[t.length-1].toFixed(1),CdMean:+S.CdMean.toFixed(3),ClAmp:+S.ClAmp.toFixed(4),St:S.St&&+S.St.toFixed(4),LrD:(S.state==='steady'&&lr)?+lr.LrD.toFixed(3):null,state:S.state,cycles:S.cycles,sec:(Date.now()-t0)/1000}));
