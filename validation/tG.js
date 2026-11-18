const P=require('./pof.js');
console.log('n  w        |true|      openrcs_err   fixedrec_err   stable_err');
for(const w of [-1e-4,-1e-3,-0.05,-0.8,-3]) for(const n of [1,2,3,5]){
 const t=P.G_true(n,w),e=g=>Math.hypot(g[0]-t[0],g[1]-t[1])/Math.hypot(...t);
 console.log(n,String(w).padEnd(8),Math.hypot(...t).toExponential(3),e(P.G_openrcs(n,w)).toExponential(1).padStart(12),e(P.G_rec(n,w)).toExponential(1).padStart(13),e(P.G_stable(n,w)).toExponential(1).padStart(12));}
