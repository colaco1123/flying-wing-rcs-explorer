eval(require('fs').readFileSync('core3d.js','utf8'));
function mk(v){const [a,b,c]=v;const u=[b[0]-a[0],b[1]-a[1],b[2]-a[2]],w=[c[0]-a[0],c[1]-a[1],c[2]-a[2]];
 const n=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]];const L=Math.hypot(...n);
 return Float64Array.from([n[0]/L,n[1]/L,n[2]/L,...a,...b,...c,L/2]);}
const v=[[0.3,-0.2,0.5],[1.1,0.4,0.2],[-0.1,0.9,1.0]];const T=mk(v);
for(const q of [[5,-3,2],[0.2,7,-4],[12,1,9]]){
 const out=[0,0];triFT(T,0,...q,out);
 let re=0,im=0,N=1200;for(let i=0;i<N;i++)for(let j=0;j<N-i;j++){const u=(i+1/3)/N,w=(j+1/3)/N;
  const r=[0,1,2].map(k=>v[0][k]+u*(v[1][k]-v[0][k])+w*(v[2][k]-v[0][k]));const ph=q[0]*r[0]+q[1]*r[1]+q[2]*r[2];re+=Math.cos(ph);im+=Math.sin(ph);}
 const W=2*T[12]/(N*N);console.log(out.map(x=>x.toFixed(5)).join(' '),'|',(re*W).toFixed(5),(im*W).toFixed(5));}
// specular
const n=[T[0],T[1],T[2]];const out=[0,0];triFT(T,0,n[0]*4,n[1]*4,n[2]*4,out);console.log(Math.hypot(...out).toFixed(6),T[12].toFixed(6));
