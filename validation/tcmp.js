const {add,sub,mul,div,sc,cexp,fact,cpow}=require('./pof.js');
var P=require('./pof.js');eval(require('fs').readFileSync('./core3d.js','utf8'));
const src=require('fs').readFileSync('tIc.js','utf8');eval(src.slice(src.indexOf('function IcPOF'),src.indexOf('function brute')));eval(src.slice(src.indexOf('function mk'),src.indexOf('const dot')));
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
let seed=7;const rnd=()=>((seed=(seed*16807)%2147483647)/2147483647)-0.5;
console.log('phase span   worst |edge-sinc - POF|/|I|   worst |POF(openrcsG) - POF(stable)|/|I|');
for(const span of [1e-7,1e-5,1e-3,1e-1,1,10,100]){
 let we=0,wo=0;
 for(let t=0;t<3000;t++){
  const v=[0,1,2].map(()=>[rnd(),rnd(),rnd()]);const q0=[rnd(),rnd(),rnd()];const L=Math.hypot(...q0);
  // scale q so in-triangle phase differences ~ span; also sometimes put q near-normal
  const T=mk(v);let q=q0.map(x=>x/L);if(t%3===0){q=[T[0]+q[0]*0.01,T[1]+q[1]*0.01,T[2]+q[2]*0.01];}
  const Dp0=dot(q,v[0].map((x,k)=>x-v[2][k])),Dq0=dot(q,v[1].map((x,k)=>x-v[2][k]));const s=span/Math.max(Math.abs(Dp0),Math.abs(Dq0),1e-30);
  q=q.map(x=>x*s);
  const Dp=Dp0*s,Dq=Dq0*s,Do=dot(q,v[2]);
  const a=IcPOF(Dp,Dq,Do,T[12],P.G_stable),b=IcPOF(Dp,Dq,Do,T[12],P.G_openrcs);const o=[0,0];triFT(T,0,...q,o);
  const m=Math.hypot(...a);we=Math.max(we,Math.hypot(o[0]-a[0],o[1]-a[1])/m);wo=Math.max(wo,Math.hypot(b[0]-a[0],b[1]-a[1])/m);
 }
 console.log(String(span).padEnd(12),we.toExponential(1).padStart(10),wo.toExponential(1).padStart(30));
}
