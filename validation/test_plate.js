eval(require('fs').readFileSync('core.js','utf8'));
// rectangle a x b centered: analytic F = ab sinc(px a/2) sinc(py b/2)
const a=1.2,b=0.7,P=[[-a/2,-b/2],[a/2,-b/2],[a/2,b/2],[-a/2,b/2]];
const sinc=x=>x?Math.sin(x)/x:1;
for(const [px,py] of [[3,1],[0.5,7],[40,-13],[0,0.001]]){
 const F=polyFT(P,px,py,1);const an=a*b*sinc(px*a/2)*sinc(py*b/2);
 console.log(Math.hypot(...F).toFixed(6), Math.abs(an).toFixed(6));}
// triangle vs brute force
const T=[[0,0],[-1.5,1],[-1,-0.3]];
let re=0,im=0,N=1500;const px=9,py=-4;
for(let i=0;i<N;i++)for(let j=0;j<N-i;j++){const u=(i+1/3)/N,v=(j+1/3)/N;const x=T[0][0]+u*(T[1][0]-T[0][0])+v*(T[2][0]-T[0][0]),y=T[0][1]+u*(T[1][1]-T[0][1])+v*(T[2][1]-T[0][1]);re+=Math.cos(px*x+py*y);im+=Math.sin(px*x+py*y);}
const A=Math.abs(polyArea(T));const w=2*A/(N*N);
console.log(Math.hypot(re*w,im*w).toFixed(5),Math.hypot(...polyFT(T,px,py,1)).toFixed(5));
// normal incidence square: 4pi a^4/lambda^2
const Q=[[0,0],[1,0],[1,1],[0,1]];console.log(rcs(Q,0.03,Math.PI/2,0,1),4*Math.PI/0.0009);
