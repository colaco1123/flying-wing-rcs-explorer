const P=require('./pof.js');const {add,sub,mul,div,sc,cexp,fact,cpow}=P;
eval(require('fs').readFileSync('./core3d.js','utf8'));
function IcPOF(Dp,Dq,Do,A,G,Lt=1e-5,Nt=5){
 const DD=Dq-Dp,eDo=cexp(Do),eDp=cexp(Dp),eDq=cexp(Dq);let s=[0,0];
 if(Math.abs(Dp)<Lt&&Math.abs(Dq)>=Lt){for(let n=0;n<=Nt;n++)s=add(s,mul(sc(cpow([0,Dp],n),1/fact(n)),add([-1/(n+1),0],mul(eDq,G(n,-Dq)))));return sc(div(mul(s,eDo),[0,Dq]),2*A);}
 if(Math.abs(Dp)<Lt&&Math.abs(Dq)<Lt){for(let n=0;n<=Nt;n++)for(let m=0;m<Nt;m++)s=add(s,sc(mul(cpow([0,Dp],n),cpow([0,Dq],m)),1/fact(n+m+2)));return sc(mul(s,eDo),2*A);}
 if(Math.abs(Dp)>=Lt&&Math.abs(Dq)<Lt){for(let n=0;n<=Nt;n++)s=add(s,sc(mul(cpow([0,Dq],n),G(n+1,-Dp)),1/(fact(n)*(n+1))));return sc(mul(mul(s,eDo),eDp),2*A);}
 if(Math.abs(DD)<Lt){for(let n=0;n<=Nt;n++)s=add(s,mul(sc(cpow([0,DD],n),1/fact(n)),add(sc(G(n,Dq),-1),sc(eDq,1/(n+1)))));return sc(div(mul(s,eDo),[0,Dq]),2*A);}
 return sc(mul(eDo,sub(sub(div(eDp,[Dp*DD,0]),div(eDq,[Dq*DD,0])),[1/(Dp*Dq),0])),2*A);
}
function brute(v,q,N=1500){let re=0,im=0;for(let i=0;i<N;i++)for(let j=0;j<N-i;j++){const u=(i+1/3)/N,w=(j+1/3)/N;let ph=0;for(let k=0;k<3;k++)ph+=q[k]*(v[0][k]+u*(v[1][k]-v[0][k])+w*(v[2][k]-v[0][k]));re+=Math.cos(ph);im+=Math.sin(ph);}
 const a=[0,1,2].map(k=>v[1][k]-v[0][k]),b=[0,1,2].map(k=>v[2][k]-v[0][k]);const A=Math.hypot(a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])/2;return [re*2*A/(N*N),im*2*A/(N*N)];}
function mk(v){const [a,b,c]=v;const u=[b[0]-a[0],b[1]-a[1],b[2]-a[2]],w=[c[0]-a[0],c[1]-a[1],c[2]-a[2]];const n=[u[1]*w[2]-u[2]*w[1],u[2]*w[0]-u[0]*w[2],u[0]*w[1]-u[1]*w[0]];const L=Math.hypot(...n);return Float64Array.from([n[0]/L,n[1]/L,n[2]/L,...a,...b,...c,L/2]);}
const dot=(a,b)=>a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
const cases=[
 ['generic',[[0.3,-0.2,0.5],[1.1,0.4,0.2],[-0.1,0.9,1.0]],[5,-3,2]],
 ['small Dp (edge 1-3 ⟂ q)',[[0,0,0.002],[0.01,0.003,0],[0,0,0]],[0,40,3]],
 ['both small',[[1e-7,0,0],[0,1e-7,0],[0,0,0]],[30,20,10]],
 ['Dq small, Dp mid',[[0.02,0,0],[0,0.00001,0.0004],[0,0,0]],[200,10,0]],
 ['near specular',[[0,0,0],[0.01,0,0],[0,0.01,0]],[1e-6,-2e-6,400]],
 ['DD small',[[0.01,0.0,0],[0.01,0.00000001,0.001],[0,0,0]],[300,7,0]],
];
for(const [name,v,q] of cases){
 const A=mk(v)[12];const Dp=dot(q,v[0].map((x,k)=>x-v[2][k])),Dq=dot(q,v[1].map((x,k)=>x-v[2][k])),Do=dot(q,v[2]);
 const t=brute(v,q),rel=g=>(Math.hypot(g[0]-t[0],g[1]-t[1])/Math.hypot(...t)).toExponential(1);
 const o=[0,0];triFT(mk(v),0,...q,o);
 console.log(name.padEnd(26),'POF+openrcsG',rel(IcPOF(Dp,Dq,Do,A,P.G_openrcs)).padStart(8),' POF+stableG',rel(IcPOF(Dp,Dq,Do,A,P.G_stable)).padStart(8),' edge-sinc',rel(o).padStart(8));
}
console.log('\n--- direct comparison on inputs that land in each special-case branch (reference: edge-sinc in quad-ish via series check) ---');
// series reference: exact integral via 2D Gauss-Legendre high order on triangle (Duffy), accurate to ~1e-13 for smooth phase
function gl(n){const x=[],w=[];for(let i=1;i<=n;i++){let z=Math.cos(Math.PI*(i-0.25)/(n+0.5)),pp;for(let it=0;it<100;it++){let p1=1,p2=0;for(let j=1;j<=n;j++){const p3=p2;p2=p1;p1=((2*j-1)*z*p2-(j-1)*p3)/j;}pp=n*(z*p1-p2)/(z*z-1);const z1=z;z=z1-p1/pp;if(Math.abs(z-z1)<1e-15)break;}x.push(z);w.push(2/((1-z*z)*pp*pp));}return {x,w};}
const Q=gl(60);
function ref(Dp,Dq,Do,A){let re=0,im=0;for(let i=0;i<60;i++)for(let j=0;j<60;j++){const s=(Q.x[i]+1)/2,t=(Q.x[j]+1)/2;const u=s,v=(1-s)*t;const jac=(1-s)/4*Q.w[i]*Q.w[j];const ph=Do+Dp*u+Dq*v;re+=jac*Math.cos(ph);im+=jac*Math.sin(ph);}return [re*2*A,im*2*A];}
for(const [nm,Dp,Dq] of [['case1 |Dp|<Lt, Dq=2e-4',3e-6,2e-4],['case1 |Dp|<Lt, Dq=0.02',3e-6,0.02],['case1 |Dp|<Lt, Dq=1.5',3e-6,1.5],['case3 |Dq|<Lt, Dp=3e-3',3e-3,4e-6],['case3 |Dq|<Lt, Dp=0.3',0.3,4e-6],['case4 DD<Lt, Dq=5e-4',5e-4+2e-6,5e-4],['case4 DD<Lt, Dq=2',2+3e-6,2],['general',0.7,-2.1]]){
 const A=0.37,Do=1.3,t=ref(Dp,Dq,Do,A),rel=g=>(Math.hypot(g[0]-t[0],g[1]-t[1])/Math.hypot(...t)).toExponential(1);
 // build triangle with these phases for edge-sinc: r3=0, choose q=(1,0,0)... set r1=(Dp,a,0), r2=(Dq,b,0), r3=(Do,0,0)? phases = x-coord with q=x-hat; need area A
 const r3=[Do,0,0],r1=[Do+Dp,0.5,0],r2=[Do+Dq,0,0.5];const T=mk([r1,r2,r3]);const s=A/T[12];const T2=mk([r1.map((x,k)=>k?x*s:x),r2.map((x,k)=>k?x*s:x),r3]);
 const o=[0,0];triFT(T2,0,1,0,0,o);
 console.log(nm.padEnd(28),'openrcsG',rel(IcPOF(Dp,Dq,Do,A,P.G_openrcs)).padStart(8),' stableG',rel(IcPOF(Dp,Dq,Do,A,P.G_stable)).padStart(8),' edge-sinc',rel(o).padStart(8), ' A-check',(T2[12]/A).toFixed(6));
}
