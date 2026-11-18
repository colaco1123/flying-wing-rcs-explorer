// complex helpers
const C=(r,i=0)=>[r,i],add=(a,b)=>[a[0]+b[0],a[1]+b[1]],sub=(a,b)=>[a[0]-b[0],a[1]-b[1]],
mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]],div=(a,b)=>{const d=b[0]*b[0]+b[1]*b[1];return [(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d];},
sc=(a,s)=>[a[0]*s,a[1]*s],cexp=t=>[Math.cos(t),Math.sin(t)],fact=n=>{let f=1;for(let i=2;i<=n;i++)f*=i;return f;},
cpow=(a,n)=>{let r=[1,0];for(let i=0;i<n;i++)r=mul(r,a);return r;};
// G(n,w)=∫0^1 x^n e^{jwx}dx
function G_openrcs(n,w){const jw=[0,w];let g=div(sub(cexp(w),[1,0]),jw);for(let m=1;m<=n;m++){g=div(sub(cexp(w),sc(g,n)),jw);}return g;}
function G_rec(n,w){const jw=[0,w];let g=div(sub(cexp(w),[1,0]),jw);for(let m=1;m<=n;m++){g=div(sub(cexp(w),sc(g,m)),jw);}return g;}
function G_stable(n,w){ if(Math.abs(w)>=1) return G_rec(n,w); let s=[0,0],t=[1,0];for(let k=0;k<40;k++){s=add(s,sc(t,1/(n+k+1)));t=sc(mul(t,[0,w]),1/(k+1));}return s;}
function G_true(n,w){let re=0,im=0,N=200000;for(let i=0;i<N;i++){const x=(i+.5)/N;const xn=Math.pow(x,n);re+=xn*Math.cos(w*x);im+=xn*Math.sin(w*x);}return [re/N,im/N];}
module.exports={G_openrcs,G_rec,G_stable,G_true,C,add,sub,mul,div,sc,cexp,fact,cpow};
