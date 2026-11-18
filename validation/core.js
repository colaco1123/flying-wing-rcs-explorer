function polyArea(P){let a=0;for(let i=0;i<P.length;i++){const p=P[i],q=P[(i+1)%P.length];a+=p[0]*q[1]-q[0]*p[1];}return a/2;}
function polyFT(P,px,py,scale){
  const p2=px*px+py*py;
  if(p2*scale*scale<1e-10) return [polyArea(P),0];
  let tr=0,ti=0;
  for(let i=0;i<P.length;i++){
    const a=P[i],b=P[(i+1)%P.length];
    const ax=b[0]-a[0],ay=b[1]-a[1];
    if(ax===0&&ay===0) continue;
    const cx=(a[0]+b[0])/2,cy=(a[1]+b[1])/2;
    const s=px*ay-py*ax, u=(px*ax+py*ay)/2;
    const sc=Math.abs(u)<1e-9?1:Math.sin(u)/u;
    const arg=px*cx+py*cy, m=s*sc;
    tr+=m*Math.cos(arg); ti+=m*Math.sin(arg);
  }
  return [ti/p2,-tr/p2];
}
function rcs(P,lambda,el,phi,scale){
  const k=2*Math.PI/lambda, c=2*k*Math.cos(el);
  const F=polyFT(P,c*Math.cos(phi),c*Math.sin(phi),scale);
  const s=Math.sin(el);
  return 4*Math.PI/(lambda*lambda)*s*s*(F[0]*F[0]+F[1]*F[1]);
}
