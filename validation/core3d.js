// T: Float64Array per tri 13 values: nx,ny,nz, x0,y0,z0, x1,y1,z1, x2,y2,z2, area
function triFT(T,o,qx,qy,qz,out){
  const nx=T[o],ny=T[o+1],nz=T[o+2];
  const qn=qx*nx+qy*ny+qz*nz, p2=qx*qx+qy*qy+qz*qz-qn*qn;
  const A=T[o+12];
  const qq=qx*qx+qy*qy+qz*qz;
  if(p2<1e-10*qq+1e-14){ // specular: constant in-plane phase
    const cx=(T[o+3]+T[o+6]+T[o+9])/3,cy=(T[o+4]+T[o+7]+T[o+10])/3,cz=(T[o+5]+T[o+8]+T[o+11])/3;
    const ar=qx*cx+qy*cy+qz*cz; out[0]=A*Math.cos(ar); out[1]=A*Math.sin(ar); return;
  }
  let tr=0,ti=0;
  for(let e=0;e<3;e++){
    const a=o+3+3*e,b=o+3+3*((e+1)%3);
    const ax=T[b]-T[a],ay=T[b+1]-T[a+1],az=T[b+2]-T[a+2];
    // s = n·(q×a)
    const s=nx*(qy*az-qz*ay)+ny*(qz*ax-qx*az)+nz*(qx*ay-qy*ax);
    const u=(qx*ax+qy*ay+qz*az)/2;
    const sc=Math.abs(u)<1e-9?1:Math.sin(u)/u;
    const arg=qx*(T[a]+T[b])/2+qy*(T[a+1]+T[b+1])/2+qz*(T[a+2]+T[b+2])/2;
    const m=s*sc; tr+=m*Math.cos(arg); ti+=m*Math.sin(arg);
  }
  out[0]=ti/p2; out[1]=-tr/p2;
}
