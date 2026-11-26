<script>
(function(){
// FIGURE MODE. ?fig=depth hides the chrome and fixes the depth axis at 70 m;
// the selected gate model becomes the node nearest the CENTRE of the admissible
// range, which is what the paper's panel shows. Nothing else changes.
const FIGMODE = new URLSearchParams(location.search).get('fig');
if(FIGMODE){
  document.body.classList.add('figmode');
  // headless Chrome reports a LIGHT colour scheme, so the dark palette has to
  // be asked for rather than inherited: a figure is always on the dark ground.
  document.documentElement.setAttribute('data-theme','dark');
  const zq = new URLSearchParams(location.search).get('zmax');
  if(zq) D.zmax = parseFloat(zq);
}
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const F = (x,n) => x.toFixed(n===undefined?3:n);
const clamp=(a,b,c)=>Math.max(b,Math.min(c,a));

// ---- depth grids -> cell edges
const aemTop = D.aem.top.slice();
const aemBot = D.aem.bot.slice();
const smrEdge = (()=>{const z=D.snmr.z, e=[Math.max(0,z[0]-(z[1]-z[0])/2)];
  for(let i=0;i<z.length-1;i++) e.push((z[i]+z[i+1])/2);
  e.push(z[z.length-1]+(z[z.length-1]-z[z.length-2])/2); return e;})();

const PAD={l:42,r:10,t:10,b:30};

// ---- inferno, the matplotlib ramp, as 9 stops linearly interpolated.
// Density is drawn as an image on its own dark ground rather than as ink on
// the page, so the panel reads the same in either theme and the overlay lines
// have one background to be legible against instead of two.
const INFERNO=[[0,0,4],[31,12,72],[85,15,109],[136,34,106],[186,54,85],
               [227,89,51],[249,142,9],[249,212,45],[252,255,164]];
function ramp(t){
  t=Math.max(0,Math.min(1,t));
  const u=t*(INFERNO.length-1), i=Math.min(INFERNO.length-2,Math.floor(u)), f=u-i;
  const a=INFERNO[i], b=INFERNO[i+1];
  return 'rgb('+Math.round(a[0]+(b[0]-a[0])*f)+','+Math.round(a[1]+(b[1]-a[1])*f)+','
               +Math.round(a[2]+(b[2]-a[2])*f)+')';
}
const GROUND=ramp(0);
// Inferno runs black -> purple -> red -> orange -> near-white yellow, so the
// overlays are cold hues the ramp never reaches, and the truth carries a dark
// halo: the p50 sits exactly where density peaks, which is the palest part of
// the ramp, and a plain white line would disappear into it.
//   white  the true model, in every panel
  //   cyan   the ENSEMBLE -- p50 and the p05/p95 pair, and nothing else
  //   green  the selected gate model. It had been cyan too, which read as a
  //          third percentile line: one model forced to a conductance is not
  //          a posterior summary and must not share the ensemble's colour.
//   amber  the lowest-misfit model, a reference and nothing else
//   green  a cell resting on a stated limit
const C_P50='#00e5ff', C_PTL='rgba(0,229,255,.72)', C_TRUE='#ffffff',
      C_FAINT='rgba(255,255,255,.22)',
        // The admissible family is CONTEXT, not lines to read one by one. 17
        // strokes at 20% opacity accumulate to near-white where they overlap,
        // which is why it still competed with the truth; 12% and a thinner
        // line keep the eye on the truth and the selected model.
        C_FAM='rgba(255,255,255,.12)', C_BEST='#ffb703',
      C_SEL='#4ade80', C_LIM='#ffffff', C_GRID='rgba(255,255,255,.13)',
      // magenta: a node where no smoother model fits, so the min-chi2 model
      // stands. Kept well away from amber, which is the lowest-misfit
      // reference and would otherwise be confused with it in this panel.
      C_FB='#f25fd0', C_FB_FAINT='rgba(242,95,208,.26)',
      C_HALO='rgba(0,0,0,.6)';

function setup(cv){
  const r=cv.getBoundingClientRect(), d=window.devicePixelRatio||1;
  cv.width=Math.max(1,Math.round(r.width*d)); cv.height=Math.max(1,Math.round(r.height*d));
  const g=cv.getContext('2d'); g.setTransform(d,0,0,d,0,0);
  g.clearRect(0,0,r.width,r.height);
  return {g,w:r.width,h:r.height};
}
// axis furniture sits on the PAGE and is themed; the plot interior is the
// plasma ground and everything drawn inside it is light
function frame(g,w,h,xlab,xticks,xof,ylab){
  const x0=PAD.l,x1=w-PAD.r,y0=PAD.t,y1=h-PAD.b;
  g.fillStyle=GROUND; g.fillRect(x0,y0,x1-x0,y1-y0);
  g.font='9.5px ui-monospace,monospace';
  g.strokeStyle=C_GRID; g.lineWidth=1;
  for(const t of xticks){const X=xof(t); if(X<x0+1||X>x1-1) continue;
    g.beginPath(); g.moveTo(X,y0); g.lineTo(X,y1); g.stroke();}
  const step = D.zmax<=30?5:(D.zmax<=60?10:20);
  for(let d=step; d<D.zmax; d+=step){
    const Y=y0+(d/D.zmax)*(y1-y0);
    g.beginPath(); g.moveTo(x0,Y); g.lineTo(x1,Y); g.stroke();}
  g.fillStyle=css('--mute'); g.textAlign='center'; g.textBaseline='top';
  for(const t of xticks){const X=xof(t); if(X<x0-1||X>x1+1) continue;
    g.fillText(String(t.lab!==undefined?t.lab:t.v),X,y1+4);}
  g.textAlign='right'; g.textBaseline='middle';
  for(let d=0; d<=D.zmax; d+=step){
    const Y=y0+(d/D.zmax)*(y1-y0); g.fillText(String(d),x0-5,Y);}
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText(ylab,0,0); g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y0,x1-x0,y1-y0);
  g.fillStyle=css('--mute'); g.textAlign='right'; g.textBaseline='top';
  g.fillText(xlab,x1,y1+15);
  return {x0,x1,y0,y1};
}
const yOf=(A,d)=>A.y0+(clamp(d,0,D.zmax)/D.zmax)*(A.y1-A.y0);
function band(g,A){
  if(!D.win || D.win[0]===null || !isFinite(D.win[0])) return;   // FIELD: no interval detected
  const a=yOf(A,D.win[0]), b=yOf(A,D.win[1]);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  // the tint competes with the density underneath in print; the two dashed
  // lines carry the window on their own
  if(!FIGMODE){ g.fillStyle='rgba(255,255,255,.07)'; g.fillRect(A.x0,a,A.x1-A.x0,b-a); }
  g.strokeStyle='rgba(255,255,255,.55)'; g.lineWidth=1; g.setLineDash([4,3]);
  g.beginPath(); g.moveTo(A.x0,a); g.lineTo(A.x1,a); g.moveTo(A.x0,b); g.lineTo(A.x1,b);
  g.stroke(); g.setLineDash([]); g.restore();
}
// Below the sensitivity depth the models are drawn, then dimmed: they are
// still what the solver returned, but they are not a statement about the
// ground, so they should not compete with the part that is.
function veilAt(g,A,z,label){
  if(z===null||z===undefined||z>=D.zmax) return;
  const Y=yOf(A,z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.fillStyle='rgba(0,0,0,.52)'; g.fillRect(A.x0,Y,A.x1-A.x0,A.y1-Y);
  g.strokeStyle='rgba(255,255,255,.65)'; g.lineWidth=1.2; g.setLineDash([2,3]);
  g.beginPath(); g.moveTo(A.x0,Y); g.lineTo(A.x1,Y); g.stroke(); g.setLineDash([]);
  g.fillStyle='rgba(255,255,255,.8)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom';
  g.fillText(label+' '+z.toFixed(0)+' m', A.x0+4, Y-2);
  g.restore();
}
function note(g,A,txt){
  g.save(); g.fillStyle='rgba(255,255,255,.62)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom'; g.fillText(txt,A.x0+4,A.y1-4); g.restore();
}
function veil(g,A){
  if(!D.sens || D.sens.z>=D.zmax) return;
  const Y=yOf(A,D.sens.z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.fillStyle='rgba(0,0,0,.52)'; g.fillRect(A.x0,Y,A.x1-A.x0,A.y1-Y);
  g.restore();
}
function sensline(g,A){
  if(!D.sens || D.sens.z>=D.zmax) return;
  const Y=yOf(A,D.sens.z);
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.strokeStyle='rgba(255,255,255,.65)'; g.lineWidth=1.2; g.setLineDash([2,3]);
  g.beginPath(); g.moveTo(A.x0,Y); g.lineTo(A.x1,Y); g.stroke(); g.setLineDash([]);
  g.fillStyle='rgba(255,255,255,.8)'; g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='bottom';
  g.fillText('sensitivity depth '+D.sens.z.toFixed(0)+' m', A.x0+4, Y-2);
  g.restore();
}
// piecewise-constant: for a layered model, which is what it is
function steps(g,A,xof,v,top,bot,col,wd,dash,halo,alpha){
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  if(alpha!==undefined) g.globalAlpha=alpha;
  const path=()=>{g.beginPath(); let started=false;
    for(let i=0;i<v.length;i++){
      if(top[i]>D.zmax) break;
      if(!isFinite(v[i])) continue;
      const X=xof(v[i]), ya=yOf(A,top[i]), yb=yOf(A,Math.min(bot[i],D.zmax));
      if(!started){g.moveTo(X,ya);started=true;} else {g.lineTo(X,ya);}
      g.lineTo(X,yb);}};
  if(halo){g.strokeStyle=C_HALO; g.lineWidth=wd+2.2; g.setLineDash([]); path(); g.stroke();}
  g.strokeStyle=col; g.lineWidth=wd; g.setLineDash(dash||[]); path(); g.stroke();
  g.setLineDash([]); g.restore();
}
// smooth: a Catmull-Rom through the per-cell values at their cell centres.
// The percentiles are a summary of a distribution that varies continuously
// with depth, not a layered model, so they are drawn as a curve. The truth
// and the gate's models stay as steps, because those ARE layered models.
function smooth(g,A,xof,v,top,bot,col,wd,dash,halo){
  const pts=[];
  for(let i=0;i<v.length;i++){
    const zc=(top[i]+bot[i])/2;
    if(!isFinite(v[i])) continue;
    pts.push({x:xof(v[i]),y:yOf(A,zc)});
    if(zc>D.zmax) break;
  }
  if(pts.length<2) return;
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  g.lineJoin='round';
  // MONOTONE between points: the cardinal tangents are kept, but every control
  // point is clamped into the box spanned by the two points it lies between, so
  // the curve cannot overshoot either of them.
  const cl=(q,a,b)=>Math.max(Math.min(a,b),Math.min(Math.max(a,b),q));
  const trace=()=>{
    g.beginPath(); g.moveTo(pts[0].x,pts[0].y);
    for(let i=0;i<pts.length-1;i++){
      const p0=pts[i>0?i-1:0], p1=pts[i], p2=pts[i+1],
            p3=pts[i+2<pts.length?i+2:pts.length-1];
      g.bezierCurveTo(cl(p1.x+(p2.x-p0.x)/6,p1.x,p2.x), cl(p1.y+(p2.y-p0.y)/6,p1.y,p2.y),
                      cl(p2.x-(p3.x-p1.x)/6,p1.x,p2.x), cl(p2.y-(p3.y-p1.y)/6,p1.y,p2.y),
                      p2.x, p2.y);
    }
  };
  // a dark halo first, as steps() does, so a white line stays legible on the ramp
  if(halo){ g.strokeStyle=C_HALO; g.lineWidth=wd+2.2; g.setLineDash([]); trace(); g.stroke(); }
  g.strokeStyle=col; g.lineWidth=wd; g.setLineDash(dash||[]); trace(); g.stroke();
  g.setLineDash([]); g.restore();
}
// per-depth normalisation: each row is the marginal at that depth scaled to
// its own maximum, so a well-determined depth and a poorly determined one are
// both readable. A single global scale would wash out everything below the
// depth of investigation.
function hist(g,A,H,rmax,lo,hi,nb,top,bot){
  const bw=(hi-lo)/nb;
  g.save(); g.beginPath(); g.rect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0); g.clip();
  for(let i=0;i<H.length;i++){
    if(top[i]>D.zmax) break;
    const ya=yOf(A,top[i]), yb=yOf(A,Math.min(bot[i],D.zmax)), hgt=Math.max(1,yb-ya);
    const row=H[i], mx=rmax[i];
    for(let j=0;j<nb;j++){
      const c=row[j]; if(!c) continue;
      const xa=A.x0+((j*bw))/(hi-lo)*(A.x1-A.x0);
      const xb=A.x0+(((j+1)*bw))/(hi-lo)*(A.x1-A.x0);
      g.fillStyle=ramp(Math.sqrt(c/mx));
      g.fillRect(xa,ya,Math.max(1,xb-xa)+0.4,hgt+0.4);
    }
  }
  g.restore();
}

const LOGT=[{v:-3,lab:'.001'},{v:-2,lab:'.01'},{v:-1,lab:'0.1'},{v:0,lab:'1'}];
const sxLog = A => (s => A.x0+((Math.log10(s)-D.aem.lo)/(D.aem.hi-D.aem.lo))*(A.x1-A.x0));

function drawAEM(){
  const {g,w,h}=setup(document.getElementById('c1')); if(w<8) return;
  const xof=t=>PAD.l+((t.v-D.aem.lo)/(D.aem.hi-D.aem.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'conductivity S/m',LOGT,xof,'depth m');
  const X=sxLog(A);
  hist(g,A,D.aem.h,D.aem.rmax,D.aem.lo,D.aem.hi,D.aem.nb,aemTop,aemBot);
  band(g,A);
  smooth(g,A,X,D.aem.p05,aemTop,aemBot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.aem.p95,aemTop,aemBot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.aem.truth,aemTop,aemBot,C_TRUE,2,null,true);
  // p50 ON TOP of the truth: where the two coincide the question is whether the inversion found it, so the inversion's line is the one that must be visible.
  smooth(g,A,X,D.aem.p50,aemTop,aemBot,C_P50,2.1);
}
function drawSNMR(){
  const {g,w,h}=setup(document.getElementById('c2')); if(w<8) return;
  const T=[{v:0,lab:'0'},{v:0.1,lab:'0.1'},{v:0.2,lab:'0.2'},{v:0.3,lab:'0.3'},{v:0.4,lab:'0.4'}];
  const xof=t=>PAD.l+((t.v-D.snmr.lo)/(D.snmr.hi-D.snmr.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'water content v/v',T,xof,'depth m');
  const X=s=>A.x0+((s-D.snmr.lo)/(D.snmr.hi-D.snmr.lo))*(A.x1-A.x0);
  const top=smrEdge.slice(0,-1), bot=smrEdge.slice(1);
  hist(g,A,D.snmr.h,D.snmr.rmax,D.snmr.lo,D.snmr.hi,D.snmr.nb,top,bot);
  band(g,A);
  smooth(g,A,X,D.snmr.p05,top,bot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.snmr.p95,top,bot,C_PTL,1.3,[5,3]);
  smooth(g,A,X,D.snmr.truth,top,bot,C_TRUE,2,null,true);
  // p50 ON TOP of the truth: where the two coincide the question is whether the inversion found it, so the inversion's line is the one that must be visible.
  smooth(g,A,X,D.snmr.p50,top,bot,C_P50,2.1);
  if(!FIGMODE && (D.snmr.dlim===null||D.snmr.dlim===undefined)) note(g,A,'depth limit not stored');
  else veilAt(g,A,D.snmr.dlim,'depth limit');
}
function drawFit(){
  const {g,w,h}=setup(document.getElementById('c2')); if(w<8) return;
  const F=D.fit;
  if(!F){ const A=frame(g,w,h,'',[],()=>-99,'');
    note(g,A,'SNMR data file not found'); return; }
  const qs=F.q, all=F.V.concat(F.sig).concat(F.haspred?F.p05.concat(F.p95):[]).filter(v=>v>0);
  const vlo=Math.log10(Math.min(...all))-0.08, vhi=Math.log10(Math.max(...all))+0.08;
  const qlo=Math.log10(Math.min(...qs))-0.05, qhi=Math.log10(Math.max(...qs))+0.05;
  const T=[]; for(let d=Math.ceil(qlo);d<=Math.floor(qhi);d++)
    T.push({v:d,lab:(d<0?'1e'+d:(d<4?String(Math.pow(10,d)):'1e'+d))});
  const xof=t=>PAD.l+((t.v-qlo)/(qhi-qlo))*(w-PAD.l-PAD.r);
  // Amplitude spans three decades while the noise is 0.07-2.5% of signal, so
  // a single log axis shows the sounding curve and hides the fit completely.
  // The panel is therefore split: the curve above, and below it the residuals
  // in units of the noise, which is where a 5-sigma misfit is actually visible.
  const x0=PAD.l,x1=w-PAD.r,y0=PAD.t,yB=h-PAD.b, SPLIT=0.60;
  const y1=y0+(yB-y0)*SPLIT-8, y2=y0+(yB-y0)*SPLIT+8;
  const X=q=>x0+((Math.log10(q)-qlo)/(qhi-qlo))*(x1-x0);
  const Y=v=>y1-((Math.log10(Math.max(v,1e-30))-vlo)/(vhi-vlo))*(y1-y0);
  g.fillStyle=GROUND; g.fillRect(x0,y0,x1-x0,y1-y0); g.fillRect(x0,y2,x1-x0,yB-y2);
  g.font='9.5px ui-monospace,monospace';
  g.strokeStyle=C_GRID; g.lineWidth=1;
  for(const t of T){const px=xof(t); if(px<x0+1||px>x1-1)continue;
    g.beginPath();g.moveTo(px,y0);g.lineTo(px,y1);g.moveTo(px,y2);g.lineTo(px,yB);g.stroke();}
  for(let d=Math.ceil(vlo);d<=Math.floor(vhi);d++){
    const py=Y(Math.pow(10,d)); if(py<y0||py>y1)continue;
    g.beginPath();g.moveTo(x0,py);g.lineTo(x1,py);g.stroke();}
  g.fillStyle=css('--mute'); g.textAlign='center'; g.textBaseline='top';
  for(const t of T){const px=xof(t); if(px<x0-1||px>x1+1)continue; g.fillText(t.lab,px,yB+4);}
  g.textAlign='right'; g.textBaseline='middle';
  for(let d=Math.ceil(vlo);d<=Math.floor(vhi);d++){
    const py=Y(Math.pow(10,d)); if(py<y0||py>y1)continue;
    g.fillText(d<=0?'1e'+d:(d<4?String(Math.pow(10,d)):'1e'+d),x0-5,py);}
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText('amplitude V',0,0); g.restore();

  // ---- upper: the sounding curve
  g.save(); g.beginPath(); g.rect(x0,y0,x1-x0,y1-y0); g.clip();
  if(F.haspred){
    g.strokeStyle=C_P50; g.lineWidth=1.8; g.beginPath();
    qs.forEach((q,i)=>{const px=X(q),py=Y(F.p50[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
    g.stroke();
  }
  g.strokeStyle=C_BEST; g.lineWidth=1.4; g.setLineDash([5,3]); g.beginPath();
  qs.forEach((q,i)=>{const px=X(q),py=Y(F.sig[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
  g.stroke(); g.setLineDash([]);
  g.fillStyle=C_TRUE;
  qs.forEach((q,i)=>{g.beginPath(); g.arc(X(q),Y(F.V[i]),2.6,0,6.284); g.fill();});
  g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y0,x1-x0,y1-y0);
  g.fillStyle=css('--mute'); g.font='9px ui-monospace,monospace';
  g.textAlign='left'; g.textBaseline='top'; g.fillText('noise',x0+4,Y(F.sig[0])+3);

  // ---- lower: residuals in units of the noise
  const res = F.haspred ? qs.map((q,i)=>(F.V[i]-F.p50[i])/F.sig[i]) : qs.map(()=>0);
  const rhi = F.haspred
    ? Math.max(3, Math.ceil(Math.max(...res.map(Math.abs),
        ...qs.map((q,i)=>Math.abs((F.V[i]-F.p05[i])/F.sig[i])),
        ...qs.map((q,i)=>Math.abs((F.V[i]-F.p95[i])/F.sig[i])))))
    : 3;
  const R=r=>y2+(yB-y2)/2-(r/rhi)*((yB-y2)/2-3);
  g.save(); g.beginPath(); g.rect(x0,y2,x1-x0,yB-y2); g.clip();
  g.fillStyle='rgba(255,255,255,.13)'; g.fillRect(x0,R(1),x1-x0,R(-1)-R(1));
  if(F.haspred){
    g.fillStyle='rgba(0,229,255,.22)'; g.beginPath();
    qs.forEach((q,i)=>{const px=X(q),py=R((F.V[i]-F.p95[i])/F.sig[i]); i?g.lineTo(px,py):g.moveTo(px,py);});
    for(let i=qs.length-1;i>=0;i--) g.lineTo(X(qs[i]),R((F.V[i]-F.p05[i])/F.sig[i]));
    g.closePath(); g.fill();
  }
  g.strokeStyle='rgba(255,255,255,.55)'; g.lineWidth=1;
  g.beginPath(); g.moveTo(x0,R(0)); g.lineTo(x1,R(0)); g.stroke();
  if(F.haspred){
    g.strokeStyle=C_TRUE; g.lineWidth=1.2; g.fillStyle=C_TRUE;
    qs.forEach((q,i)=>{const px=X(q);
      g.beginPath(); g.moveTo(px,R(0)); g.lineTo(px,R(res[i])); g.stroke();
      g.beginPath(); g.arc(px,R(res[i]),2.4,0,6.284); g.fill();});
  }
  g.restore();
  g.strokeStyle=css('--line'); g.lineWidth=1; g.strokeRect(x0,y2,x1-x0,yB-y2);
  g.fillStyle=css('--mute'); g.font='9.5px ui-monospace,monospace';
  g.textAlign='right'; g.textBaseline='middle';
  g.fillText('+'+rhi,x0-5,R(rhi*0.92)); g.fillText('0',x0-5,R(0)); g.fillText('-'+rhi,x0-5,R(-rhi*0.92));
  g.save(); g.translate(11,(y2+yB)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle';
  g.fillText('residual / '+String.fromCharCode(963),0,0); g.restore();
  g.textAlign='right'; g.textBaseline='top';
  g.fillText('pulse moment q',x1,yB+15);
  // top of the upper panel: the bottom-left corner already carries 'noise'
  g.fillStyle='rgba(255,255,255,.72)'; g.font='9.5px ui-monospace,monospace';
  g.textAlign='right'; g.textBaseline='top';
  g.fillText(F.haspred
    ? ('median '+String.fromCharCode(967)+String.fromCharCode(178)+'/N = '+F.chi2n.toFixed(2))
    : 'no kernel stored: measurements and noise only', x1-5, y0+4);
}
function drawGate(){
  const {g,w,h}=setup(document.getElementById('c3')); if(w<8) return;
  const xof=t=>PAD.l+((t.v-D.aem.lo)/(D.aem.hi-D.aem.lo))*(w-PAD.l-PAD.r);
  const A=frame(g,w,h,'conductivity S/m',LOGT,xof,'depth m');
  const X=sxLog(A);
  band(g,A);
  for(const k of admIdx){
    if(k===sel) continue;
    // a node that fell back is still drawn, in magenta, rather than passed off
    // as a smoothest model
    steps(g,A,X,modelOf(k).map(v=>Math.pow(10,v)),aemTop,aemBot,
          (useSmooth&&fellback(k))?C_FB_FAINT:C_FAM,0.8);
  }
  const bi=D.gate.C.indexOf(D.gate.best);
  if(bi>=0 && bi!==sel)
    // the best fit is the default selection now, so the amber line would sit
    // under the green one; draw it only when the reader has moved off it
    steps(g,A,X,modelOf(bi).map(v=>Math.pow(10,v)),aemTop,aemBot,C_BEST,1.4,[5,3]);
  smooth(g,A,X,D.aem.truth,aemTop,aemBot,C_TRUE,2,null,true);
  const m=modelOf(sel).map(v=>Math.pow(10,v));
  steps(g,A,X,m,aemTop,aemBot,C_SEL,3.2,null,true);
  for(const i of obOf(sel)){
    if(aemTop[i]>D.zmax) continue;
    const Y=(yOf(A,aemTop[i])+yOf(A,Math.min(aemBot[i],D.zmax)))/2;
    g.beginPath(); g.arc(X(m[i]),Y,4.4,0,6.284);
    g.fillStyle='rgba(0,0,0,.8)'; g.fill();
    g.lineWidth=1.8; g.strokeStyle=C_LIM; g.stroke();
  }
  // the sensitivity veil goes with the sensitivity line, which figure mode
  // does not draw; at a 100 m axis it would otherwise dim test one below
  // 84.6 m with no line to explain why
  if(!FIGMODE) veil(g,A);
  if(useSmooth && fellback(sel)){
    g.save(); g.fillStyle=C_FB; g.font='9.5px ui-monospace,monospace';
    g.textAlign='left'; g.textBaseline='top';
    g.fillText('no smoother model fits', A.x0+4, A.y0+4); g.restore();
  }
  g.strokeStyle=css('--line'); g.lineWidth=1;
  g.strokeRect(A.x0,A.y0,A.x1-A.x0,A.y1-A.y0);
}
function drawChi(){
  const {g,w,h}=setup(document.getElementById('c4')); if(w<8) return;
  const C=D.gate.C, lC=C.map(Math.log10);
  const lo=Math.min(...lC)-0.08, hi=Math.max(...lC)+0.08;
  const y0=8,y1=h-26,x0=40,x1=w-8;
  const chis=D.gate.C.map((c,i)=>chiOf(i)).filter(isFinite);
  // The three rejected nodes reach chi2 49, 234 and 932. Scaling to them
  // flattens everything that matters; CLAMPING them to the top -- which this
  // did -- draws three different misfits as one horizontal line at a value
  // none of them has. Instead the axis covers the decision region, the curve
  // is clipped so it genuinely leaves the top, and each node above the axis
  // is marked with its own value.
  const cmin=Math.min(...chis)-1.5, cmax=D.tabs+8;
  const X=c=>x0+((Math.log10(c)-lo)/(hi-lo))*(x1-x0);
  const Y=v=>y1-((v-cmin)/(cmax-cmin))*(y1-y0);
  g.font='9.5px ui-monospace,monospace'; g.fillStyle=css('--mute');
  g.strokeStyle=css('--grid'); g.lineWidth=1;
  const TK=[0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1,2,5].filter(v=>Math.log10(v)>=lo&&Math.log10(v)<=hi);
  for(const v of TK){const x=X(v);
    g.beginPath();g.moveTo(x,y0);g.lineTo(x,y1);g.stroke();
    g.textAlign='center';g.textBaseline='top';
    g.fillText(v<1?String(v).replace(/^0/,''):String(v),x,y1+4);}
  g.strokeStyle=css('--hot'); g.setLineDash([5,3]); g.beginPath();
  g.moveTo(x0,Y(D.tabs)); g.lineTo(x1,Y(D.tabs)); g.stroke();
  g.strokeStyle=css('--acc2'); g.beginPath();
  g.moveTo(x0,Y(D.trel)); g.lineTo(x1,Y(D.trel)); g.stroke(); g.setLineDash([]);
  g.textAlign='left'; g.textBaseline='bottom'; g.fillStyle=css('--hot');
  g.fillText('N + 3.84 = '+F(D.tabs,2), x0+3, Y(D.tabs)-1);
  g.fillStyle=css('--acc2');
  g.fillText('best fit + 3.84 = '+F(D.trel,2), x0+3, Y(D.trel)-1);
  g.strokeStyle=css('--grid'); g.textAlign='right'; g.textBaseline='middle';
  g.fillStyle=css('--mute');
  for(const v of [30,40,50]){ if(v<cmin||v>cmax) continue;
    g.beginPath(); g.moveTo(x0,Y(v)); g.lineTo(x1,Y(v)); g.stroke();
    g.fillText(String(v),x0-4,Y(v)); }
  g.strokeStyle=css('--line'); g.beginPath();
  g.moveTo(x0,y0); g.lineTo(x0,y1); g.lineTo(x1,y1); g.stroke();
  g.fillStyle=css('--mute'); g.textAlign='right'; g.textBaseline='top';
  g.fillText(D.field?'window conductivity S/m (geometric mean)':'window conductance S',x1,y1+14);
  g.save(); g.translate(11,(y0+y1)/2); g.rotate(-Math.PI/2);
  g.textAlign='center'; g.textBaseline='middle'; g.fillText('min '+String.fromCharCode(967)+String.fromCharCode(178),0,0); g.restore();
  g.save(); g.beginPath(); g.rect(x0,y0-2,x1-x0,y1-y0+2); g.clip();
  g.strokeStyle=css('--acc'); g.lineWidth=1.8; g.beginPath();
  C.forEach((c,i)=>{const px=X(c),py=Y(chiOf(i)); i?g.lineTo(px,py):g.moveTo(px,py);});
  g.stroke();
  C.forEach((c,i)=>{ if(chiOf(i)>cmax) return;
    const px=X(c), py=Y(chiOf(i));
    if(useSmooth && fellback(i)){            // marked, not silently substituted
      g.beginPath(); g.moveTo(px,py-4); g.lineTo(px+4,py); g.lineTo(px,py+4); g.lineTo(px-4,py);
      g.closePath(); g.fillStyle=C_FB; g.fill();
      g.lineWidth=1; g.strokeStyle=css('--panel'); g.stroke();
    } else {
      g.beginPath();g.arc(px,py,2.6,0,6.284);
      g.fillStyle=D.gate.adm[i]?css('--acc'):css('--mute');g.fill();
    }});
  g.restore();
  // nodes off the top, each with its own value
  g.fillStyle=css('--mute'); g.font='9px ui-monospace,monospace';
  g.textAlign='center'; g.textBaseline='top';
  C.forEach((c,i)=>{ if(chiOf(i)<=cmax) return;
    const px=X(c);
    g.beginPath(); g.moveTo(px,y0+1); g.lineTo(px-3.4,y0+6); g.lineTo(px+3.4,y0+6);
    g.closePath(); g.fill();
    g.fillText(String(Math.round(chiOf(i))),px,y0+7);});
  if(D.truthC!==null && D.truthC!==undefined && isFinite(D.truthC)){
  g.strokeStyle=css('--hot'); g.lineWidth=1.3; g.setLineDash([2,2]); g.beginPath();
  g.moveTo(X(D.truthC),y0); g.lineTo(X(D.truthC),y1); g.stroke(); g.setLineDash([]);}
  // the slider's position, as a marker that moves along the curve
  const sc=D.gate.C[sel], sy=Y(chiOf(sel)), sx=X(sc);
  g.strokeStyle=css('--ink'); g.globalAlpha=.45; g.lineWidth=1.2;
  g.beginPath(); g.moveTo(sx,y0); g.lineTo(sx,y1); g.stroke();
  g.beginPath(); g.moveTo(x0,sy); g.lineTo(sx,sy); g.stroke(); g.globalAlpha=1;
  g.beginPath(); g.arc(sx,sy,7.5,0,6.284);
  g.fillStyle=C_SEL; g.fill();
  g.lineWidth=2.4; g.strokeStyle=css('--panel'); g.stroke();
  g.lineWidth=1.2; g.strokeStyle=css('--ink'); g.beginPath(); g.arc(sx,sy,8.8,0,6.284); g.stroke();
  g.font='10px ui-monospace,monospace'; g.fillStyle=css('--ink');
  g.textAlign = sx>(x0+x1)/2 ? 'right' : 'left';
  g.textBaseline='bottom';
  g.fillText(String.fromCharCode(967)+String.fromCharCode(178)+' '+F(chiOf(sel),2),
             sx+(sx>(x0+x1)/2?-13:13), sy+4);
}

// ---- conductance over a depth range, for one model (S/m per cell)
function cond(m,a,b){let s=0;
  for(let i=0;i<m.length;i++){
    const ov=Math.max(0,Math.min(aemBot[i],b)-Math.max(aemTop[i],a));
    if(ov>0) s+=m[i]*ov;}
  return s;}

// FIELD: a site with no detected interval has no profile and no gate. Draw the
// two ensembles and the data fit, say why the rest is missing, and stop here.
if(!D.gate){
  document.querySelector('.strip').hidden=true;
  document.querySelectorAll('.pan')[2].hidden=true;
  document.querySelector('.panels').classList.add('two');
  document.getElementById('prov').innerHTML=D.prov||'';
  if(D.summ){const sm=document.getElementById('summ'); sm.innerHTML=D.summ; sm.hidden=false;}
  document.getElementById('n1').textContent=D.aem.ndraw.toLocaleString()+' draws';
  let tb='wc';
  const d2=()=>{ tb==='wc' ? drawSNMR() : drawFit(); };
  const st=w=>{ tb=w;
    document.getElementById('t2a').classList.toggle('on',w==='wc');
    document.getElementById('t2b').classList.toggle('on',w==='fit');
    document.getElementById('n2').textContent = w==='wc'
      ? D.snmr.ndraw.toLocaleString()+' draws' : (D.fit ? D.fit.q.length+' gates, '+D.fit.src : 'no data file');
    d2(); };
  document.getElementById('t2a').addEventListener('click',()=>st('wc'));
  document.getElementById('t2b').addEventListener('click',()=>st('fit'));
  if(!D.fit) document.getElementById('t2b').disabled=true;
  const rd=()=>{drawAEM(); d2();};
  st('wc'); rd();
  let tt=null; addEventListener('resize',()=>{clearTimeout(tt);tt=setTimeout(rd,120);});
  return;
}
const admIdx = D.gate.C.map((c,i)=>i).filter(i=>D.gate.adm[i]);
let useSmooth=false;
// A node whose smoothest search fell back carries the min-chi2 model and is
// MARKED, not quietly substituted: the fan would otherwise mix two kinds of
// model with nothing to say which is which.
const fellback = k => !!(D.gate.fb && D.gate.fb[k]);
const hasSm    = k => !!(D.gate.ms && D.gate.ms[k] && !fellback(k));
const modelOf  = k => (useSmooth && hasSm(k)) ? D.gate.ms[k] : D.gate.m[k];
const chiOf    = k => (useSmooth && hasSm(k) && D.gate.mschi2 && D.gate.mschi2[k]!==null)
  ? D.gate.mschi2[k] : D.gate.chi2[k];
// the stored `ob` lists belong to the min-chi2 models; for the smoothest ones
// recompute with the same predicate the gate uses
const obOf = k => (useSmooth && hasSm(k))
  ? D.gate.ms[k].map((v,i)=>[v,i]).filter(([v])=>v<=D.mb[0]+1e-9||v>=D.mb[1]-1e-9).map(([,i])=>i)
  : D.gate.ob[k];
// OPEN ON THE LOWEST-MISFIT ADMISSIBLE NODE. Chosen 19 Nov after the previous
// default -- the node nearest the geometric centre of the admissible range --
// put the confined column on 0.5060 S, a node with 15 of 65 cells resting on a
// stated limit, sitting between two clean ones. That rule picks on the range's
// arithmetic alone and never looks at the model it is about to show; on a band
// spanning a factor of 126 it lands wherever the ends happen to be. It looked
// right on the unconfined column only by luck.
//
// The lowest-misfit node is defined by the DATA, not by the range: 1.6000 S on
// the confined column and 0.2014 S on the unconfined one, and neither has a
// cell on a limit. It is a representative of the admissible set, not the
// answer -- the slider moves through every other node.
let sel = admIdx.reduce((a,b)=> D.gate.chi2[b] < D.gate.chi2[a] ? b : a, admIdx[0]);

// Conductance above and below the window, as a factor of the lowest-misfit
// model: the compensation the forced conductance is paid for with. Exported
// per node rather than shown, so the table stays still.
const ZB=Math.min(D.zmax,(D.zb!==undefined&&D.zb!==null)?D.zb:80.7);
function ratios(k){
  const bi=D.gate.C.indexOf(D.gate.best);
  const mk=modelOf(k).map(v=>Math.pow(10,v)), bm=modelOf(bi).map(v=>Math.pow(10,v));
  return [cond(mk,0,D.win[0])/cond(bm,0,D.win[0]),
          cond(mk,D.win[1],ZB)/cond(bm,D.win[1],ZB)];
}
function say(){
  const C=D.gate.C[sel], chi=chiOf(sel);
  const mm=modelOf(sel).map(v=>Math.pow(10,v));
  const th=+document.getElementById('sP').value, ex=+document.getElementById('sM').value;
  const kf=+document.getElementById('sK').value/100;
  const sigb=C/D.thick, clay=kf*sigb;
  // sigma_w = (sigma_bulk - clay) / theta^m, the relation Stage 5 uses. Clay
  // can only ADD conductance, so clay = 0 gives the MOST saline reading the
  // data allow; every step of this slider lowers it.
  const tds=(sigb-clay)/Math.pow(th,ex)*D.pet.mgl;
  const ceilK=D.pet.scmax*D.pet.c0*Math.max(0,1-th/D.pet.phiref);
  const overK=(D.pet.scmax===null)?false:clay>Math.min(ceilK,D.pet.cp[1]);
  const U=D.field?' S/m':' S', TV=v=>(v===null||v===undefined)?'&middot;':v;
  const ob=obOf(sel);
  const row=(k,v,t)=>'<tr><td class="k">'+k+'</td><td class="v">'+v+
    '</td><td class="t">'+(t===undefined?'&middot;':t)+'</td></tr>';
  document.getElementById('say').innerHTML =
    '<colgroup><col class="k"><col class="v"><col class="t"></colgroup>'+
    '<thead><tr><th class="k"></th><th>selected</th><th>'+(D.field?'':'true')+'</th></tr></thead><tbody>'+
    row(D.field?'window conductivity':'window conductance', F(C,4)+U, D.field?undefined:F(D.truthC,4)+' S')+
    row('misfit', String.fromCharCode(967)+String.fromCharCode(178)+' '+F(chi,2)+
        (useSmooth&&fellback(sel)?' <em>no smoother model fits</em>':''))+
    row('on a limit', ob.length===0 ? '0 <em>of '+mm.length+'</em>'
        : ob.length+' <em>of '+mm.length+', from '+F(aemTop[ob[0]],1)+' m</em>')+
    row('clay term', F(clay,4)+' S/m <em>'+Math.round(kf*100)+'%'+
        (overK?', past the ceiling':'')+'</em>',
        D.claytrue?F(D.claytrue.term,4)+' S/m  '+Math.round(D.claytrue.frac*100)+'%':undefined)+
    row('pore water', Math.round(tds).toLocaleString()+' mg/L'+
        (kf===0?' <em>most saline allowed</em>':''),
        D.claytrue?Math.round(D.claytrue.tds).toLocaleString()+' mg/L':undefined)+
    '</tbody>';
  document.getElementById('oC').textContent=F(C,4)+U;
  document.getElementById('oP').textContent=F(th,3);
  document.getElementById('oM').textContent=F(ex,2);
  document.getElementById('oK').textContent=Math.round(kf*100)+'% / '+F(clay,4);
}

function csv(){
  const L=['block,key1,key2,value'];
  const FT=x=>(x===undefined||x===null)?'':F(x,6);
  const bw1=(D.aem.hi-D.aem.lo)/D.aem.nb;
  for(let i=0;i<D.aem.h.length;i++)for(let j=0;j<D.aem.nb;j++){const c=D.aem.h[i][j];
    if(c) L.push('aem_hist,'+F(D.aem.z[i],3)+','+F(D.aem.lo+(j+0.5)*bw1,4)+','+c);}
  for(let i=0;i<D.aem.z.length;i++) L.push('aem_lines,'+F(D.aem.z[i],3)+',p05_p50_p95_truth,'+
    F(D.aem.p05[i],6)+' '+F(D.aem.p50[i],6)+' '+F(D.aem.p95[i],6)+' '+FT(D.aem.truth[i]));
  const bw2=(D.snmr.hi-D.snmr.lo)/D.snmr.nb;
  for(let i=0;i<D.snmr.h.length;i++)for(let j=0;j<D.snmr.nb;j++){const c=D.snmr.h[i][j];
    if(c) L.push('snmr_hist,'+F(D.snmr.z[i],3)+','+F(D.snmr.lo+(j+0.5)*bw2,4)+','+c);}
  for(let i=0;i<D.snmr.z.length;i++) L.push('snmr_lines,'+F(D.snmr.z[i],3)+',p05_p50_p95_truth,'+
    F(D.snmr.p05[i],6)+' '+F(D.snmr.p50[i],6)+' '+F(D.snmr.p95[i],6)+' '+FT(D.snmr.truth[i]));
  for(let k=0;k<D.gate.C.length;k++){
    L.push('gate_node,'+F(D.gate.C[k],6)+',chi2_admissible,'+F(D.gate.chi2[k],6)+' '+D.gate.adm[k]);
    if(D.gate.adm[k]){ const r=ratios(k);
      L.push('gate_ratio,'+F(D.gate.C[k],6)+',above_0_'+F(D.win[0],1)+'_below_'+
             F(D.win[1],1)+'_'+F(ZB,1)+',' +F(r[0],6)+' '+F(r[1],6)); }
    const ob=new Set(D.gate.ob[k]);
    for(let i=0;i<D.gate.m[k].length;i++)
      L.push('gate_model,'+F(D.gate.C[k],6)+','+F(D.aem.z[i],3)+','+
        F(Math.pow(10,D.gate.m[k][i]),6)+' '+(ob.has(i)?1:0));}
  return L.join('\n');
}

let tab2='wc';
function draw2(){ tab2==='wc' ? drawSNMR() : drawFit(); }
function redraw(){drawAEM();draw2();drawGate();drawChi();say();}

const sC=document.getElementById('sC');
sC.max=String(admIdx.length-1);
sC.value=String(admIdx.indexOf(sel));
sC.addEventListener('input',()=>{sel=admIdx[+sC.value];drawGate();drawChi();say();});
(function(){
  const P=document.getElementById('sP'), M=document.getElementById('sM');
  P.min=D.pet.phi[0]; P.max=D.pet.phi[1]; P.value=(D.pet.phi[0]+D.pet.phi[1])/2;
  M.min=D.pet.m[0];  M.max=D.pet.m[1];  M.value=(D.pet.m[0]+D.pet.m[1])/2;
})();
for(const id of ['sP','sM','sK']) document.getElementById(id).addEventListener('input',say);
(function(){
  const b=document.getElementById('bM'), h=document.getElementById('hM');
  const have=!!D.gate.ms;
  const nfb=have?D.gate.fb.reduce((a,b)=>a+b,0):0;
  b.textContent='lowest misfit';
  b.disabled=!have;
  if(have){
    b.addEventListener('click',()=>{
      useSmooth=!useSmooth; b.setAttribute('aria-pressed',String(useSmooth));
      b.textContent=useSmooth?'smoothest passing':'lowest misfit';
      h.textContent=useSmooth&&nfb ? nfb+' fell back' : '';
      drawGate(); drawChi(); say();});
  } else {
    b.title='the smoothest admissible model per node is not stored yet';
    h.textContent='not stored yet';
  }
})();

document.getElementById('bC').addEventListener('click',()=>{
  const text=csv(), ta=document.getElementById('ta');
  ta.value=text; ta.hidden=false; document.getElementById('bS').hidden=false;
  try{
    const u=URL.createObjectURL(new Blob([text],{type:'text/csv'}));
    const a=document.createElement('a'); a.href=u; a.download=(D.csvname||'report_3P_data.csv');
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(()=>URL.revokeObjectURL(u),2000);
  }catch(e){}
});
document.getElementById('bS').addEventListener('click',()=>{
  const ta=document.getElementById('ta'); ta.focus(); ta.select();});

if(D.field){ document.getElementById('prov').innerHTML = D.prov;
  const sm=document.getElementById('summ'); sm.innerHTML=D.summ; sm.hidden=false;
  document.getElementById('lC').textContent='window conductivity';
} else document.getElementById('prov').innerHTML =
  'AEM '+D.aem.src+' &times;'+D.aem.nens+', '+D.aem.ndraw.toLocaleString()+' draws &middot; '+
  'SNMR '+D.snmr.src+', '+D.snmr.ndraw.toLocaleString()+' draws <span class="flag">'+D.snmr.note+'</span> &middot; '+
  'gate '+D.gate.src+', '+Math.pow(10,D.mb[0]).toPrecision(2)+'&ndash;'+
  Math.pow(10,D.mb[1]).toPrecision(2)+' S/m &middot; '+D.gen;
document.getElementById('n1').textContent=D.aem.ndraw.toLocaleString()+' draws';
function setTab(which){
  tab2=which;
  const a=document.getElementById('t2a'), b=document.getElementById('t2b');
  a.classList.toggle('on',which==='wc'); b.classList.toggle('on',which==='fit');
  a.setAttribute('aria-pressed',String(which==='wc'));
  b.setAttribute('aria-pressed',String(which==='fit'));
  document.getElementById('n2').textContent = which==='wc'
    ? D.snmr.ndraw.toLocaleString()+' draws'
    : (D.fit ? D.fit.q.length+' gates, '+D.fit.src : 'no data file');
  draw2();
}
document.getElementById('t2a').addEventListener('click',()=>setTab('wc'));
document.getElementById('t2b').addEventListener('click',()=>setTab('fit'));
if(!D.fit) document.getElementById('t2b').disabled=true;
document.getElementById('n2').textContent=D.snmr.ndraw.toLocaleString()+' draws';
document.getElementById('n3').textContent=admIdx.length+' of '+D.gate.C.length+' nodes';
document.getElementById('colname').textContent = D.colname || '';
['a','b','c'].forEach((t,n)=>{const p=document.querySelectorAll('.pan')[n];
  if(p) p.setAttribute('data-letter',t);});
document.getElementById('key').innerHTML=
  '<span><i style="background:'+C_P50+';outline:1px solid '+css('--rule')+'"></i>p50</span>'+
  '<span><i style="background:'+C_PTL+';outline:1px solid '+css('--rule')+'"></i>p05 p95</span>'+
  (D.field?'':'<span><i style="background:'+C_TRUE+';outline:1px solid '+css('--rule')+'"></i>true</span>')+'<span><i style="background:'+C_SEL+'"></i>gate model shown</span>'+(D.gate.C.indexOf(D.gate.best)!==sel ? '<span><i style="background:'+C_BEST+'"></i>best fit</span>' : '')+
  '<span><i style="background:transparent;border:1.6px solid '+C_LIM+';outline:1px solid '+css('--rule')+';height:8px;width:8px;border-radius:50%"></i>on a stated limit</span>'+
  '<span>window '+D.win[0]+'&ndash;'+D.win[1]+' m</span>';

redraw();
let t=null;
addEventListener('resize',()=>{clearTimeout(t);t=setTimeout(redraw,120);});
try{
  const mq=matchMedia('(prefers-color-scheme:dark)');
  const onTheme=()=>setTimeout(redraw,30);
  mq.addEventListener ? mq.addEventListener('change',onTheme) : mq.addListener(onTheme);
}catch(e){}
})();
</script>
