"""Frozen pilot geometry methods; no import-time map processing.

See pilot_methods_provenance.json. New campaign parameters are in protocol.json.
"""
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from scipy.signal import find_peaks
VOXEL=.05
REF=2.74
STATIONS=[3.,7.,20.]

def xyz(m):
    f={v.name:v for v in m.fields};types={7:'f4',8:'f8'};end='>' if m.is_bigendian else '<'
    dt=np.dtype({'names':['x','y','z'],'formats':[end+types[f[k].datatype] for k in ['x','y','z']],'offsets':[f[k].offset for k in ['x','y','z']],'itemsize':m.point_step})
    a=np.ndarray((m.height,m.width),dtype=dt,buffer=m.data,strides=(m.row_step,m.point_step))
    p=np.column_stack([a[k].reshape(-1) for k in ['x','y','z']])
    return p[np.isfinite(p).all(axis=1)].astype(np.float64)


def pcd(path):
    with path.open('rb') as f:
        hd={}
        while True:
            s=f.readline().decode().strip()
            if not s or s.startswith('#'):continue
            k,*v=s.split();hd[k]=v
            if k=='DATA':break
        names=hd['FIELDS']
        if hd['DATA']==['ascii']:
            a=np.loadtxt(f);p=a[:,[names.index(k) for k in ['x','y','z']]]
        else:
            assert hd['DATA']==['binary'],hd['DATA']
            counts=hd.get('COUNT',['1']*len(names))
            dt=np.dtype([(n,{'F':'f','I':'i','U':'u'}[t]+s,(int(c),)) if int(c)!=1 else (n,{'F':'f','I':'i','U':'u'}[t]+s) for n,t,s,c in zip(names,hd['TYPE'],hd['SIZE'],counts)])
            a=np.frombuffer(f.read(),dt);p=np.column_stack([a[k] for k in ['x','y','z']])
    return p[np.isfinite(p).all(axis=1)].astype(float)


def voxel(points):
    if not len(points):return points
    cells=np.floor(points/VOXEL).astype(np.int64)
    _,idx,inv=np.unique(cells,axis=0,return_index=True,return_inverse=True)
    # Preserve one original point per cell: no smoothing of the geometry.
    return points[idx]


def tls(p):
    center=np.mean(p,axis=0);_,v=np.linalg.eigh((p-center).T@(p-center)/len(p));n=v[:,0]
    if n[2]<0:n=-n
    d=(p-center)@n
    return center,n,dict(n=len(p),rmse_m=float(np.sqrt(np.mean(d*d))),mean_abs_m=float(np.mean(abs(d))),p95_abs_m=float(np.quantile(abs(d),.95)))


def normals(p):
    ix=np.linspace(0,len(p)-1,min(len(p),35000),dtype=int);s=p[ix].astype(float);tree=cKDTree(p)
    n=[];quality=[];radius=[]
    for b in range(0,len(s),3000):
        distances,ids=tree.query(s[b:b+3000],k=20)
        q=p[ids].astype(float);q-=q.mean(axis=1,keepdims=True);cov=np.einsum('nki,nkj->nij',q,q)/20
        vals,vecs=np.linalg.eigh(cov);n.append(vecs[:,:,0]);quality.append((vals[:,0]<.15*vals[:,1])&(vals[:,1]>.08*vals[:,2]));radius.append(distances[:,-1])
    return s,np.vstack(n),np.concatenate(quality)&(np.concatenate(radius)<.8)


def floor_seed(s,n,good,route):
    # First isolate the clear corridor interior in XY. Horizontal scan rings on
    # vertical walls must never be treated as a horizontal ground surface.
    xy=route[:,:2];origin=xy[0];_,v=np.linalg.eigh(np.cov((xy-origin).T));u=v[:,-1]
    far=xy[np.argmax(np.linalg.norm(xy-origin,axis=1))]-origin
    if np.dot(far,u)<0:u=-u
    R=np.array([[u[0],u[1],0],[-u[1],u[0],0],[0,0,1]])
    q=s@R.T;nn=n@R.T;x0=float((route[0]@R.T)[0])
    along=(q[:,0]>x0-5)&(q[:,0]<x0+12)
    wall=q[along&(abs(nn[:,1])>.65)]
    edges=np.arange(np.quantile(wall[:,1],.01)-.1,np.quantile(wall[:,1],.99)+.11,.05)
    hist,_=np.histogram(wall[:,1],edges);peaks,_=find_peaks(hist,distance=12,prominence=max(2,hist.max()*.025));opts=[]
    for i in peaks:
        for j in peaks:
            if 1.5<edges[j]-edges[i]<4.5:opts.append((min(hist[i],hist[j])+.3*(hist[i]+hist[j]),i,j))
    assert opts,'no corridor side walls for floor selection'
    _,i,j=max(opts);center=(edges[i]+edges[i+1]+edges[j]+edges[j+1])/4
    cand=s[along&(abs(q[:,1]-center)<.65)]
    assert len(cand)>50,('not enough floor interior candidates',len(cand))
    # Lower horizontal structures only. Random plane seeds are deterministic.
    cand=cand[cand[:,2]<np.quantile(cand[:,2],.45)]
    rng=np.random.RandomState(7301);best=None
    centerxy=np.median(cand[:,:2],axis=0)
    for _ in range(500):
        p=cand[rng.choice(len(cand),3,replace=False)];v=np.cross(p[1]-p[0],p[2]-p[0]);norm=np.linalg.norm(v)
        if norm<1e-6:continue
        v/=norm
        if abs(v[2])<.94:continue
        residual=abs((cand-p[0])@v);mask=residual<.08
        if mask.sum()<50:continue
        q=cand[mask];spread=np.ptp(q[:,:2],axis=0).max()
        if spread<8:continue
        score=int(mask.sum())
        if best is None or score>best[0]:best=(score,q)
    assert best is not None,'floor seed not found'
    c,v,st=tls(best[1]);axis=np.cross(v,[0.,0.,1.]);angle=np.arccos(np.clip(v[2],-1,1));rot=Rotation.from_rotvec(axis/max(np.linalg.norm(axis),1e-15)*angle).as_matrix()
    return rot,float((rot@c)[2]),dict(seed_count=best[0],seed_plane_normal=v.tolist(),seed_plane_center=c.tolist(),seed_rmse=st['rmse_m'],clear_interior_half_width_m=.65)


def wall_pair(s,n,good,route):
    # Trajectory provides an approximate longitudinal direction, never a metric reference.
    xy=route[:,:2];origin=xy[0];_,v=np.linalg.eigh(np.cov((xy-origin).T));u=v[:,-1]
    far=xy[np.argmax(np.linalg.norm(xy-origin,axis=1))]-origin
    if np.dot(far,u)<0:u=-u
    R=np.array([[u[0],u[1],0],[-u[1],u[0],0],[0,0,1]])
    q=s@R.T;nn=n@R.T
    # Locate two long parallel side walls from vertical surface candidates.
    x0=float((route[0]@R.T)[0]);xlo=x0-5.;xhi=x0+12.
    for _ in range(3):
        q=s@R.T;nn=n@R.T
        mask=(abs(nn[:,1])>.65)&(q[:,0]>xlo)&(q[:,0]<xhi)&(q[:,2]>.35)&(q[:,2]<2.15)
        w=q[mask];assert len(w)>80,('no wall candidates',len(w))
        edges=np.arange(np.quantile(w[:,1],.01)-.1,np.quantile(w[:,1],.99)+.11,.05)
        hist,_=np.histogram(w[:,1],edges);peaks,_=find_peaks(hist,distance=12,prominence=max(2,hist.max()*.025))
        opts=[]
        for i in peaks:
            for j in peaks:
                gap=edges[j]-edges[i]
                if 1.5<gap<4.5:opts.append((min(hist[i],hist[j])+.3*(hist[i]+hist[j]),i,j))
        assert opts,('wall pair not found',len(peaks))
        _,i,j=max(opts);lower=float((edges[i]+edges[i+1])/2);upper=float((edges[j]+edges[j+1])/2)
        slopes=[]
        for y in [lower,upper]:
            sel=w[abs(w[:,1]-y)<.28]
            A=np.column_stack([sel[:,0]-sel[:,0].mean(),np.ones(len(sel))]);a,_=np.linalg.lstsq(A,sel[:,1],rcond=None)[0];slopes.append(a)
        angle=np.arctan(np.mean(slopes));c=np.cos(angle);ss=np.sin(angle);adjust=np.array([[c,ss,0],[-ss,c,0],[0,0,1]])
        R=adjust@R
    q=s@R.T;nn=n@R.T
    mask=(abs(nn[:,1])>.65)&(q[:,0]>xlo)&(q[:,0]<xhi)&(q[:,2]>.35)&(q[:,2]<2.15)
    w=q[mask];centers=[]
    for old in [lower,upper]:
        sel=w[abs(w[:,1]-old)<.38];centers.append(float(np.median(sel[:,1])))
    ty=REF/2-np.mean(centers)
    q[:,1]+=ty
    # Longitudinal anchor: the end wall spanning the clear corridor interior.
    endmask=(q[:,1]>.35)&(q[:,1]<REF-.35)&(q[:,0]<x0-2)&(q[:,0]>x0-18)&(q[:,2]>.3)&(q[:,2]<2.1)
    e=q[endmask];assert len(e)>8,('end wall anchor absent',len(e))
    edges=np.arange(e[:,0].min()-.1,e[:,0].max()+.11,.08);hist,_=np.histogram(e[:,0],edges);k=int(np.argmax(hist));ex=(edges[k]+edges[k+1])/2
    e=e[abs(e[:,0]-ex)<.3];ex=float(np.median(e[:,0]));tx=-ex
    return R,np.array([tx,ty,0.]),dict(side_wall_centers_before_translation=centers,side_wall_seed_gap_m=float(centers[1]-centers[0]),end_wall_anchor_x_before_translation=ex,end_wall_anchor_point_count=len(e),initial_long_axis_from='own trajectory PCA; direction only; no position accuracy reference')


def assess(p, corridor_margin=.6, branch_margin=.70, probe_z=1.2):
    result=[];floor_all=[];wall_all=[]
    for x in STATIONS:
        station={'center_x_m':x,'longitudinal_half_width_m':.5,'walls':[]}
        for y in [0.,REF]:
            mask=(abs(p[:,0]-x)<=.5)&(abs(p[:,1]-y)<=corridor_margin)&(p[:,2]>=.35)&(p[:,2]<=2.15)
            q=p[mask];bins=np.unique(np.floor((q[:,0]-(x-.5))/.2)).size if len(q) else 0
            if len(q)<12 or bins<4:
                station['walls'].append({'valid':False,'point_count':len(q),'longitudinal_20cm_bins':bins});continue
            c,n,st=tls(q)
            valid=abs(n[1])>.7
            est=float(c[1]-(n[0]*(x-c[0])+n[2]*(probe_z-c[2]))/n[1]) if valid else None
            station['walls'].append(dict(valid=bool(valid),point_count=len(q),longitudinal_20cm_bins=bins,plane_center=c.tolist(),plane_normal=n.tolist(),y_at_station_and_1p2m_height=est,**st))
            if valid:wall_all.append(q)
        valid=all(w['valid'] for w in station['walls'])
        if valid:
            width=station['walls'][1]['y_at_station_and_1p2m_height']-station['walls'][0]['y_at_station_and_1p2m_height']
            station.update(width_m=width,signed_width_error_m=width-REF,abs_width_error_m=abs(width-REF),percent_error=100*abs(width-REF)/REF)
        else:station.update(width_m=None,signed_width_error_m=None,abs_width_error_m=None,percent_error=None)
        f=p[(abs(p[:,0]-x)<=1.5)&(p[:,1]>=.65)&(p[:,1]<=REF-.65)&(abs(p[:,2])<=.6)]
        if len(f)>=20:
            c,n,st=tls(f);station['floor']=dict(valid=bool(abs(n[2])>.8),plane_center=c.tolist(),plane_normal=n.tolist(),**st)
            if abs(n[2])>.8:floor_all.append(f)
        else:station['floor']=dict(valid=False,n=len(f))
        result.append(station)
    err=np.array([s['signed_width_error_m'] for s in result if s['width_m'] is not None]);widths=[s['width_m'] for s in result if s['width_m'] is not None]
    # Broad continuous region includes every point within a fixed geometric window,
    # without rejecting points based on their fit residual.
    f=p[(p[:,0]>=3)&(p[:,0]<=20)&(p[:,1]>=.65)&(p[:,1]<=REF-.65)&(abs(p[:,2])<=.6)]
    floor=None
    if len(f)>=80:
        c,n,st=tls(f);bins=np.unique(np.floor((f[:,0]-3))).size
        floor=dict(valid=bool(abs(n[2])>.8 and bins>=12),plane_center=c.tolist(),plane_normal=n.tolist(),longitudinal_1m_bins=bins,**st)
    patches=[w for s in result for w in s['walls'] if w['valid']]
    wall_rms=float(np.sqrt(np.mean([w['rmse_m']**2 for w in patches]))) if patches else None
    # Lower alcove: three additional dimensions already present in the sketch.
    branch={};branch_def=[('left',0,21.99,1,(-2.65,-.55)),('right',0,27.69,1,(-2.65,-.55)),('back',1,-3.10,0,(22.59,27.09))]
    for key,axis,center,span_axis,span in branch_def:
        q=p[(abs(p[:,axis]-center)<branch_margin)&(p[:,span_axis]>span[0])&(p[:,span_axis]<span[1])&(p[:,2]>.35)&(p[:,2]<2.15)]
        bins=np.unique(np.floor((q[:,span_axis]-span[0])/.5)).size if len(q) else 0
        if len(q)<20 or bins<3:branch[key]=dict(valid=False,n=len(q),support_bins=bins);continue
        c,n,st=tls(q);valid=abs(n[axis])>.7
        at=np.array([24.84,-1.55,probe_z]);at[axis]=center
        estimate=float(c[axis]-(np.dot(n,at-c)-n[axis]*(at[axis]-c[axis]))/n[axis]) if valid else None
        branch[key]=dict(valid=bool(valid),plane_center=c.tolist(),plane_normal=n.tolist(),coordinate_at_probe=estimate,support_bins=bins,**st)
    dimensions=[]
    def dimension(key,label,reference,value):
        dimensions.append(dict(key=key,label=label,reference_m=reference,measured_m=value,signed_error_m=value-reference if value is not None else None,abs_error_m=abs(value-reference) if value is not None else None,percent_error=100*abs(value-reference)/reference if value is not None else None))
    for st in result:dimension('width_'+str(int(st['center_x_m'])),'Rộng hành lang tại x='+str(int(st['center_x_m']))+' m',REF,st['width_m'])
    left=branch['left'].get('coordinate_at_probe');right=branch['right'].get('coordinate_at_probe');back=branch['back'].get('coordinate_at_probe')
    dimension('straight_length','Dài cạnh thẳng tới nhánh dưới',21.99,left)
    dimension('branch_width','Rộng nhánh dưới',5.70,right-left if right is not None and left is not None else None)
    dimension('branch_depth','Tọa độ chiếu mặt lưng nhánh sau căn (chẩn đoán)',3.10,-back if back is not None else None)
    errs=np.array([d['signed_error_m'] for d in dimensions if d['signed_error_m'] is not None])
    return dict(stations=result,valid_width_stations=len(err),width_mean_m=float(np.mean(widths)) if len(err) else None,width_mae_m=float(np.mean(abs(err))) if len(err) else None,width_rmse_m=float(np.sqrt(np.mean(err*err))) if len(err) else None,width_max_abs_m=float(max(abs(err))) if len(err) else None,width_mean_percent=float(np.mean(abs(err))/REF*100) if len(err) else None,wall_patch_rmse_m=wall_rms,floor=floor,branch_planes=branch,dimensions=dimensions,dimension_mae_m=float(np.mean(abs(errs))) if len(errs) else None,dimension_rmse_m=float(np.sqrt(np.mean(errs**2))) if len(errs) else None,dimension_max_abs_m=float(np.max(abs(errs))) if len(errs) else None,valid_dimensions=len(errs)),f

