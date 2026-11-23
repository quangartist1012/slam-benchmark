#!/usr/bin/env python3
"""Five compact multi-panel research figures from retained, verified observations."""
from pathlib import Path
import csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle,FancyBboxPatch
from mpl_toolkits.mplot3d import Axes3D

O=Path(__file__).resolve().parent;R=O.parent;C=O/'data';F=O/'generated';F.mkdir(exist_ok=True)
REP='manual_20261113_164656';SEG=[((0,0),(21.99,0)),((0,2.74),(21.96,2.74)),((0,0),(0,2.74)),((21.99,0),(21.99,-3.1)),((21.99,-3.1),(27.69,-3.1)),((27.69,-3.1),(27.69,0))]
plt.rcParams.update({'font.size':8.5,'axes.titlesize':9,'axes.labelsize':8.5,'xtick.labelsize':7.5,'ytick.labelsize':8,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def read(p):return json.loads(p.read_text())
def csvread(p):
    with p.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def save(fig,stem):
    for ext in ['png','pdf']:fig.savefig(str(F/(stem+'.'+ext)),dpi=300,facecolor='white')
    plt.close(fig)
def reference(ax):
    for a,b in SEG:ax.plot([a[0],b[0]],[a[1],b[1]],color='#2b759b',lw=.8,zorder=5)
    ax.set_xlim(-1.2,32.5);ax.set_ylim(-5.1,4.4);ax.set_aspect('equal')
def config_axis(ax,names):
    ax.set_yticks(range(len(names)));ax.set_yticklabels(names);ax.set_ylim(17.6,-.6);ax.grid(axis='x',alpha=.18)
def main():
    s=read(C/'geometry_summary.json');names=s['configurations'];bags=s['bag_ids'];metric={(r['name'],r['bag_id']):r for r in s['bag_metrics']}
    # Platform + one protocol + reference geometry + one actual 3D cloud.
    fig=plt.figure(figsize=(7,4.35));grid=fig.add_gridspec(2,2,width_ratios=[1,1],height_ratios=[1,.97],left=.04,right=.99,top=.95,bottom=.08,hspace=.32,wspace=.22)
    ax=fig.add_subplot(grid[0,0]);im=plt.imread(str(C/'robot_platform.png'));ax.imshow(im[:,625:]);ax.axis('off');ax.set_title('(a) FoodMate sensing platform',loc='left',fontsize=8.5)
    ax=fig.add_subplot(grid[0,1]);ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
    for y,text,color in [(.72,'10 physical recordings → 18 configurations','#eaf1f6'),(.40,'180 map/pose/log outputs at 1× replay','#e8f3ed'),(.08,'3 groups: geometry · computation · continuity','#fff1df')]:
        ax.add_patch(FancyBboxPatch((.015,y),.97,.23,boxstyle='round,pad=.012',facecolor=color,edgecolor='#466578',lw=.7));ax.text(.5,y+.115,text,ha='center',va='center',fontsize=7.7)
    for y in [.70,.38]:ax.annotate('',xy=(.5,y-.05),xytext=(.5,y+.01),arrowprops=dict(arrowstyle='->',lw=.8))
    ax.set_title('(b) Replay and evaluation',loc='left',fontsize=8.5)
    ax=fig.add_subplot(grid[1,0]);reference(ax)
    for x,label in [(3,'S1'),(7,'S2'),(20,'S3')]:
        for y in [0,2.74]:ax.add_patch(Rectangle((x-.5,y-.6),1,1.2,facecolor='#b14f30',alpha=.25,edgecolor='none'))
        ax.text(x,3.55,label,ha='center',fontsize=7)
    for a,b,text,loc in [((0,-.8),(21.99,-.8),'21.99 m',(11,-1.9)),((11,0),(11,2.74),'2.74 m',(14,1.35)),((21.99,-3.8),(27.69,-3.8),'5.70 m',(24.8,-4.65))]:
        ax.annotate('',xy=a,xytext=b,arrowprops=dict(arrowstyle='<->',color='#2b759b',lw=.65));ax.text(*loc,text,ha='center',va='center',color='#245775',fontsize=7.7)
    ax.axis('off');ax.set_title('(c) Reference dimensions and wall stations',loc='left',fontsize=8.5)
    ax=fig.add_subplot(grid[1,1],projection='3d');p=np.load(str(C/'map_previews.npz'))['LIO-SAM_3d'];p=p[(p[:,0]>=-1)&(p[:,0]<=33)&(p[:,1]>=-5)&(p[:,1]<=4)&(p[:,2]>=-.2)&(p[:,2]<=3)]
    p=p[np.linspace(0,len(p)-1,min(14000,len(p)),dtype=int)];ax.scatter(p[:,0],p[:,1],p[:,2],c=p[:,2],s=.16,cmap='viridis',vmin=0,vmax=3,rasterized=True);ax.view_init(elev=32,azim=-72);ax.set_axis_off();ax.set_title('(d) LIO-SAM 3D cloud: 164656',loc='left',fontsize=8.5,pad=0)
    fig.text(.5,.015,'Resource comparison: one common recording, 1 s samples. Timing: 18 × 3 further replays of recording 164656.\nDrawing dimensions have no retained survey-uncertainty budget; cloud colour indicates height, not error.',ha='center',va='bottom',fontsize=7)
    save(fig,'01_setup_method')
    # Compact all-configuration map plate. These are actual aligned points, no scaling.
    fig,axs=plt.subplots(6,3,figsize=(7,6.4));fig.subplots_adjust(left=.015,right=.995,top=.965,bottom=.055,hspace=.25,wspace=.10)
    for ax,name in zip(axs.flat,names):
        p=np.load(str(C/'map_previews.npz'))[name];p=p[(p[:,2]>=.35)&(p[:,2]<=2.15)];p=p[np.linspace(0,len(p)-1,min(14000,len(p)),dtype=int)]
        ax.scatter(p[:,0],p[:,1],s=.08,color='#1b242c',alpha=.6,rasterized=True,linewidths=0);reference(ax);ax.axis('off');ax.set_title(name,fontsize=8.4,pad=2)
    fig.text(.5,.017,'Same recording 164656 · same extent and scale · black: wall-height points; blue: project reference.\nAll 180 maps, full-size 2D/3D views and measurement support are available in the public supplement.',ha='center',fontsize=7)
    save(fig,'02_maps_all_18')
    # Three geometry measures, each with bag-level points and n displayed for lengths.
    fig,axs=plt.subplots(1,3,figsize=(7,5.45),sharey=True)
    for i,n in enumerate(names):
        for ax,key in zip(axs,['core_dimension_mae_m','core_surface_rmse_m','mpv_m']):
            vv=[metric[n,b][key]*100 for b in bags if metric[n,b][key] is not None];ax.scatter(vv,i+np.linspace(-.15,.15,len(vv)),s=10,color='#297ca2',alpha=.78,linewidths=0)
            if vv:ax.plot(np.mean(vv),i,'D',ms=3.2,color='#263d4b')
        count=sum(metric[n,b]['core_dimension_mae_m'] is not None for b in bags);axs[0].text(.99,i,str(count),transform=axs[0].get_yaxis_transform(),va='center',ha='right',fontsize=6.7)
    for ax,title in zip(axs,['(a) Two-length MAE (cm)','(b) Four-patch RMSE (cm)','(c) MPV, r = 0.6 m (cm)']):ax.set_title(title,fontsize=8.3);ax.grid(axis='x',alpha=.18)
    config_axis(axs[0],names);axs[0].set_xlim(0,65);axs[1].set_xlim(left=0);axs[2].set_xlim(left=0)
    fig.subplots_adjust(left=.23,right=.99,top=.94,bottom=.10,wspace=.16)
    fig.text(.58,.015,'Dot: one recording; diamond: mean of available recordings. Numbers in (a): measured maps / 10.\nThe three measures capture different properties; low surface dispersion does not establish absolute accuracy.',ha='center',fontsize=7)
    save(fig,'03_geometry_combined')
    # CPU, RSS and timing are adjacent but explicitly from their respective campaigns.
    cpu={x['name']:x for x in read(C/'resources_common_recording.json')};lat={x['name']:x for x in read(C/'latency_summary.json')['results']}
    fig,axs=plt.subplots(1,3,figsize=(7,5.45),sharey=True)
    for i,n in enumerate(names):
        for ax,mean,peak in [(axs[0],'cpu_mean_percent_one_core','cpu_peak_percent_one_core'),(axs[1],'rss_mean_mib','rss_peak_mib')]:
            ax.plot([0,cpu[n][mean]],[i,i],lw=3,color='#5a98ba');ax.plot(cpu[n][peak],i,'|',color='#b94a31',ms=6)
        means=lat[n]['mean_ms']['values'];p95=lat[n]['p95_ms']['values'];yy=i+np.array([-.14,0,.14]);axs[2].scatter(means,yy,s=10,color='#297ca2',zorder=3);axs[2].scatter(p95,yy,s=14,facecolors='none',edgecolors='#b94a31',linewidths=.7,zorder=3)
    for ax,t in zip(axs,['(a) CPU: % of one core','(b) Summed RSS (MiB)','(c) Latency (ms, log scale)']):ax.set_title(t,fontsize=8.2);ax.grid(axis='x',alpha=.18)
    config_axis(axs[0],[n+('*' if n.startswith('VoxelMap-') else '') for n in names]);axs[2].set_xscale('log');axs[2].set_xticks([2,10,100,300]);axs[2].set_xticklabels(['2','10','100','300'])
    fig.subplots_adjust(left=.23,right=.99,top=.94,bottom=.10,wspace=.18)
    fig.text(.58,.015,'(a,b) Mean bar; observed peak tick, same recording / 1 s resource sampling. * Mapping core only.\n(c) Three timing replays: filled = mean, open = P95. VoxelMap timing includes live preparation.',ha='center',fontsize=7)
    save(fig,'04_computation_combined')
    pose=csvread(C/'poses_and_stability_per_run.csv');fig,axs=plt.subplots(1,2,figsize=(7,4.6),sharey=True)
    for i,n in enumerate(names):
        pp=sorted([p for p in pose if p['name']==n],key=lambda x:x['bag_id']);yy=i+np.linspace(-.15,.15,len(pp))
        axs[0].scatter([float(p['record_rate_mean_hz']) for p in pp],yy,s=11,color='#297ca2',linewidths=0);axs[1].scatter([1000*float(p['record_gap_max_s']) for p in pp],yy,s=11,color='#b36d36',linewidths=0)
    config_axis(axs[0],names);axs[0].set_xscale('log');axs[0].set_xticks([5,10,100]);axs[0].set_xticklabels(['5','10','100']);axs[0].set_title('(a) Valid pose output rate (Hz)');axs[1].set_title('(b) Largest arrival gap / replay (ms)');axs[1].grid(axis='x',alpha=.18)
    fig.subplots_adjust(left=.23,right=.99,top=.93,bottom=.12,wspace=.20)
    fig.text(.58,.02,'Ten points per configuration, one per recording. Arrival times use simulated ROS time.\nRate and gaps do not measure position accuracy or wall-clock processing latency.',ha='center',fontsize=7)
    save(fig,'05_stability_combined')
    manifest=dict(configurations=names,representative_bag=REP,figure_count=5,files=[dict(file=str(p.relative_to(O)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sorted(F.iterdir())],sampling='Plots subsample retained map points for display only; all numeric statistics retain their original definitions.')
    (F/'figure_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print('Built five combined figures, PNG and PDF.',flush=True)
if __name__=='__main__':main()
