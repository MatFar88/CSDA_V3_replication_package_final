"""Deterministic conceptual redraw of the three geometries on RSS 2026 slide 7.
No simulated or empirical observations are represented. Existing empirical
figures are carried over unchanged from the supplied V7.12 package.
Run from any directory; requires numpy and matplotlib.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

root = Path(__file__).resolve().parents[1]
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':12,
                     'pdf.fonttype':42,'ps.fonttype':42})
centers=np.array([[-1.45,-1.05],[1.45,-1.05],[0.,1.35]])
offsets=np.array([[-.22,-.11],[-.13,.18],[.08,.28],[.24,.04],[.07,-.22],[-.10,-.26]])
colors=['#4263AD','#D17825','#35855B']
fig,axes=plt.subplots(1,3,figsize=(8.4,3.35))
fig.subplots_adjust(left=.025,right=.98,bottom=.25,top=.83,wspace=.12)
for j,ax in enumerate(axes):
    for mu,col in zip(centers,colors):
        xy=mu+offsets
        ax.scatter(xy[:,0],xy[:,1],s=24,c=col,edgecolors='white',linewidths=.3,zorder=3)
    ax.set(xlim=(-3.0,3.0),ylim=(-2.45,2.6),xticks=[],yticks=[],aspect='equal')
    for spine in ax.spines.values():spine.set_color('#D0D4D8')
    ax.set_title(['(a) Clean','(b) Bridge','(c) Diffuse leverage'][j],fontsize=12,pad=11)
bridges=np.array([[-.25,-.94],[.2,-1.1],[-.9,.12],[-.6,.35],[.57,.23],[.85,-.02]])
leverage=np.array([[-2.6,2.05],[2.5,2.15],[-2.55,-2.05],[2.55,-1.95],[2.62,.38],[-2.5,.48]])
for ax,xy in [(axes[1],bridges),(axes[2],leverage)]:
    ax.scatter(xy[:,0],xy[:,1],s=34,c='#BA2632',marker='x',linewidths=1.5,zorder=4)
for ax,label in zip(axes,['Genuine observations only','Ambiguity between groups','Large ambient perturbations']):
    ax.text(.5,-.13,label,transform=ax.transAxes,ha='center',va='top',fontsize=10)
fig.legend(handles=[Line2D([0],[0],marker='o',linestyle='',color='#6E7681',markersize=5,label='Genuine groups (colors)'),Line2D([0],[0],marker='x',linestyle='',color='#BA2632',markersize=5,label='Contaminants')],loc='lower center',ncol=2,frameon=False,fontsize=10,bbox_to_anchor=(.5,.005))
fig.savefig(root/'figures/fig_contamination_geometries.pdf',bbox_inches='tight',pad_inches=.03)
fig.savefig(root/'figures/fig_contamination_geometries.png',dpi=180,bbox_inches='tight',pad_inches=.03)
plt.close(fig)
