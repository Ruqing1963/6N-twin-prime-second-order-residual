#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the two-panel showdown figure (obs / first-order / second-order) from
../data/so_plot_S10.csv (produced by second_order.py with default MAXK=10).
Left panel 6dN=42: second order helps but a two-centre (shield-gain) residual
remains (shaded). Right panel 6dN=210: second order closes the residual.
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/so_plot_S10.csv')))
def series(sixd):
    om=[];o=[];f=[];s=[]
    for r in rows:
        if int(r['6dN'])==sixd:
            om.append(int(r['omega'])); o.append(float(r['obs']))
            f.append(float(r['first_order'])); s.append(float(r['second_order']))
    return np.array(om),np.array(o),np.array(f),np.array(s)
fig,(axL,axR)=plt.subplots(1,2,figsize=(13.5,5.2))
om,o,f,s=series(210)
axR.plot(om,o,'o-',color='#c0392b',lw=2.4,ms=8,label='observed',zorder=4)
axR.plot(om,f,'^--',color='#e8845b',lw=1.6,ms=7,alpha=.8,label='first-order $\\mathfrak{S}_{2,\\omega}$',zorder=3)
axR.plot(om,s,'s-',color='#185FA5',lw=1.8,ms=7,label='second-order (+compression)',zorder=3)
axR.axhline(1,color='gray',ls=':',lw=1)
axR.set_title('$6\\Delta N=210$: second order matches\n(spatial-compression penalty closes the residual)',fontsize=11)
axR.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11); axR.set_ylabel(r'relative preference $r(d\mid\omega)$',fontsize=11)
axR.legend(fontsize=9,loc='lower left'); axR.grid(alpha=.25); axR.set_xticks(range(1,7))
om,o,f,s=series(42)
axL.plot(om,o,'o-',color='#c0392b',lw=2.4,ms=8,label='observed',zorder=4)
axL.plot(om,f,'^--',color='#e8845b',lw=1.6,ms=7,alpha=.8,label='first-order $\\mathfrak{S}_{2,\\omega}$',zorder=3)
axL.plot(om,s,'s-',color='#185FA5',lw=1.8,ms=7,label='second-order (+compression)',zorder=3)
axL.fill_between(om,s,o,where=(o>s),color='#2ca25f',alpha=.16,label='two-centre (shield-gain) residual',zorder=2)
axL.axhline(1,color='gray',ls=':',lw=1)
axL.set_title('$6\\Delta N=42$: second order helps,\nhigh-$\\omega$ shield-gain residual remains',fontsize=11)
axL.set_xlabel(r'$\omega_{>3}(N)$',fontsize=11); axL.set_ylabel(r'relative preference $r(d\mid\omega)$',fontsize=11)
axL.legend(fontsize=9,loc='upper left'); axL.grid(alpha=.25); axL.set_xticks(range(1,7))
plt.suptitle('Observed vs first-order vs second-order conditional singular series in $S_{10}$ (23,988,173 twin centres)',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper4_second_order.pdf',bbox_inches='tight')
plt.savefig('fig_paper4_second_order.png',dpi=160,bbox_inches='tight')
print("figure saved")
