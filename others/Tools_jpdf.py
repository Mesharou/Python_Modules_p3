import numpy as np
import matplotlib.pyplot as plt

binx = np.arange(-6,6.05,0.05)
biny = np.arange(0,6.025,0.025)


#biny_size = 0.1
#binx_size = biny_size
#binx = np.arange(-6,6+binx_size,binx_size)
#biny = np.arange(0,6+biny_size,biny_size)

xedges = (binx[1:]-binx[:-1]) / 2
yedges = (biny[1:]-biny[:-1]) / 2



def add_contour_per(H_sum_tot,proba_list,c):
    #proba_list = [99.99]
    n = 0

    level = np.zeros(len(proba_list))
    p = np.zeros(len(proba_list))

    for proba in proba_list:

        i = dichotomie(f,-4,6,0.01,proba,H_sum_tot)
        #print(i)

        H_tot = np.sum(H_sum_tot)

        H_filter = np.where(H_sum_tot > 10**i, H_sum_tot,0)
        p[n] = np.sum(H_filter)*100/H_tot

        level[n] = 10**i
        n=n+1

    fmt = {}
    for l, s in zip(level, p):
        fmt[l] = str(round(s))+'%'
        #fmt[l] = ''

    #CS = plt.contour(xedges,yedges,H_sum_tot,level,linewidths=2,alpha=1,colors=c,linestyles='-')
    #plt.clabel(CS,level, inline=1, fontsize=18,fmt = fmt)
    
    return level,fmt

def f(x,p,H_sum_tot):
    
    H_filter = np.where(H_sum_tot > 10**x, H_sum_tot,0)
    H_tot = np.sum(H_sum_tot)
    return(np.sum(H_filter)*100/H_tot-p) 

def dichotomie(f,a,b,e,p,H_sum_tot):
    delta = 1

    while delta > e:
        m = a + (b - a) / 2
        delta = abs(b - a)
        #print("{:15} , {:15} , {:15} , {:15} , {:15} , {:15} , {:15} ".format(a,b,m,f(a),f(b),f(m),delta) )
        if f(m,p,H_sum_tot) == 0:
            return m
        elif f(a,p,H_sum_tot) * f(m,p,H_sum_tot)  > 0:
            a = m
        else:
            b = m
    return m


def add_contour_per_wc(H_sum_tot,proba_list,c):
    
    n = 0

    level = np.zeros(len(proba_list))
    p = np.zeros(len(proba_list))

    
    H_sum_tot_neg = abs(np.where(H_sum_tot < 0, H_sum_tot,0))
    H_sum_tot_pos = abs(np.where(H_sum_tot > 0, H_sum_tot,0))


    for proba in proba_list:


        i = dichotomie(f,-10,10,0.001,proba,H_sum_tot_neg)
        #print(i)

        H_tot = np.sum(H_sum_tot_neg)

        H_filter = np.where(H_sum_tot_neg > 10**i, H_sum_tot_neg,0)
        p[n] = np.sum(H_filter)*100/H_tot

        level[n] = 10**i
        n=n+1

    fmt = {}
    for l, s in zip(level, p):
        #fmt[l] = str(round(s,1))+'%'
        fmt[l] = ''
    
    CS = plt.contour(xedges[1:],yedges[1:],H_sum_tot_neg,level,linewidths=2,alpha=0.8,colors=c)
    ax.clabel(CS,level, inline=1, fontsize=16,fmt = fmt)
    
    n = 0
    level = np.zeros(len(proba_list))
    p = np.zeros(len(proba_list))

    return


from matplotlib.cm import get_cmap
from cycler import cycler
from matplotlib.pyplot import cm

def compute_fraction(H,H_wc):
        
    H_wc_neg = np.where(H_wc<0,H_wc,0)

    xx, yy = np.meshgrid(xedges[1:], yedges[1:])
    #H_80 £*= np.where(H < 10**x_80, H,0)
    H_f = np.where(xx > 0.5, H,0)
    H_f = np.where(yy > xx, H_f,0)
    
    H_f_wc = np.where(xx > 0.5, H_wc,0)
    H_f_wc = np.where(yy > xx, H_f_wc,0)

    H_f_per = np.sum(H_f)*100/np.sum(H)
    H_f_wc_per = np.sum(H_f_wc)*100/np.sum(H_wc_neg)
    
    return(H_f_per,H_f_wc_per)


def compute_fraction_sf(H,H_wc):
        
    H_wc_neg = np.where(H_wc<0,H_wc,0)

    xx, yy = np.meshgrid(xedges[1:], yedges[1:])
    #H_80 = np.where(H < 10**x_80, H,0)
    H_f = np.where(xx > 0.5, H,0)
    H_f = np.where(yy > xx, H_f,0)
    
    H_f_wc = np.where(xx > 0.5, H_wc,0)
    H_f_wc = np.where(yy > xx, H_f_wc,0)

    H_f_per = np.sum(H_f)*100/np.sum(H)
    H_f_wc_per = np.sum(H_f_wc)*100/np.sum(H_wc_neg)
    
    return(H_f_per,H_f_wc_per)


def compute_fraction_sf_ac(H,H_wc):
        
    H_wc_neg = np.where(H_wc<0,H_wc,0)

    xx, yy = np.meshgrid(xedges[1:], yedges[1:])
    #H_80 = np.where(H < 10**x_80, H,0)
    H_f = np.where(xx < -0.5, H,0)
    H_f = np.where(yy > -xx, H_f,0)
    
    H_f_wc = np.where(xx < -0.5, H_wc,0)
    H_f_wc = np.where(yy > -xx, H_f_wc,0)

    H_f_per = np.sum(H_f)*100/np.sum(H)
    H_f_wc_per = np.sum(H_f_wc)*100/np.sum(H_wc_neg)
    
    return(H_f_per,H_f_wc_per)

def compute_fraction_sf_wc(H_wc):
        
    H_wc_neg = np.where(H_wc<0,H_wc,0)

    xx, yy = np.meshgrid(xedges[1:], yedges[1:])
    #H_80 = np.where(H < 10**x_80, H,0)
    H_f = np.where(xx < -0.5, H,0)
    H_f = np.where(yy > -xx, H_f,0)
    
    H_f_wc = np.where(xx < -0.5, H_wc,0)
    H_f_wc = np.where(yy > -xx, H_f_wc,0)

    H_f_per = np.sum(H_f)*100/np.sum(H)
    H_f_wc_per = np.sum(H_f_wc)*100/np.sum(H_wc_neg)
    
    return(H_f_per,H_f_wc_per)


def add_contour_front_h(H_sum_tot,color,per,proba_list):
    
    n = 0

    level = np.zeros(len(proba_list))
    p = np.zeros(len(proba_list))



    for proba in proba_list:


        i = dichotomie(f,-8,2,0.001,proba,H_sum_tot)
        #print(i)

        H_tot = np.sum(H_sum_tot)

        H_filter = np.where(H_sum_tot > 10**i, H_sum_tot,0)
        p[n] = np.sum(H_filter)*100/H_tot

        level[n] = 10**i
        n=n+1

    fmt = {}
    #for l, s in zip(level, p):
    #    fmt[l] = str(ceil(per))+'%'
        #fmt[l] = '3 %'
    
    CS = plt.contour(xedges[:],yedges[:],H_sum_tot,level,linewidths=2,alpha=0.8,colors=color)
    plt.contourf(xedges[:],yedges[:],H_sum_tot,level,colors='none',
                  hatches=['xx'],
                  extend='lower')
    
    #ax.clabel(CS,level, inline=2, fontsize=26,fmt = fmt)
    
    n = 0
    level = np.zeros(len(proba_list))
    p = np.zeros(len(proba_list))

    return
