import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from numpy import exp
from matplotlib import cm

from scipy.optimize import curve_fit
from scipy.integrate import trapz, cumtrapz

# get list of directories
cwd = os.getcwd()
dlist = np.sort(glob.glob('sim_*'))

# allocate for data
T = []
C = []
IF = []
TB = []
FRC = []

# grab data in each study subfloder
for d in dlist:

    # laod primary results
    print('...loading data in study '+d)
    res = np.genfromtxt(cwd+'/'+d+'/data/primary_results')

    T.append(res[:,0])
    C.append(res[:,1])

    # load time step resolution
    iF = np.genfromtxt(cwd+'/'+d+'/data/inletFlow')
    IF.append(iF)

    # load breath periods
    TBTVN = np.genfromtxt(cwd+'/'+d+'/data/TBTVN')
    period = TBTVN[:,0]
    TB.append(period)

    # il = np.genfromtxt(cwd+'/'+d+'/data/inputLung', comments='//')
    il = np.genfromtxt(cwd+'/'+d+'/constant/systemProperties', comments='//')
    FRC.append(il[0,1]*1e03) # in liter

T  = np.array(T)
C  = np.array(C)
IF = np.array(IF)
TB = np.array(TB)
FRC = np.array(FRC)

# mask for end-exp. concentration
end_exp_t = 10 # time substraction for end-exp. concentration
M = []
for tb, t in zip(TB, T):
    mask = np.zeros_like(t).astype(np.bool)
    tb_c = 0
    for p in tb:
        tb_c += p
        mask += (1e03*t).astype(np.int32) == (1e03*tb_c-end_exp_t).astype(np.int32)

    M.append(mask)

M = np.array(M)

# volume, cumulative expired volume, turn over
CEV = []
FRC_F = []
VOL = []
TO = []

for t, c, iF, frc in zip(T, C, IF, FRC):

    vol = cumtrapz(iF, t, initial=0)
    VOL.append(vol)

    exp_flow = iF.copy()
    exp_flow[iF>0] = 0.

    cev = cumtrapz(-exp_flow, t, initial=0)*1e03 # in liter
    frc_f = cumtrapz(-exp_flow * c, t, initial=0)*1e03  # in liter

    CEV.append(cev)
    FRC_F.append(frc_f)
    TO.append(cev/frc)

VOL = np.array(VOL)
CEV = np.array(CEV)
FRC_F = np.array(FRC_F)
TO  = np.array(TO)

print('initial FRC = ' + str(FRC))
print('actual FRC = ' + str (FRC_F))

# plot - FRC
lw = 0.7; step = 5000; k = 0
fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(111)

for t, f in zip(T, FRC_F):
    plt.text(t[-1], f[-1], round(FRC[k],2), {'fontsize': 8, 'ha': 'center', 'va': 'center'})
    plt.text(t[-1]-20, f[-1], round(f[-1],2), {'fontsize': 8, 'ha': 'center', 'va': 'center'})
    ax.plot(t[::step], f[::step], lw=lw, color=cm.jet(k*100), linewidth=2, label=f[-1])
    k += 1

ax.set_xlabel('time')
ax.set_ylabel('FRC')
plt.savefig('FRC_vs_time.pdf')
plt.show()

# lung clearance index (LCI)
threshold = 0.025
for t, tb, c, cev, m, frc in zip(T, TB, C, CEV, M, FRC):
    end_exp_c = c[m]
    LCI_ind = -1
    for k in range(tb.size-3):
        eec    = end_exp_c[k]
        if (eec < threshold):
            LCI_ind = k
            break

    if LCI_ind > 0:
        LCI_time = np.cumsum(tb)[LCI_ind] # in seconds
        LCI = cev[t > LCI_time - 1e-03*end_exp_t].min()/frc
        LCI_F = cev[t > LCI_time - 1e-03*end_exp_t].min()/max(frc_f)
        # print('LCI = '+str(LCI))
    else:
        print('Threshold concentration for LCI calculation not reached in given washout.')

print('initial LCI = ' + str(LCI))
print('actual LCI = ' + str(LCI_F))

# slope phase III
def S_III(cev, cN2, bounds):

    # lower and upper bounds
    lb, ub = bounds
    ev = cev - cev[0]
    mask = (ev > lb*ev[-1]) & (ev < ub*ev[-1])

    # slope
    p = np.polyfit(cev[mask], cN2[mask], 1)
    slope = p[0]

    # normalized slope
    slope_norm = slope/cN2[mask].mean()

    return slope, slope_norm, p

SIII = []
P = []
for t, tb, c, cev in zip(T, TB, C, CEV):
    cN2 = 80.*c
    p = []
    sIII = []
    tb_c = 0
    for tb_i in tb:
        mask = (t >= tb_c) & (t < tb_c+tb_i)
        sIII_i, snIII_i, p_i = S_III(cev[mask], cN2[mask], [0.50, 0.90])
        sIII.append([sIII_i, snIII_i])
        p.append(p_i)

        tb_c += tb_i

    SIII.append(sIII)
    P.append(p)


SIII = np.array(SIII)
P = np.array(P)

# Scond, Sacin
Scond = []
Sacin = []
for to, sIII, m in zip(TO, SIII, M):
    to_ee = to[m]
    mask = (to_ee > 1.5 ) & (to_ee < 6.0)
    p = np.polyfit(to_ee[mask], sIII[mask,1], 1)

    s_cond = p[0]
    s_acin = sIII[:3,1].mean() - s_cond*to_ee[:3].mean()

    Scond.append(s_cond)
    Sacin.append(s_acin)

# plot - LCI
fig = plt.figure(figsize=(8, 5))
L = ['$\sigma=0.1$','$\sigma=0.2$','$\sigma=0.4$']
#L = ['$\hat\sigma=0.05$','$\hat\sigma=0.25$','$\hat\sigma=0.5$']

ax = fig.add_subplot(111)
k = 0
for t, c in zip(T, C):
    cN2 = 80.*c
    ax.plot(t, cN2, color=cm.jet(k*100), lw=lw, label=L[k])
    k += 1
ax.plot(t, threshold*80*np.ones_like(t), '--k', lw=lw, label='LCI threshold')

ax.set_xlabel('time $t$ $[s]$')
ax.set_ylabel('$N_2$   $c/c_{max}$')
plt.legend(fontsize=8)

plt.savefig('phi_global_LN_LCI.pdf')

# plot - Scond
fig = plt.figure(figsize=(8, 5))

ax = fig.add_subplot(111)
k = 0
for t, tb, cev, c, p in zip(T, TB, CEV, C, P):
    cN2 = 80.*c
    ax.plot(cev, cN2, color=cm.jet(k*100), lw=lw, label=L[k])
    tb_c = 0
    for tb_i, p_i in zip(tb, p):
        mask = (t >= tb_c) & (t < tb_c+tb_i)
        f = np.polyval(p_i, cev[mask])
        ax.plot(cev[mask], f, color=cm.jet(k*100), ls='-.', lw=lw)
        tb_c += tb_i

    k += 1

ax.set_xlabel('$CEV$   $[l]$')
ax.set_ylabel('$N_2$   $c/c_{max}$')
plt.legend(fontsize=8)

plt.savefig('phi_global_LN_SIII.pdf')

# plot - LCI
fig = plt.figure(figsize=(8, 5))
lw = 0.7

ax = fig.add_subplot(111, ylim=(0, 1))
k = 0
for to, sIII, m in zip(TO, SIII, M):
    to_ee = to[m]
    mask = (to_ee > 1.5 ) & (to_ee < 6.0)
    p = np.polyfit(to_ee[mask], sIII[mask,1], 1)
    f = np.polyval(p, to_ee)
    ax.plot(to_ee, sIII[:,1], marker='s', ms=2, color=cm.jet(k*100), linestyle='none', label=L[k])
    ax.plot(to_ee, f, color=cm.jet(k*100), ls='-.', lw=lw)

    k += 1

ax.set_xlabel('turnover $TO$')
ax.set_ylabel('${S_n}_{III}$')
plt.legend(fontsize=8)

plt.savefig('phi_global_LN_Scond.pdf')
plt.show()
