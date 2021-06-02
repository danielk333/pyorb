'''
Orbit transfer
===============
'''

import pyorb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

init = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True, a=1*pyorb.AU, e=0, i=0, omega=0, Omega=0, anom=0)
target = pyorb.Orbit(M0 = pyorb.M_sol, degrees=True, a=2*pyorb.AU, e=0, i=23, omega=0, Omega=90, anom=0)

sliders = []

def plot_transfer(init, target):
    
    def gen_r(dv1, dv2):
        r = np.empty((3,600))
        ind = 0

        orb = init.copy()
        orb.vy += dv1
        for i in range(100):
            orb.anom += 90.0/100.0
            r[:,ind] = orb.r[:,0]
            ind += 1
        orb.vz += dv2
        for i in range(100):
            orb.anom += 90.0/100.0
            r[:,ind] = orb.r[:,0]
            ind += 1
        orb.vy -= dv1
        for i in range(400):
            orb.anom += 90.0/100.0
            r[:,ind] = orb.r[:,0]
            ind += 1
        return r

    r = gen_r(5e3, 5e3)

    rr = np.empty((3,400))
    ind = 0
    ref = target.copy()
    for i in range(400):
        ref.anom += 90.0/100.0
        rr[:,ind] = ref.r[:,0]
        ind += 1
    fig = plt.figure(figsize=(15,15))
    ax = fig.add_subplot(111, projection='3d')
    plt.subplots_adjust(left=0.25, bottom=0.25)
    l, = ax.plot(r[0,:], r[1,:], r[2,:],  '-b')
    ax.plot(rr[0,:], rr[1,:], rr[2,:],  '-r')
    ax.set_title('Orbits', fontsize=22)
    ax.set_xlabel('X-position $x$ [m]', fontsize=20)
    ax.set_ylabel('Y-position $y$ [m]', fontsize=20)
    ax.set_zlabel('Z-position $z$ [m]', fontsize=20)

    axcolor = 'lightgoldenrodyellow'
    axdv1 = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)
    axdv2 = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    sdv1 = Slider(axdv1, 'dV1', 1e3, 50e3, valinit=5e3)
    sdv2 = Slider(axdv2, 'dV2', 1e3, 50e3, valinit=5e3)

    sliders.append(sdv1)
    sliders.append(sdv2)

    def update(val):
        r = gen_r(sdv1.val, sdv2.val)
        l.set_xdata(r[0,:])
        l.set_ydata(r[1,:])
        l.set_3d_properties(r[2,:])
        fig.canvas.draw_idle()

    sdv1.on_changed(update)
    sdv2.on_changed(update)



def transfer(dv1, dv2, init):
    orb = init.copy()
    orb.vy += dv1
    orb.anom += 90
    orb.vz += dv2
    orb.anom += 90
    orb.vy -= dv1
    return orb

def optim(params):
    orb = transfer(params[0], params[1], init)
    d = ((orb.kepler[0] - target.kepler[0])/pyorb.AU)**2
    d += (orb.kepler[1] - target.kepler[1])**2
    d += np.sum((np.radians(orb.kepler[2:5] - target.kepler[2:5]))**2)
    return np.sqrt(d)

v1 = np.linspace(1e3, 10e3, num=100)
v2 = np.linspace(1e3, 10e3, num=100)

plot_transfer(init, target)
plt.show()