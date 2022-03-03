'''
Orbit transfer
===============
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import pyorb

num = 1000

orb = pyorb.Orbit(
    M0 = 1.0,
    G = pyorb.get_G(length='AU', mass='Msol', time='y'),
    num = num,
    a = 1.0, 
    e = 0, 
    i = 0, 
    omega = 0, 
    Omega = 0, 
    anom = np.linspace(0,360,num=num),
    degrees = True,
    type = 'true',
)

#for some optimization
orb.direct_update = False


r = orb.r

fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.25)
l, = ax.plot(r[0,:], r[1,:], r[2,:],  '-b')
dot, = ax.plot([r[0,0]], [r[1,0]], [r[2,0]],  'ob')
ax.plot([0], [0], [0], 'or')
ax.set_title('Orbit', fontsize=22)
ax.set_xlabel('X-position [AU]', fontsize=20)
ax.set_ylabel('Y-position [AU]', fontsize=20)
ax.set_zlabel('Z-position [AU]', fontsize=20)
ax_lims = 2
ax.set_xlim([-ax_lims,ax_lims])
ax.set_ylim([-ax_lims,ax_lims])
ax.set_zlim([-ax_lims,ax_lims])

axcolor = 'lightgoldenrodyellow'
ax_a = plt.axes([0.25, 0.05, 0.2, 0.03], facecolor=axcolor)
ax_e = plt.axes([0.25, 0.1, 0.2, 0.03], facecolor=axcolor)
ax_i = plt.axes([0.25, 0.15, 0.2, 0.03], facecolor=axcolor)

ax_omega = plt.axes([0.6, 0.05, 0.2, 0.03], facecolor=axcolor)
ax_Omega = plt.axes([0.6, 0.1, 0.2, 0.03], facecolor=axcolor)
ax_nu = plt.axes([0.6, 0.15, 0.2, 0.03], facecolor=axcolor)

s_a = Slider(ax_a, 'a [AU]', 0, 2, valinit=1)
s_e = Slider(ax_e, 'e [1]', 0, 2, valinit=0)
s_i = Slider(ax_i, 'i [deg]', 0, 180, valinit=0)

s_omega = Slider(ax_omega, 'omega [deg]', 0, 360, valinit=0)
s_Omega = Slider(ax_Omega, 'Omega [deg]', 0, 360, valinit=0)
s_nu = Slider(ax_nu, 'nu [deg]', 0, 360, valinit=0)

sliders = [s_a, s_e, s_i, s_omega, s_Omega, s_nu]

def update(val):
    orb.a = s_a.val
    orb.e = s_e.val
    orb.i = s_i.val
    orb.omega = s_omega.val
    orb.Omega = s_Omega.val
    orb._kep[5,0] = s_nu.val

    r = orb.r

    l.set_xdata(r[0,1:])
    l.set_ydata(r[1,1:])
    l.set_3d_properties(r[2,1:])

    dot.set_xdata([r[0,0]])
    dot.set_ydata([r[1,0]])
    dot.set_3d_properties([r[2,0]])
    fig.canvas.draw_idle()

for slide in sliders:
    slide.on_changed(update)


plt.show()