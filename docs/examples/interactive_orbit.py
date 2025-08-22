'''
Interactive orbit
===================
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

import pyorb

num = 1000
ax_lims = 2
asymtote_limit = 0.99

orb = pyorb.Orbit(
    M0 = 1.0,
    G = pyorb.get_G(length='AU', mass='Msol', time='y'),
    num = num,
    a = 1.0, 
    e = 0, 
    i = 0, 
    omega = 0, 
    Omega = 0, 
    anom = np.linspace(0, 360, num=num),
    degrees = True,
    type = 'true',
)

# for some optimization
orb.direct_update = False


def add_vel(ax):
    r = orb._cart[:3, 0]
    v = orb._cart[3:, 0]
    vel = ax.quiver(
        r[0], r[1], r[2], 
        v[0], v[1], v[2], 
        length=ax_lims*0.05,
    )
    return vel


def draw():
    r = orb.r

    l.set_xdata(r[0, 1:])
    l.set_ydata(r[1, 1:])
    l.set_3d_properties(r[2, 1:])

    dot.set_xdata([r[0, 0]])
    dot.set_ydata([r[1, 0]])
    dot.set_3d_properties([r[2, 0]])

    global vel
    vel.remove()
    vel = add_vel(ax)

    fig.canvas.draw_idle()


def update_orb(val):
    orb.a = s_a.val
    orb.i = s_i.val
    orb.omega = s_omega.val
    orb.Omega = s_Omega.val
    orb._kep[5, 0] = s_nu.val
    draw()


r = orb.r

fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')
plt.subplots_adjust(left=0.25, bottom=0.25)
l, = ax.plot(r[0, :], r[1, :], r[2, :], '-b')
dot, = ax.plot([r[0, 0]], [r[1, 0]], [r[2, 0]], 'ob')
vel = add_vel(ax)
ax.plot([0], [0], [0], 'or')

axc_b = fig.add_axes([0.15, 0.05, 0.05, 0.02])
axp_b = fig.add_axes([0.15, 0.10, 0.05, 0.02])
axh_b = fig.add_axes([0.15, 0.15, 0.05, 0.02])
axr_b = fig.add_axes([0.05, 0.10, 0.05, 0.02])
c_b = Button(axc_b, 'Circular')
p_b = Button(axp_b, 'Parabolic')
h_b = Button(axh_b, 'Hyperbolic')
r_b = Button(axr_b, 'Reset')

ax.set_title('Orbit', fontsize=22)
ax.set_xlabel('X-position [AU]', fontsize=20)
ax.set_ylabel('Y-position [AU]', fontsize=20)
ax.set_zlabel('Z-position [AU]', fontsize=20)

ax.set_xlim([-ax_lims, ax_lims])
ax.set_ylim([-ax_lims, ax_lims])
ax.set_zlim([-ax_lims, ax_lims])

axcolor = 'lightgoldenrodyellow'
ax_a = plt.axes([0.25, 0.05, 0.2, 0.03], facecolor=axcolor)
ax_e = plt.axes([0.25, 0.1, 0.2, 0.03], facecolor=axcolor)
ax_i = plt.axes([0.25, 0.15, 0.2, 0.03], facecolor=axcolor)

ax_omega = plt.axes([0.6, 0.05, 0.2, 0.03], facecolor=axcolor)
ax_Omega = plt.axes([0.6, 0.1, 0.2, 0.03], facecolor=axcolor)
ax_nu = plt.axes([0.6, 0.15, 0.2, 0.03], facecolor=axcolor)

s_a = Slider(ax_a, 'a [AU]', 0, 2, valinit=1)
s_e = Slider(ax_e, 'e [1]', 0, 2, valinit=0)
s_e.is_hyp = False
s_i = Slider(ax_i, 'i [deg]', 0, 180, valinit=0)

s_omega = Slider(ax_omega, 'omega [deg]', 0, 360, valinit=0)
s_Omega = Slider(ax_Omega, 'Omega [deg]', 0, 360, valinit=0)
s_nu = Slider(ax_nu, 'nu [deg]', 0, 360, valinit=0)

sliders = [s_a, s_i, s_omega, s_Omega, s_nu]


def set_state(event, source):
    if source == 'Circular':
        s_e.set_val(0)
    elif source == 'Parabolic':
        s_e.set_val(1)
    elif source == 'Hyperbolic':
        s_e.set_val(2)
    elif source == 'Reset':
        s_a.set_val(1)
        s_e.set_val(0)
        s_i.set_val(0)
        s_omega.set_val(0)
        s_Omega.set_val(0)
        s_nu.set_val(0)
    update_orb(None)


c_b.on_clicked(lambda x: set_state(x, 'Circular'))
p_b.on_clicked(lambda x: set_state(x, 'Parabolic'))
h_b.on_clicked(lambda x: set_state(x, 'Hyperbolic'))
r_b.on_clicked(lambda x: set_state(x, 'Reset'))


def update_ecc(val):
    orb.e = s_e.val

    if s_e.val > 1:
        s_e.is_hyp = True

        theta_inf = pyorb.kepler.true_of_the_asymptote(s_e.val, degrees=True)
        orb._kep[5, :] = np.linspace(
            -asymtote_limit*theta_inf, 
            asymtote_limit*theta_inf, 
            num=num,
        )

        s_nu.set_val(np.mod(s_nu.val, asymtote_limit*theta_inf))
        s_nu.valmin = -asymtote_limit*theta_inf
        s_nu.valmax = asymtote_limit*theta_inf
        s_nu.ax.set_xlim(s_nu.valmin, s_nu.valmax)
        orb._kep[5, 0] = s_nu.val

    elif s_e.val == 1:
        s_e.is_hyp = True
        orb._kep[5, :] = np.linspace(-180, 180, num=num)

        s_nu.set_val(np.mod(s_nu.val, 180))
        s_nu.valmin = -180
        s_nu.valmax = 180
        s_nu.ax.set_xlim(s_nu.valmin, s_nu.valmax)
        orb._kep[5, 0] = s_nu.val
    else:
        if s_e.is_hyp:
            orb._kep[5, :] = np.linspace(0, 360, num=num)
            s_nu.set_val(np.mod(s_nu.val + 360, 360))
            s_nu.valmin = 0
            s_nu.valmax = 360
            s_nu.ax.set_xlim(s_nu.valmin, s_nu.valmax)    
            orb._kep[5, 0] = s_nu.val
        s_e.was_hyp = False
    draw()


for slide in sliders:
    slide.on_changed(update_orb)

s_e.on_changed(update_ecc)

plt.show()
