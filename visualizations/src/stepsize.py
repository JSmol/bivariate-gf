import pyglet as pg
import numpy as np
import sympy as sp

from helpers.spaces import ComplexPlane
from helpers.trackers import ValueTracker

def yroots(H, x0):
  coeffs = list(map(complex, sp.Poly(H.subs(x, x0), domain=sp.CC).all_coeffs()))
  return list(map(complex, np.roots(coeffs))) # np is faster

x, y = sp.symbols('x, y')
# H = 1 - x**2 - y**2
# H = x**5*y**2 + 2*x**2*y - 2*x**3*y + 4*y + x - 2
H = (x**2 - x + 1) * y**5 + (x**4+x) * y**4 + (x+1) * y**2 + (x**2 - 1) * y + (x**5+1)

ltc = sp.LC(H, y)
res = sp.resultant(H, sp.diff(H, y), y)

ZERO = 1e-9
Hy = sp.diff(H, y)
Hy = sp.lambdify([x, y], Hy)

DyDx = - sp.diff(H, x) / sp.diff(H, y)
DyDx = sp.lambdify([x, y], DyDx) # this is a lot faster

""" get_delta(x0, y0, roots)
  Try to calculate a radius around x0 such that if |x0 - x1| < delta and si are such that |yi(x0) - yi(x0+delta)| < si we have:
    y0(x0) - yi(x0+delta)| > s0 whenever i != y0
"""
BINARY_SEARCH_ITERS = 12
ERROR_TERM = 0 # should be a function of H, x0, y0, y1, delta, something..... ?
MIN_DELTA = 1e-9
def get_delta(x0, y0, roots):
  delta = 1/2
  for i in range(2, BINARY_SEARCH_ITERS):
    # if every root changes less than the seperation amount, delta is small enough
    if all(abs(y0 - y1) / 2 > abs(delta * DyDx(x0, y1) + ERROR_TERM) for y1 in roots):
      delta += 1 / (2**i)
    else:
      delta -= 1 / (2**i)
  return delta

window = pg.window.Window(fullscreen=True, resizable=False)
HEIGHT = window.height
WIDTH = window.width

batch = pg.graphics.Batch()
plane = ComplexPlane(window, batch)

x0 = ValueTracker(plane)
window.push_handlers(x0.on_mouse_motion)

pg.text.Label(f'{H = }', font_size=16, x=10, y=HEIGHT-(26), color=(200, 50, 50, 255), batch=batch)

pg.text.Label('Resultant(H, diff(H, y)) Roots', font_size=16, x=10, y=10+26, color=(200, 50, 200, 255), batch=batch)
res_circles = [
  pg.shapes.Circle(x=plane.c2p(complex(r))[0], y=plane.c2p(complex(r))[1], radius=8, color=(200, 50, 200), batch=batch)
  for r in sp.roots(res)
]

pg.text.Label('Leading Term Roots', font_size=16, x=10, y=10, color=(200, 200, 50, 255), batch=batch)
ltc_circles = [
  pg.shapes.Circle(x=plane.c2p(complex(r))[0], y=plane.c2p(complex(r))[1], radius=5, color=(200, 200, 50), batch=batch)
  for r in sp.roots(ltc)
]

MAXRADIUS = 200
@window.event
def on_draw():

  window.clear()
  batch.draw()

  x0pos = plane.c2p(x0.z)
  pg.shapes.Circle(x=x0pos[0], y=x0pos[1], radius=5, color=(50, 50, 200)).draw()
  pg.text.Label(f'{x0.z:.2f}', font_size=16, x=x0pos[0]+10, y=x0pos[1]+10, color=(50, 50, 200, 255)).draw()

  roots = yroots(H, x0.z)
  deltas = []
  for r in roots:

    pnt = plane.c2p(r)
    pg.shapes.Circle(x=pnt[0], y=pnt[1], radius=10, color=(50, 200, 50)).draw()

    # dont divide by 0
    if abs(Hy(x0.z, r)) < ZERO: continue

    deltas.append(get_delta(x0.z, r, list(filter(lambda x: x != r, roots))))

  delta = min([1, *deltas])
  pg.text.Label(f'{delta = }', font_size=16, x=10, y=HEIGHT - (2*26), color=(200, 50, 50, 255)).draw()
  pg.shapes.Arc(x=x0pos[0], y=x0pos[1], radius=min(MAXRADIUS, plane.scale * delta), color=(200, 200, 200)).draw()

  # draw circles around each root
  for r in roots:

    # non-smooth point, or very close to one
    if abs(Hy(x0.z, r)) < ZERO: continue

    pnt = plane.c2p(r)
    ych = abs(delta * DyDx(x0.z, r) + ERROR_TERM)
    pg.shapes.Arc(x=pnt[0], y=pnt[1], radius=min(MAXRADIUS, plane.scale * ych), color=(50, 200, 200)).draw()

###########################################################################################################################################
@window.event
def on_key_press(symbol, _):
  if symbol == 113: # press q to close window and exit python
    window.close()
    exit()

pg.app.run()

