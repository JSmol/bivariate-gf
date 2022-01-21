from sage.all import *
from sage.plot.plot3d.base import Graphics3d

def plot_paths(H, paths, rng=10, color='red'):
  return sum([
  point3d(step, size=10, color=color) if abs(step[0]) < rng and abs(step[1]) < rng else None
    for critical_pnt in paths for path in critical_pnt for step in path
  ])

def mag(z):
  return n(sqrt(real(z)**2 + imag(z)**2))

SMALL = 1e-5
def h(x, y, r, s):
  if abs(x) < SMALL or abs(y) < SMALL:
    return 10000
  return -r*log(abs(x)) -s*log(abs(y))

def plot_params(H, x, y, r, s, rng=4, maxval=4):
  params = solve(SR(H), SR(y))
  return sum([
    plot3d(
      lambda u, v: min(h(u+v*1j, g.rhs()(x=u+v*1j), r, s), maxval),
      # lambda u, v: h(u+v*1j, g.rhs()(x==u+v*1j), r, s),
      urange=(-rng, rng),
      vrange=(-rng, rng),
    )
    for g in params
  ])
