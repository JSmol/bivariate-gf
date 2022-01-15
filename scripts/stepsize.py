from sage.all import *
from helpers import *

# WARNING: ERROR TERM NEEDS TO BE WORKED OUT!
ERROR_TERM = 0.05
def step_params(H, x, y, x0, y0):
  yunivariate = PolynomialRing(BASEFIELD, y)(H(x=x0))
  yroots = [r[0] for r in yunivariate.roots()]
  yroots = list(filter(lambda r: not zero(r - y0), yroots))
  DyDx = - diff(H, x) / diff(H, y)

  xunivariate = PolynomialRing(BASEFIELD, x)(H(y=0))
  xroots = [r[0] for r in xunivariate.roots()]
  delta = min(1, abs(x0) / 2, *[abs(x0 - r) / 2 for r in xroots])
  while delta > ZERO:

    valid = [
      abs(y0 - y1) / 2 > abs(delta * DyDx(x=x0, y=y1) + ERROR_TERM)
      for y1 in yroots
    ]

    if all(valid):
      eps = min(abs(y0 - y1) / 2 + ZERO for y1 in yroots)
      return delta, eps
    else: delta /= 2

  # TODO: deal with extremely small delta better, possibly return a FAIL signal?
  return 0, ZERO

def step(H, x, y, x0, y0, direction):
  print(f'-- taking step from {x0, y0} --')
  delta, eps = step_params(H, x, y, x0, y0)
  x1 = x0 + delta * direction

  univariate = PolynomialRing(BASEFIELD, y)(H(x=x1))
  yroots = [r[0] for r in univariate.roots()]
  yroots = list(filter(lambda r: abs(y0 - r) <= eps, yroots))

  print(f'{delta, eps = }')
  print(f'{yroots = }')

  assert(len(yroots) == 1)
  return x1, yroots[0]

# TESTING:
if __name__ == '__main__':

  STEPS = 10
  MAX_HEIGHT = 10
  def walk(H, x, y, r, s, x0, y0):
    DyDx = - diff(H, x) / diff(H, y)
    Dh = - r/x - (s/y)*DyDx
    path = []
    for i in range(STEPS):
      direction = conjugate(Dh(x=x0, y=y0))
      direction /= abs(direction)
      x0, y0 = step(H, x, y, x0, y0, direction)
      path.append((x0, y0))
    return path

  R = PolynomialRing(BASEFIELD, 'x, y')
  x, y = R.gens()
  r, s = 1, 1
  H = 1 - x**2 - y**2

  C, NS = get_points(H, x, y, r, s)
  print(*C, sep='\n')

  for p in C:

    Dh, degree = degeneracy(H, p, x, y, r, s)
    print(f'{p = }')
    print(f'{Dh = }')
    print(f'{degree = }')

    # x**(1/n) does not return an algebraiac number in sage
    dir0 = conjugate(Dh(x=p[x], y=p[y]))
    dir0 = (PolynomialRing(BASEFIELD, x)(x**degree - dir0)).roots()[0][0]
    dir0 /= abs(dir0)

    for k in range(degree):

      rotate = exp(k*2*pi*i/degree)
      direction = rotate * dir0
      print(f'{direction = }')

      x1, y1 = step(H, x, y, p[x], p[y], direction)
      print(f'{x1, y1 = }')

      steps = walk(H, x, y, r, s, x1, y1)
      
