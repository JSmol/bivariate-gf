from sage.all import *
from helpers import get_points, degeneracy, h

SMALL = 1e-14
def step(H, x, y, r, s, x0, y0, direction):
  yroots = H(x=x0).univariate_polynomial().roots(multiplicities=False)
  yroots = list(filter(lambda r: r != y0, yroots))
  # print(y0, yroots)

  delta = QQbar(min(1, abs(x0)/2))
  while True:
    while True:
      try:
        dx = CBF(x0).add_error(delta)
        yroots = H(x=dx).univariate_polynomial().roots(multiplicities=False)
      except ValueError:
        delta /= 2
      else:
        break

      if delta < SMALL:
        assert(False, 'delta is very small')

    # print(f'{delta, yroots = }')

    # NOTE: y0 is only in 1 ball
    assert(len(list(filter(lambda r: y0 in r, yroots))) == 1)

    y1ball = list(filter(lambda r: y0 in r, yroots))[0]

    x1 = QQbar(x0 + delta * direction)
    algebraic_yroots = H(x=x1).univariate_polynomial().roots(multiplicities=False)

    # NOTE: only 1 root is in y1ball
    assert(len(list(filter(lambda r: r in y1ball, algebraic_yroots))) == 1)

    y1 = list(filter(lambda r: r in y1ball, algebraic_yroots))[0]

    # NOTE: y1 is only in 1 ball
    assert(list(filter(lambda r: y1 in r, yroots)))

    if RR(h(x0, y0, r, s)) > RR(h(x1, y1, r, s)):
      # print(RR(h(x0, y0, r, s)), RR(h(x1, y1, r, s)))
      print('Not Ascending')
      delta /= 2
    else:
      break

    print(f'{delta = }')
    break

  # returned 2 roots in the base field
  return x1, y1

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
      x0, y0 = step(H, x, y, r, s, x0, y0, direction)
      path.append((x0, y0))
    return path

  R = PolynomialRing(QQbar, 'x, y')
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
    dir0 = (PolynomialRing(QQbar, x)(x**degree - dir0)).roots()[0][0]
    dir0 /= abs(dir0)

    for k in range(degree):

      rotate = exp(k*2*pi*i/degree)
      direction = rotate * dir0
      print(f'{direction = }')

      x1, y1 = step(H, x, y, r, s, p[x], p[y], direction)
      print(f'{x1, y1 = }')

      steps = walk(H, x, y, r, s, x1, y1)
