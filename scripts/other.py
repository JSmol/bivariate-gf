# WARNING: these functions are no longer being used in the main algorithm, so they may no longer work as advertised.
from sage.all import *
from helpers import *

# returns eps such that the following is true:
# H(p[x], y) = 0 and |p[y] - y| < eps implies y = p[y]
def isolate_root(H, p, x, y):
  Hy = H.diff(y)

  # should be true (smooth point p with Hy(p) non zero)
  assert(abs(Hy(x=p[x], y=p[y])) > 0)

  E = var('E')
  lhs = abs(H(x=p[x], y=p[y]+E) - E*Hy(x=p[x], y=p[y])) / abs(E)
  rhs = abs(Hy(x=p[x], y=p[y]))

  eps = 16
  while lhs(E=eps) >= rhs:
    eps /= 2
  return eps

# returns delta and eps such that V is parameterizable by x on the set:
# V intersect (B_delta(p[x]) cross B_eps(p[y]))
def param_nbd(H, p, x, y):
  assert(abs(H(x=p[x], y=p[y])) < ZERO)

  isolate = isolate_root(H, p, x, y)
  delta = abs(p[x]) / 2
  eps = min(isolate, abs(p[y]) / 2)

  Hy = H.diff(y)
  rhs = abs(Hy(x=p[x], y=p[y]))
  assert(rhs > 0)

  # lhs of equations 4.10.1, 4.10.2 and 4.10.3
  D, E = var('D, E')
  P1 = abs(Hy(x=p[x]+D, y=p[y]) - Hy(x=p[x], y=p[y]))
  P2 = abs(H(x=p[x]+D, y=p[y]+E) - H(x=p[x]+D, y=p[y]) - E*Hy(x=p[x]+D, y=p[y])) / E
  P3 = abs(H(x=p[x]+D, y=p[y]) - H(x=p[x], y=p[y]))

  while P1(D=delta) >= rhs / 2:
    delta /= 2

  while P2(D=delta, E=2*eps) >= rhs / 4:
    eps /= 2

  while P3(D=delta) >= (eps / 4) * rhs:
    delta /= 2

  # dealing with approximates.
  assert(abs(P1(D=0)) < ZERO)
  assert(abs((E * P2)(D=0, E=0)) < ZERO)
  assert(abs(P3(D=0, E=0)) < ZERO)

  # inequalities of Lemma 4.10.1
  assert(P1(D=delta) < rhs / 2)
  assert(P2(D=delta, E=2*eps) < rhs / 4)
  assert(P3(D=delta) < (eps / 4) * rhs)

  # doesn't contain (0, 0)
  assert(abs(p[x]) > delta and abs(p[y]) > eps)

  return delta, eps

# returns the radius of a neighborhood of p on which parameterization is well-understood geometrically
def ascent_nbd(H, p, x, y, r, s, delta, eps, n, Dh, by_x, **meta):
  x1 = abs(p[x]) - delta
  y1 = abs(p[y]) - eps

  log_bound = r*log(x1) + s*log(y1) if by_x else r*log(y1) + s*log(x1)
  M = sqrt(log_bound**2 + (2*pi*(r+s))**2)

  P = var('P')
  lhs = ((n + 1 - n*P)*P) / (1 - P)**2
  rhs = delta**n * abs(Dh(x=p[x], y=p[y])) / (factorial(n-1) * sqrt(2) * M)

  rho = 1/2
  while lhs(P=rho) >= rhs:
    rho /= 2
  return rho * delta
  
