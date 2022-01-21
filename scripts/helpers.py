from sage.all import *
# TODO: refactor basically everything...

# height function
def h(x, y, r, s):
  return -r*log(abs(x)) -s*log(abs(y))

# finite height function implemented exactly as is on page 85.
def finite_height(H, r, s):
  P = H.polynomial(QQ).newton_polytope()
  vertices = P.vertices()
  vertices = list(map(vector, vertices))
  for j, vj in enumerate(vertices):
    on_line = negative = positive = False
    for k, vk in enumerate(vertices):
      if j == k: continue
      position = (vj - vk).dot_product(vector([-s, r]))
      if position == 0: on_line = True
      if position < 0: negative = True
      if position > 0: positive = True
    # there is a vk such that the line vj -> vk has slope -s/r and
    # the remainder of the polygon is on one side of the line
    if on_line and not (negative and positive):
      return True
  return False

# TODO: rewrite this function
def termination_criteria(H, x, y, r, s):
  z = var('z')
  poly = sum([abs(t[0]) * z**t[1] for t in H(x=z, y=z).coefficients()])
  c0 = poly(z=0)
  poly = poly - c0
  eps = 16
  while poly(z=eps) >= c0: eps /= 2
  return -(r*log(eps) + s*log(eps)), eps

# returns the saddle/critical and non-smooth points
def get_points(H, x, y, r, s):
  saddle_points = Ideal(
    PolynomialRing(QQbar, [x, y]),
    [H, s*x*diff(H, x) - r*y*diff(H, y)]
  ).variety()
  non_smooth_points = Ideal(
    PolynomialRing(QQbar, [x, y]),
    [H, diff(H, x), diff(H, y)]
  ).variety()
  return saddle_points, non_smooth_points

# computes smallest n and d^n/dx^n h(x, y(x)) in terms of x and y such that Dh(p[x], p[y]) != 0
def degeneracy(H, p, x, y, r, s):
  # p is always a critical point here (otherwise degeneracy is 1)
  assert(diff(H, y)(x=p[x], y=p[y]) != 0)
  DyDx = - diff(H, x) / diff(H, y)
  Dh = - r/x - (s/y)*DyDx
  while n := 1:
    n += 1
    Dh = diff(Dh, x) + diff(Dh, y)*DyDx
    if Dh(x=p[x], y=p[y]) != 0:
      return Dh, n
  assert(False)

# returns true if we can parameterize V as (x, y(x)) in a neighborhood N of p (saddle point)
def param_by_x(H, p, x, y):
  if abs(diff(H, y)(x=p[x], y=p[y])) > 0:
    return True
  else:
    # NOTE: if p is critical we are never here, maybe this function is entirely avoidable?
    return False

