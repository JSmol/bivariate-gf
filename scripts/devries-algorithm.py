from sage.all import *
from helpers import *

from stepsize import step

LOGGING_FREQ = 100
PRINT_STEPS = False
def path(H, p, j, x, y, r, s, h_stop, radius, critical_points, meta):
  ITERS = 0
  steps = []
  x0, y0 = p[x], p[y]
  near_point = high = False
  # TODO: rewrite this loop
  while not high and not near_point:
    
    if ITERS > 100000:
      print('--- too many iters ---')
      return True, steps

    for k, q in enumerate(critical_points):
      if k == j: break # do not process above j'th point

      # TODO: rewrite so R, eps are not used in this calculation
      R = meta[k]['R']
      Rd = meta[k]['eps']
      if abs(p[x] - x0) < R: # close in independant variable
        if abs(q[y] - y0) < Rd: # close in dependant variable
          assert(not near_point) # can't be close to two saddle points?
          near_point = True
          xpole = meta[j]['xpole']

    if not near_point:
      height = -r*log(abs(x0)) -s*log(abs(y0))
      if height > h_stop:
        high = True
        if abs(x0) < radius:
          xpole = True
        else:
          assert(abs(y0) < radius)
          xpole = False

    if not high and not near_point:
      by_x = True # TODO: work out the param by x thing, param_by_x(H, {x:x0, y:y0}, x, y)
      if by_x:
        Dh = -r/x + (s/y)*(H.diff(x)/H.diff(y))
        dir0 = conjugate(Dh(x=x0, y=y0))
        dir0 /= abs(dir0)
        x0, y0 = step(H, x, y, x0, y0, dir0)

      else: # TODO: figure this out evenually (param by x not possible at x0, y0)
        assert(False) # none of my examples enter this, but I am curious to see one that does...

    ITERS += 1
    if PRINT_STEPS and ITERS % LOGGING_FREQ == 0:
      print('step direction:', dir0)
      print('x0, y0:', x0.n(), y0.n())
      steps.append((real(x0), imag(x0), h(x0, y0, r, s)))

  print('path step iterations:', ITERS)
  return xpole, steps

def MAIN(G, H, x, y, r, s):
  print(f'Calculating F = {G / H} in direction {(r, s)}.')

  # WARNING: edge cases return strings (REWRITE)
  if r == 0 or s == 0: return "Univariate problem."
  if finite_height(H, x, y, r, s): return "There are points at infinity of finite height."

  # sort the points by the height function in decreasing order
  C, NS = get_points(H, x, y, r, s)
  C = sorted(C, reverse=True, key=lambda point: h(point[x], point[y], r, s))
  NS = sorted(NS, reverse=True, key=lambda point: h(point[x], point[y], r, s))

  print('critical points:', *C, sep='\n')
  print('non-smooth points:', *NS, sep='\n')

  heights = [h(point[x], point[y], r, s) for point in C]
  print('heights:', *heights, sep='\n')

  # The height of the highest non-smooth point
  c0 = h(NS[0][x], NS[0][y], r, s) if len(NS) > 0 else -infinity
  print(f'c0: {c0}')

  # TODO: rewrite for better paths storage (use kwarg)
  # store critical point meta data and paths for plotting
  meta = [dict()] * len(C)
  paths = [list()] * len(C)

  h_stop, radius = termination_criteria(H, x, y, r, s)
  print(f'termination_criteria: {h_stop, radius}')

  print('\n\n')
  print(f'--- processing points ---')
  contributing_points = []
  for j, p in enumerate(C):

    # c0 is updated to the highest contributing point, lower points do not contribute to main term asymptotic
    if not h(p[x], p[y], r, s) >= c0: break

    print()
    print('considering point:', p)

    # TODO: rewrite this section (very very carefully)

    # NOTE: for critical points, param_by_x is always true (in bivariate case)
    # meta[j]['by_x'] = param_by_x(H, p, x, y)
    # assert(param_by_x(H, p, x, y))

    # meta[j]['delta'], meta[j]['eps'] = param_nbd(H, p, x, y)
    # meta[j]['R'] = ascent_nbd(H, p, x, y, r, s, by_x=True, **meta[j])
    meta[j]['n'], meta[j]['Dh'] = degeneracy(H, p, x, y, r, s)
    degree = meta[j]['degree']
    Dh = meta[j]['Dh']
    print('saddle point degeneracy:', degree)

    dir0 = conjugate(Dh(x=p[x], y=p[y]))
    # NOTE: taking n'th root doesn't return an algebraic number in sage?
    dir0 = (PolynomialRing(BASEFIELD, x)(x**degree - dir0)).roots()[0][0]
    dir0 /= abs(dir0)

    # For each ascent region: follow ascending path through that region
    # meta[j][path_to_x][k] = does path from point j in direction k go to x = 0?
    meta[j]['path_to_x'] = []
    for k in range(n):

      rotate = exp(k*2*pi*i/degree)
      direction = rotate * dir0
      print(f'{direction = }')

      # the first step is done here since the degeneracy may be greater than 1
      x1, y1 = step(H, x, y, p[x], p[y], direction)
      print(f'{x1, y1 = }')

      to_x, steps = path(H, {x:x1, y:y1}, j, x, y, r, s, h_stop, radius, C, meta)
      print(f'{to_x = }')

      meta[j]['path_to_x'].append(to_x)
      paths[j].append(steps)

    contributer = False
    for k in range(n):
      if meta[j]['path_to_x'][k] != meta[j]['path_to_x'][(k+1) % n]:
        contributer = True
        c0 = heights[j]

    if contributer: contributing_points.append(p)
    else: # all paths go to either x = 0 or y = 0
      meta[j]['xpole'] = meta[j]['path_to_x'][0]

    # there is too much output! (Dh is sometimes massively complicated)
    # print('point meta data:')
    # print(*meta[j].items(), sep='\n')

  # --------------------------------------------------------------------------------
  if False: # NOTE: plotting code (super trees example)

    # hard code the hieght function parameterized by y=g(x) so that we can plot the surface in 3 dim
    g1 = (x**3 - x**2 - sqrt(x**4 - 4*x**3 + 4*x**2 + 4) - 2)/x**5
    g2 = (x**3 - x**2 + sqrt(x**4 - 4*x**3 + 4*x**2 + 4) - 2)/x**5

    def h1(a, b):
      return max(- log(abs(a+b*1j)) - log(abs(g1(x=a+b*1j))), -3)
    def h2(a, b):
      return max(- log(abs(a+b*1j)) - log(abs(g2(x=a+b*1j))), -3)

    # return - log(abs(a+b*i)) - log(abs(sqrt(-1-(a+b*i)**2))) # param by x for 1+x^2+y^2
    # return - log(abs(a+b*i)) - log(abs(1-a-b*i)) # param for 1-x-y

    redim = (-4, 4)
    imdim = (-4, 4)
    surface_reso = 200
    p = plot3d(h1, redim, imdim, plot_points=surface_reso, color='purple', opacity=0.5)
    p += plot3d(h2, redim, imdim, plot_points=surface_reso, color='green', opacity=0.5)

    # where is hstop?
    p += plot3d(lambda a, b: h_stop, redim, imdim, plot_points=surface_reso, color='pink', opacity=0.5)
    
    colors = ['red', 'blue', 'orange']
    for j in range(len(saddle_points)):
      for k in range(len(paths[j])):
        for step in paths[j][k]:
          if redim[0] < step[0] and step[0] < redim[1] and imdim[0] < step[1] and step[1] < imdim[1]:
            p += point(step, size=20, color=colors[j%len(colors)])
        # p += sum([point(step, size=8, color=colors[j%len(colors)]) for step in paths[j][k]])

    p.show()

  return contributing_points

# TESTING:
if __name__ == '__main__':

  def test_case(G, H, x, y):
    print('\n\n\n')
    print('----- test output -----')
    print()
    print(f'Results: {MAIN(G, H, x, y, 1, 1)}')

  R = PolynomialRing(BASEFIELD, 'x, y')
  x, y = R.gens()

  test_case(1, 1 - x - y, x, y)
  test_case(1, 1 + x**2 + y**2, x, y)
  test_case(1, expand((1-x-y) * (1+2*x)), x, y) # non-smooth point example

  # # NOTE: these test cases are more complicated (and slower to run)
  # # critical point at infinity (but a saddle integral is possible)
  # test_case(1, 1 - x - y - x**2*y, x, y)
  # # critical point at infinity (there is no saddle integral possible)
  # test_case(1, 1 - y - x*y, x, y)
  #
  # test_case( # super trees (thesis case study)
  #   2*x**2*y * (2*x**5*y**2 - 3*x**3*y + x + 2*x**2*y - 1),
  #   x**5*y**2 + 2*x**2*y - 2*x**3*y + 4*y + x - 2,
  #   x, y
  # )
  #
  # test_case( # non-smooth super tree version
  #   2*x**2*y * (2*x**5*y**2 - 3*x**3*y + x + 2*x**2*y - 1),
  #   (x**5*y**2 + 2*x**2*y - 2*x**3*y + 4*y + x - 2) * (1 - x),
  #   x, y
  # )

