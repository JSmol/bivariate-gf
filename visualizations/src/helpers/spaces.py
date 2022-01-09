import pyglet as pg
from math import ceil

class ComplexPlane():

  def __init__(self, window, batch, scale=200):
    self.__window = window
    self.__scale = scale
    self.__make_axis(batch)

  @property
  def scale(self):
    return self.__scale
  
  def c2p(self, z):
    return self.scale * z.real + self.__window.width/2, self.scale * z.imag + self.__window.height/2

  def p2c(self, x, y):
    return complex(x - self.__window.width/2, y - self.__window.height/2) / self.scale

  def __make_axis(self, batch):

    ticksize = 1/16

    left = self.p2c(0, self.__window.height/2)
    hsx, hsy = self.c2p(left)
    htx, hty = self.c2p(-left)
    self.hline = pg.shapes.Line(hsx, hsy, htx, hty, batch=batch)

    down = self.p2c(self.__window.width/2, 0)
    vsx, vsy = self.c2p(down)
    vtx, vty = self.c2p(-down)
    self.vline = pg.shapes.Line(vsx, vsy, vtx, vty, batch=batch)

    self.ticks = []
    for p in range(ceil(left.real), -ceil(left.real) + 1):
      x1, y1 = self.c2p(complex(p, -ticksize))
      x2, y2 = self.c2p(complex(p, +ticksize))
      self.ticks.append(pg.shapes.Line(x1, y1, x2, y2, batch=batch))

    for p in range(ceil(down.imag), -ceil(down.imag) + 1):
      x1, y1 = self.c2p(complex(-ticksize, p))
      x2, y2 = self.c2p(complex(+ticksize, p))
      self.ticks.append(pg.shapes.Line(x1, y1, x2, y2, batch=batch))

