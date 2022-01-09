# tracks a complex number in a complexplane
class ValueTracker():
  def __init__(self, plane, z0=0):
    self.__plane = plane
    self.__z = z0

  def on_mouse_motion(self, x, y, dx, dy):
    self.__z = self.__plane.p2c(x, y)

  @property
  def z(self):
    return self.__z

  @property
  def re(self):
    return self.__z.real

  @property
  def im(self):
    return self.__z.imag

