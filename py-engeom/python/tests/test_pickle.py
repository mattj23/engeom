import pickle
from engeom.geom2 import Vector2, Point2, SurfacePoint2, Circle2, Segment2, Arc2
from engeom.geom3 import Vector3, Point3, SurfacePoint3, Plane3


def test_pickle_vector2():
    v = Vector2(1.0, 2.0)
    pickled_v = pickle.dumps(v)
    unpickled_v = pickle.loads(pickled_v)
    assert v == unpickled_v


def test_pickle_point2():
    p = Point2(1.0, 2.0)
    pickled_p = pickle.dumps(p)
    unpickled_p = pickle.loads(pickled_p)
    assert p == unpickled_p


def test_pickle_surface_point2():
    sp = SurfacePoint2(1.0, 2.0, 0.0, 1.0)
    pickled_sp = pickle.dumps(sp)
    unpickled_sp = pickle.loads(pickled_sp)
    assert sp == unpickled_sp


def test_pickle_circle2():
    c = Circle2(1.0, 2.0, 3.0)
    pickled_c = pickle.dumps(c)
    unpickled_c = pickle.loads(pickled_c)
    assert c == unpickled_c


def test_pickle_segment2():
    s = Segment2(1.0, 2.0, 3.0, 4.0)
    pickled_s = pickle.dumps(s)
    unpickled_s = pickle.loads(pickled_s)
    assert s == unpickled_s


def test_pickle_arc2():
    arc = Arc2(1.0, 2.0, 3.0, 0.5, 1.0)
    pickled_arc = pickle.dumps(arc)
    unpickled_arc = pickle.loads(pickled_arc)
    assert arc == unpickled_arc


def test_pickle_vector3():
    v = Vector3(1.0, 2.0, 3.0)
    pickled_v = pickle.dumps(v)
    unpickled_v = pickle.loads(pickled_v)
    assert v == unpickled_v


def test_pickle_point3():
    p = Point3(1.0, 2.0, 3.0)
    pickled_p = pickle.dumps(p)
    unpickled_p = pickle.loads(pickled_p)
    assert p == unpickled_p


def test_pickle_surface_point3():
    sp = SurfacePoint3(1.0, 2.0, 3.0, 0.0, 0.0, 1.0)
    pickled_sp = pickle.dumps(sp)
    unpickled_sp = pickle.loads(pickled_sp)
    assert sp == unpickled_sp


def test_pickle_plane3():
    plane = Plane3(0.0, 0.0, 1.0, -2.0)
    pickled_plane = pickle.dumps(plane)
    unpickled_plane = pickle.loads(pickled_plane)
    assert plane == unpickled_plane
