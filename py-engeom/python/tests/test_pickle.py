import pickle
from engeom.geom3 import Vector3


def test_pickle_vector3():
    v = Vector3(1.0, 2.0, 3.0)
    pickled_v = pickle.dumps(v)
    unpickled_v = pickle.loads(pickled_v)
    assert v == unpickled_v

