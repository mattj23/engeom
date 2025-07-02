# Python and Rust

While `engeom` is fundamentally a Rust language library, there is a set of bindings allowing some of its features to be used in Python versions `>=3.8`.

The Python bindings use `pyo3` and `maturin` to compile a wrapper around `engeom` as an extension module. While it isn't possible to expose them directly, many of the most useful `engeom` Rust `struct`s are wrapped in `pyo3` classes with the same name and wrapped versions of their functionality.  For simplicity, `numpy` is used as an interface layer, with `numpy.ndarray` objects replacing arrays/slices of `engeom` point and vector types.

Because of the wrapping, using `engeom` through the Python bindings will never be as fast or flexible as direct use of the library in a Rust program. However, it will provide a significant speed improvement over native Python code, and provide access to tested algorithms which would be difficult to replicate using `numpy`'s vectorization.

## Python Installation

```bash
pip install engeom
```
