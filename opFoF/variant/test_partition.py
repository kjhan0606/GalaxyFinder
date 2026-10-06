#!/usr/bin/env python3
"""Compare opfof_sb_labels with a direct union-find on one small set."""
import ctypes
from pathlib import Path

import numpy as np

SO = Path(__file__).resolve().parent / "libopfof_sb.so"
lib = ctypes.CDLL(str(SO))
F64 = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags="C_CONTIGUOUS")
U64 = np.ctypeslib.ndpointer(dtype=np.uint64, ndim=1, flags="C_CONTIGUOUS")
lib.opfof_sb_labels.argtypes = [ctypes.c_size_t, F64, F64, F64, F64, U64]
lib.opfof_sb_labels.restype = ctypes.c_int


def labels_of(x, y, z, link):
    labels = np.empty(x.size, dtype=np.uint64)
    status = lib.opfof_sb_labels(x.size, x, y, z, link, labels)
    if status:
        raise RuntimeError("opfof_sb_labels status %s" % status)
    return labels


def python_labels(x, y, z, link):
    n = x.size
    parent = np.arange(n)
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(n):
        for j in range(i + 1, n):
            lim = 0.5 * (link[i] + link[j])
            r2 = (x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 + (z[i] - z[j]) ** 2
            if r2 < lim * lim:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[max(ri, rj)] = min(ri, rj)
    return np.array([find(i) for i in range(n)])


def same_partition(a, b):
    _, ia = np.unique(a, return_inverse=True)
    _, ib = np.unique(b, return_inverse=True)
    return np.array_equal(ia, ib)


def main():
    rng = np.random.default_rng(20261006)
    xyz = np.vstack([
        rng.normal([0.0, 0.0, 0.0], 0.02, size=(80, 3)),
        rng.normal([0.4, -0.2, 0.1], 0.05, size=(40, 3)),
        rng.uniform(-1.0, 1.0, size=(30, 3)),
    ])
    link = np.full(xyz.shape[0], 0.03)
    x, y, z = (np.ascontiguousarray(xyz[:, i]) for i in range(3))
    if not same_partition(labels_of(x, y, z, link), python_labels(x, y, z, link)):
        raise SystemExit("partition mismatch")
    # Separation equal to the mean link stays in two components.
    exact = labels_of(np.array([0.0, 1.0]), np.zeros(2), np.zeros(2), np.ones(2))
    if len(np.unique(exact)) != 2:
        raise SystemExit("exact threshold was linked")
    print("PASS")


if __name__ == "__main__":
    main()
