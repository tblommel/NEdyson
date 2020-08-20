import numpy as np
from pyscf import gto, scf
import h5py

def int_symmetrize(arr):
    arr = np.array(arr)
    if arr.shape != (arr.shape[0],) * 2:
        raise ValueError("Matrix to symmetrize must be square")

    i, j = arr.nonzero()
    arr[j, i] = arr[i, j]
    return arr


def umatrix_symmetrize(umatrix):
    umatrix = np.array(umatrix)
    if umatrix.shape != (umatrix.shape[0],) * 4:
        raise ValueError("Matrix is not a U-matrix")
    # Generator 1
    i, k, j, l = umatrix.nonzero()
    umatrix[k, i, j, l] = umatrix[i, k, j, l]
    # Generator 2
    i, k, j, l = umatrix.nonzero()
    umatrix[i, k, l, j] = umatrix[i, k, j, l]
    # Generator 3
    i, k, j, l = umatrix.nonzero()
    umatrix[j, l, i, k] = umatrix[i, k, j, l]
    return umatrix


def save(path):
    with h5py.File(path, "w") as h5f:
        h5f["restricted"] = np.int8(restricted)
        h5f["transformed"] = np.int8(transformed)
        h5f["nelectron"] = nelectron
        h5f["spin"] = spin
        h5f["ovlp"] = ovlp
        h5f["hcore"] = hcore
        h5f["uchem"] = uchem
        h5f["enuc"] = enuc
        h5f["fock"] = fock

if __name__ == '__main__':
    mol = gto.M(
        atom = [['Ne', (0, 0, 0)]],
        basis = 'sto-6g',
        unit = 'B',
        verbose = 1,
        )
    mol.build()
    '''
    if you need dimer, use something like
    mol = gto.M(
        atom = [['Ne', (0, 0, 0)], ['Ne', (r0, 0, 0)]],
        basis = 'sto-6g',
        unit = 'B',
        verbose = 1,
        )
    mol.build()
    '''

    restricted = 1
    transformed = 1

    mf = scf.RHF(mol)
    mf.run()
    if mf.converged:
        print("HF converged, energy:", mf.e_tot)

    nelectron = mol.nelectron
    spin = mol.spin
    ovlp = int_symmetrize(mf.get_ovlp())
    hcore = int_symmetrize(mf.get_hcore())
    enuc = mf.energy_nuc()
    fock = int_symmetrize(mf.get_fock())
    uchem = umatrix_symmetrize(mol.intor('int2e_sph'))

    # transform
    rtol = 1e-7
    s_v, s_b = np.linalg.eigh(ovlp)

    # From Markus
    atol = rtol * s_v[-1]
    istart = s_v.searchsorted(atol)
    s_sqrtv = np.sqrt(s_v[istart:])

    x = s_b[:, istart:] * s_sqrtv
    xinv = np.linalg.pinv(x)

    nold, nnew = x.shape

    ovlp = np.eye(nnew)
    # hcore = int_symmetrize(xinv @ hcore @ xinv.conj().T)
    hcore = int_symmetrize(np.matmul(np.matmul(xinv, hcore), xinv.conj().T))
    if restricted:
        # fock = int_symmetrize(xinv @ fock @ xinv.conj().T)
        fock = int_symmetrize(np.matmul(np.matmul(xinv, fock), xinv.conj().T))

    u = np.einsum('iI,IKJL->iKJL', xinv, uchem)
    u = np.einsum('jJ,iKJL->iKjL', xinv, u)
    u = np.einsum('Kk,iKjL->ikjL', xinv.conj().T, u)
    u = np.einsum('Ll,ikjL->ikjl', xinv.conj().T, u)
    uchem = umatrix_symmetrize(u)

    save("./Nesto6g.h5")

