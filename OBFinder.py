import sys, io
import numpy as np
import ase
from ase import Atoms

# this is an openbabel v3.1.1 python library
sys.path.append('/home/shuhao/softwares/miniconda3/envs/ani_3.6/lib/python3.6/site-packages/openbabel')
import openbabel


def GetOBMolAtomIDList(mol):
    # OBMol._atomIds is a protected feature
    # and cannot be read from python interface
    # Have to rebuild it with this function
    ret = []
    for a in openbabel.OBMolAtomIter(mol):
        ret.append(a.GetId())
    
    return ret

###!!! TODO: 
# This function does not require ase.atom as input at all
# instead we can just take sliced xyz file like io.StringIO
def OBfind(atoms):
    # Based on openbabel python interface v3.1.1
    # load xyz file and use ConnectTheDots() to build bond connections
    # then output SMILES and count how many fragments in SMILES
    # May contain different SMILES that actually represent a same species
    cmatrix = atoms.cell.array
    #print(cmatrix)
    ucell = openbabel.OBUnitCell()
    ucell.SetData(openbabel.vector3(*cmatrix[0]), openbabel.vector3(*cmatrix[1]), openbabel.vector3(*cmatrix[2]))
    ucell.SetSpaceGroup('P 1')
    
    s = io.StringIO()
    ase.io.write(s, atoms, format='xyz')
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, s.getvalue())
    mol.CloneData(ucell)
    mol.SetPeriodicMol(True)

    openbabel.OBUnitCell.FillUnitCell(ucell,mol)
    mol.ConnectTheDots()
    obConversion.WriteString(mol)
    
    ret = []
    for m in mol.Separate():
        ret.append((obConversion.WriteString(m), GetOBMolAtomIDList(m)))
    
    return ret