import sys, ase.io, io

# this is an openbabel v3.1.1 python library
#sys.path.append('/home/shuhao/softwares/miniconda3/envs/ani_3.6/lib/python3.6/site-packages/openbabel')
from openbabel import openbabel


def GetOBMolAtomIDList(mol):
    # OBMol._atomIds is a protected feature
    # and cannot be read from python interface
    # Have to rebuild it with this function
    ret = []
    for a in openbabel.OBMolAtomIter(mol):
        ret.append(a.GetId())
    
    return ret

def OBfind(atoms, cmatrix=None):
    # Based on openbabel python interface v3.1.1
    # load xyz_io and use ConnectTheDots() to build bond connections
    # then output SMILES, InChI key and atom ID list of each fragments

    if cmatrix is None:
        cmatrix = atoms.get_cell()
        assert cmatrix, "No cell info found in the input atoms!"

    # no need to optimize this by using io as input
    # since all other engines use ase.Atoms as input 
    # if you use io as input, you also need to modify AtomsQueueGenerator
    # which requires reinvent the indices parser
    xyz_io = io.StringIO()
    ase.io.write(xyz_io, atoms, format='xyz')

    ucell = openbabel.OBUnitCell()
    ucell.SetData(openbabel.vector3(*cmatrix[0]), openbabel.vector3(*cmatrix[1]), openbabel.vector3(*cmatrix[2]))
    ucell.SetSpaceGroup('P 1')
    
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz_io.getvalue())
    mol.CloneData(ucell)
    mol.SetPeriodicMol(True)

    openbabel.OBUnitCell.FillUnitCell(ucell,mol)
    mol.ConnectTheDots()
    obConversion.WriteString(mol)

    obConversion2 = openbabel.OBConversion()
    obConversion2.SetInAndOutFormats("smi", "inchi")
    
    ret = []
    for m in mol.Separate():
        s = obConversion.WriteString(m)
        temp = openbabel.OBMol()
        obConversion2.ReadString(temp, s)
        ret.append((s.rstrip('\t\n'), obConversion2.WriteString(temp), GetOBMolAtomIDList(m)))
    
    return ret