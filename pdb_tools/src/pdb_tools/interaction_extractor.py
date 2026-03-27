from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB import Structure, Atom
from pathlib import Path

# Create selector class
class InteractionSelector(Select):
    """
    Child class to the PDBIO Select class, which redefines the accept_residue class
    method to select only those chain residues that meet inter-atom distance criteria
    between ligand atoms and molecular chain residue atoms.

    """
    def __init__(self, target_chain_id, ligand_resname, neighbor_residues):
        self.target_chain_id = target_chain_id
        self.ligand_resname = ligand_resname
        self.neighbor_residues = neighbor_residues

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        is_target_chain = chain_id == self.target_chain_id
        is_ligand = residue.get_resname() == self.ligand_resname and is_target_chain
        is_neighbor = residue in self.neighbor_residues and is_target_chain
        return is_ligand or is_neighbor


def find_ligand_atoms(structure: Structure,
                      target_chain: str,
                      ligand_resname: str):
    """
    Find and collect all ligand atoms and return them as a list.

    :param structure: BioPython PDB object containining a nested data structure for
                      biomolecule model(s), molecular chain(s), residues, and atoms.
    :type structure: Structure
    :param target_chain: name (ID) of the target molecular chain.
    :type target_chain: str
    :param ligand_resname: ligand name in PDB file.
    :type ligand_resname: str

    :return ligand_atoms: list of all atoms that comprise the ligand.
    :rtype: list[Atom]
    """
    # COMPLETE THIS FUNCTION!!!
    ligand_atoms = []

    model = structure[0]

    chain = model[target_chain]

    for residue in chain:
        if residue.get_resname() == ligand_resname:
            for atom in residue:
                ligand_atoms.append(atom)

    return ligand_atoms


def find_neighbor_residues(structure: Structure,
                           ligand_atoms: list[Atom],
                           target_chain: str,
                           ligand_resname: str,
                           distance_cutoff: float = 5.0):
    """
    Find all the residues that are within the distance cutoff of ANY ligand atoms
    (and therefore considered an interaction).

    :param structure: BioPython PDB object containining a nested data structure for
                      biomolecule model(s), molecular chain(s), residues, and atoms.
    :type structure: Structure
    :param ligand_atoms: list of all atoms that comprise the target ligand.
    :type ligand_atoms: list
    :param target_chain: name (ID) of the target molecular chain.
    :type target_chain: str
    :param ligand_resname: name (ID) of the target ligand.
    :type ligand_resname: list[Residue]
    :param distance_cutoff: radius (A) defining the distance between molecular chain residue atom
                             & ligand atom that qualifies for an atom-atom interaction.
    :type distance_cutoff: float

    :return: list of amino acid residues within the cutoff radius of
                              any ligand atom (therefore neighbor).
    :rtype: str
    """
    # COMPLETE THIS FUNCTION!!!
    neighbor_residues = []

    model = structure[0]

    chain = model[target_chain]

    for residue in chain:
        if residue.get_resname() == ligand_resname:
            continue
        for amino_acid_atom in residue:
            for ligand_atom in ligand_atoms:
                if amino_acid_atom - ligand_atom <= distance_cutoff:
                    neighbor_residues.append(residue)

    return list(set(neighbor_residues))


def get_protein_ligand_interaction_pdb(pdb_file: str,
                                       chain: str,
                                       ligand_resname: str,
                                       output_file: str = None,
                                       distance_cutoff: float = 5.0):
    """
    Parse (from PDB file) amino acids neighboring ligand within specified distance cutoff.

    :param pdb_file: Protein Data Bank file with structural coordinates for amino acids and ligand.
    :type pdb_file: str
    :param chain: protein chain ID in PDB file.
    :type chain: str
    :param ligand_resname: ligand name in PDB file.
    :type ligand_resname: str
    :param output_file: User-specified output filename for file containing PDB data corresponding
                       specifically to amino acids interacting with the specified ligand within the
                       specified distance cutoff; default "None" returns output file named as "{pdb_file}_out.pdb".
    :type output_file: str or None
    :param distance_cutoff: radius (Angstrom) from ligand coordinates within which amino acids are considered
                           "neighbors" and parsed from the pdb file for return in the output file.
    :type distance_cutoff: float
                           
    :return: path to PDB file filtered for amino acids that are structurally within the cutoff radius of 
             the ligand coordinates.
    :rtype: str
    """
    parser = PDBParser(QUIET=True)

    # ADD A CHECK TO MAKE SURE .PDB FILES ARE USED
    assert Path(pdb_file).suffix.lower() == ".pdb", "Please provide .PDB file."

    structure = parser.get_structure("complex", pdb_file)
    ligand_atoms = find_ligand_atoms(structure, chain, ligand_resname)
    neighbor_residues = find_neighbor_residues(structure, ligand_atoms, chain,
                                               ligand_resname, distance_cutoff)
    selector = InteractionSelector(target_chain_id=chain,
                                   ligand_resname=ligand_resname,
                                   neighbor_residues=neighbor_residues)
    # Write output
    io = PDBIO()
    io.set_structure(structure)
    # ADD NEW ARGUMENT TO SUPPORT CUSTOM OUTPUT NAMES
    if output_file:
        pdb_output = output_file
    else:
        pdb_output = Path(pdb_file).stem +"_out.pdb"
    io.save(pdb_output, selector)
