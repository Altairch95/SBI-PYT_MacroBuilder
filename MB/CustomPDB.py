#!/usr/bin/env python
# coding=utf-8
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
import Bio.PDB.NeighborSearch
import sys

class CustomModel(Model):
    """Custom biopython model class that allows models with more than one same id"""
    def add(self, entity):
        """Add a child to the Entity."""
        entity.set_parent(self)
        self.child_list.append(entity)

    def save_to_mmCIF(self, out_name):
        """Saves a model using the given output name in the cwd"""
        io = Bio.PDB.MMCIFIO()
        io.set_structure(self)
        try:
            io.save(out_name+".cif")
            print(out_name+".cif saved")
        except:
            sys.stderr.write("Couldn't save models to current working directory. "
                             "Make sure you have permission to write files")


class CustomChain(Chain):
    """Custom biopython's chain class with more flexibilty. The most important difference is the new attribute
    interactions, a list of interacting residues that is crucial for the macrocomplex building."""
    protein = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}
    dna = {'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T'}
    rna = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U'}

    def __init__(self, chainObject):
        self.child_list = chainObject.child_list
        self._id = chainObject._id
        self.parent = chainObject.parent
        self.interactions = []  # Here there will be stored the chain's known interactions as tuple of residues
        self.xtra = chainObject.xtra
        self.level = chainObject.level
        # The other attibutes are needed for biopython to be able to work.

    def add_interaction_lst(self, lst):
        """Sets the interactions attribute to a list of interaction tuples"""
        self.interactions = lst

    def get_interactions(self, other_chain):
        """Compares the distance between the atoms of two chains and returns a tuple (chain1, chain2) of
        interaction tuples (2,34,50...)"""
        atom_list_1 = list(self.get_atoms())
        atom_list_2 = list(other_chain.get_atoms())
        ns = Bio.PDB.NeighborSearch(atom_list_2)  # Generates a neighbour search tree to speed up distance calculations
        interaction_res_1 = set()
        interaction_res_2 = set()
        for atom in atom_list_1:
            inter_atoms = ns.search(atom.coord, 3.5)  # List of atoms closer than 3.5 A than the given atom
            if len(inter_atoms) > 0:  # If there is interacting atoms:
                interaction_res_1.add(atom.get_parent().id[1])  # Add the corresponding atom's residue to the set
                for iatom in inter_atoms:
                    interaction_res_2.add(iatom.get_parent().id[1]) # Same for the other interacting atoms
        return tuple(sorted(interaction_res_1)), tuple(sorted(interaction_res_2))  # Transform to tuples and return

    def get_sequence(self):
        """Returns the chain's sequence, it can be a protein, DNA or RNA sequence"""
        seq = ""
        flag = "prot"
        first_residue_name = self.child_list[0].resname.strip()  # Get the first residue name to see what kind
        # of sequence it is
        if first_residue_name not in self.protein:  # If residue name isn't in protein name dictionary
            if "D" in first_residue_name:  # If the residue name has letter D, it means it's DNA
                flag = "dna"
            else:
                flag = "rna"  # Otherwise it means it's a RNA sequence
        if flag == "prot":  # If the sequence is a protein, use the protein dictionary
            for res in self:
                if res.id[0] == " ":
                    seq += self.protein[res.resname]
        elif flag == "dna":  # If its DNA, use the corresponding dictionary
            for res in self:
                if res.id[0] == " ":
                    seq += self.dna[res.resname.strip()]
        else:
            for res in self:  # Same for RNA
                if res.id[0] == " ":
                    seq += self.rna[res.resname.strip()]
        return seq

    def get_common_atoms(self, other):
        """Compares the list of atoms of two chains and returns an even tuple of atoms"""
        self_atoms = sorted(self.get_atoms())  # Generates a sorted list of atoms to be able to compare them
        other_atoms = sorted(other.get_atoms())
        len_self = len(self_atoms)
        len_other = len(other_atoms)
        # Return the atom list sliced by the limitant distance
        if len_self > len_other:
            return self_atoms[:len_other], other_atoms
        elif len_other > len_self:
            return self_atoms, other_atoms[:len_self]
        else: # If they are equal, just return the atoms lists
            return self_atoms, other_atoms
