#!/usr/bin/env python
# coding=utf-8
from Bio.PDB.PDBParser import PDBParser
from os import listdir
import Bio.PDB.NeighborSearch
from Bio.pairwise2 import align
import random
import sys
from MB.CustomPDB import CustomModel, CustomChain
import os


def read_pdbs(directory, verbose=False):
    """Reads the input directory and generates pdb models"""
    if verbose:
        print("Reading pdb input files from %s" % directory)
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    if os.path.isdir(directory) and directory.endswith("/"):
        try:
            pdbmodels = [parser.get_structure("Model_pair", directory + f)[0] for
                    f in listdir(directory) if f.endswith(".pdb")] #  Generates pdb objects for files that end with .pdb
        except:
            sys.stderr.write("PDB files couldn't be opened. Please, revise that their format is correct.")
            sys.exit(1)
    else:
        sys.stderr.write("Directory %s doesn't exists, please select a valid directory." % directory)
        sys.exit(1)
    if not bool(pdbmodels):  # If no pdb instance is generated
        sys.stderr.write("No pdb files where read. Please make sure the given directory contains pdb files. ")
        sys.exit(1)
    for model in pdbmodels:
        if len(model.child_list) != 2:
            sys.stderr.write("A pdb input file doesn't contains two chains. Please, all input pdbs must only contain "
                             "two chains.")
            sys.exit(1)
    if verbose:
        print("Pdb objects stored")
    return pdbmodels


def get_new_id(iterator):
    """Returns a new id that is not in the given iterator"""
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789ñÑçÇ'
    for l in letters:
        if l not in iterator:
            return l
    sys.stderr.write("Too many different chains given. The program can only handle modeling"
                     " a maximum of 66 different sequences")
    exit(1)


def has_homolgs(target_seq, known_seqs):
    """Checks if a given sequence is an homolog of any of the known sequences and returns it"""
    for k_seq in known_seqs:
        alignment = align.globalxx(target_seq, k_seq)[0] # Generates and selects the alignment with the best score
        aln_seq_1 = alignment[0]  # Get the first sequence of the alignment (with - as gaps )
        aln_seq_2 = alignment[1]  # Get the second one
        al_length = len(alignment[0])
        ident = sum(base1 == base2 for base1, base2 in zip(aln_seq_1, aln_seq_2))  # Calculate number of identities
        if ident/al_length >= 0.95:  # If 95% of identity, return known sequence
            return k_seq


def unify_ids(pdbmodels, verbose=False):
    """Unifies chain identifiers and updates the pdb models to CustomModel class"""
    seq_dict = dict()  # Dictionary where sequences are keys and ids are values
    if verbose:
        print("Unifying Ids")
    for i in range(len(pdbmodels)):
        pdb = pdbmodels[i]
        model = CustomModel(str(i))  # Transforms model to CustomModel instance
        for chain in pdb:
            chain = CustomChain(chain)  # Transforms chain to CustomChain instance
            chain.parent = None  # Removes previous parent from chain
            chain_seq = chain.get_sequence()
            if chain_seq not in seq_dict:
                if not seq_dict:  # If the sequence dictionary is empty
                    new_id = get_new_id(seq_dict.values())  # Get first id (A)
                    seq_dict[chain_seq] = new_id  # Set the first sequence as key and A as value
                    chain.id = new_id  # Also update chain id to A
                else:  # If dictionary is not empty
                    sequences = seq_dict.keys()
                    homolog_seq = has_homolgs(chain_seq, sequences)  # Check if sequence has homology
                    if homolog_seq:
                        homolog_id = seq_dict[homolog_seq]  # Get homolog id from seq_dict
                        seq_dict[chain_seq] = homolog_id  # Set this sequence with the homolog id as value
                        chain.id = homolog_id   # Also change chain id to homolog's
                    else:
                        new_id = get_new_id(seq_dict.values())  # Otherwise generate a new id
                        seq_dict[chain_seq] = new_id
                        chain.id = new_id
            else:
                chain.id = seq_dict[chain_seq]  # If chain is already in seq_dict, update chain object id
            model.add(chain)
        pdbmodels[i] = model  # Update pdbmodels list with the updated model
    if verbose:
        print("Ids unified")
    return seq_dict


def get_interaction_dict(clean_pdbs, verbose=False):
    """Generates interaction dictionary. This is a dictionary of dictionaries. Each chain id is the key of the first
    dictionary, and as value a dictionary. The dictionary inside has tuples of interactiong resiudes (their number)
    from chain 1 to 2 as keys and for value it has a tuple of chain object 1, chain object 2 and interaction resiudes from 2 to 1:
    For instance:
    {Chain1 id : {
        (first interacting residues number from chain 1 to 2):  (Chain1 instance,
                                                                Chain2 instance,
                                                                interacting residues number from chain 2 to 1)
                                                                }
        (second interacting resiudes...) : (... , ... , ... )
        ...
     Chain2 id : {
        (first interacting residues number from chain 2 to 1):  (Chain2 instance,
                                                                (Chain1 instance,
                                                                interacting residues number from chain 1 to 1)
        (second interacting resiudes...) : (... , ... , ... )
        ...
    """
    interaction_dict = dict()
    if verbose:
        print("Generating interaction dictionary...")
    for pdb in clean_pdbs:
        chain1, chain2 = list(pdb.get_chains())
        inter1_2, inter2_1 = chain1.get_interactions(chain2)  # Generates interaction tuples,
        # from chain 1 to 2 and from 2 to 1. For instance:
        # inter1_2 = (2,40,120)
        # inter2_1=(34, 20)
        if inter1_2 != ():  # If tuple is not empty (there is an interaction)
            interaction_dict.setdefault(chain1.id, dict())[inter1_2] = (chain1, chain2, inter2_1) # Update dictionary
        if inter2_1 != ():  # Same for the other interaction
            interaction_dict.setdefault(chain2.id, dict())[inter2_1] = (chain2, chain1, inter1_2)
    if verbose:
        print("Interaction dictionary generated")
    return interaction_dict


def update_interactions_dict(interaction_dict, verbose=False):
    """Updates the interactions attribute of each chain inside the dictionary"""
    if verbose:
        print("Updating interaction dictionary...")
    for chain in interaction_dict:
        for interaction_tple in interaction_dict[chain]:
            chain1, chain2, ref_inter = interaction_dict[chain][interaction_tple] # De-packs the interaction tuple value
            # Generates a list with all interactions of that key minus the current one
            chain1_filtered_interactions_lst = [x for x in interaction_dict[chain1.id].keys() if x != interaction_tple]
            # Updates chain_interaction attribute with the list of interaction tuples
            chain1.add_interaction_lst(chain1_filtered_interactions_lst)
            # Generates a list with all interactions minus the one from chain2 to chain1
            chain2_filtered_interactions_lst = [x for x in interaction_dict[chain2.id].keys() if x != ref_inter]
            chain2.add_interaction_lst(chain2_filtered_interactions_lst)
            parent = chain1.get_parent()
            parent.child_list = [chain1, chain2] # Updates the model chains with these new updated chains
            interaction_dict[chain][interaction_tple] = chain1, chain2 # Updates the dictionary value
            # (now without the interaction from chain 2 to 1)
    if verbose:
        print("Interaction dictionary updated")


def has_clashes(move_atoms, model):
    """Compares the atoms backbone atoms of the moving chain with the backbone atoms of the model"""
    backbone = {"CA", "C1\'"}
    chain_atoms = [atom for atom in move_atoms if atom.id in backbone]  # Gets only the backbone atoms
    model_atoms = [atom for atom in model.get_atoms() if atom.id in backbone]
    ns = Bio.PDB.NeighborSearch(model_atoms)  # Generates a neigbour search tree to speed up distance calculations
    clashes = 0
    for atom in chain_atoms:
        clashes += bool(ns.search(atom.coord, 2))  # If this atom shows clashes, add 1 to the clashes counter
    if clashes/len(chain_atoms) >= 0.03:  # If more than 3% of atoms show clashes return yes
        return True
    else:  # Otherwise return no
        return False


def get_starting_model(interaction_dict, verbose=False):
    """Returns as a starting model one of the CustomModel of the chain with more recorded interactions"""
    if verbose:
        print("Selecting best interaction from where to start modeling...")
    max_len = 0
    interaction_key = None
    for key in interaction_dict:
        length = len(interaction_dict[key])
        if length > max_len:  # If a chain has more interactions
            max_len = length  # Updates maximum interaction umber
            interaction_key = key  # Updates maximum interaction chain key
    inter_tple = next(iter(interaction_dict[interaction_key]))  # Gets an interaction tuple of the chain
    # with more interactions
    if verbose:
        print("First two chains added")
    return interaction_dict[interaction_key][inter_tple][0].parent.copy()  # Returns the model of the interaction tuple


def generate_model_profile(model):
    """Generates a dictionary with the id chain as key and the number of repetitions of this chain as values"""
    profile = {}  # { "A":1, "B":4, ...}
    for chain in model:
        profile.setdefault(chain.id, 0)  # If id not in dic, set it to 0
        profile[chain.id] += 1 # Sum 1 to the counter of the id
    return profile


def save_results(out_models, output):
    """Saves the resulting models into cif files (at the current working directory"""
    i = 1
    print("Saving models...")
    path = os.getcwd()
    for model in out_models:  # Saves all models in the current working directory
        model.save_to_mmCIF(path+"/"+output + "_" + str(i))
        i += 1
    print("Done\n")


def main_loop(num_models, output, interaction_dict, verbose=False, max_chains=100, dirty=False,
              stech_dict=False):
    """Using the interaction dictionary, this function generates macrocomplex model/s. It begins with a template model
    and starts adding chains until conditions allow. Finally it returns a list of pdb instance models."""
    out_models = []
    for i in range(1, num_models + 1):
        print("Macrocomplex " + str(i) + " ...")
        macrocomplex = get_starting_model(interaction_dict, verbose).copy()  # Selects a starting model
        model_stech = generate_model_profile(macrocomplex)  # Generates the stechometry of the first two chains
        macrocomplex.id = "Model_" + str(i)
        run = True  # WHile this variable is true, the program will keep trying to add chains to the macrocomplex
        num_of_chains = 2  # The model starts with 2 chains already
        num_empty_chains = 0  # NUmber of chains that have all their interactions depleted
        while run:
            for chain in macrocomplex:  # Iterates the macrocomplex chains
                if num_of_chains < max_chains:  # If the number of chains still hasn't reached the maximum allowed
                    if chain.interactions:  # If this chain still has pending interactions
                        random.shuffle(chain.interactions)  # Shuffle the interactions list (to avoid
                        # repetitive behaviour)
                        for inter_tple in chain.interactions:
                            if stech_dict:  # If there is stechometry input (either as stirng or template pdb)
                                target_chain_id = interaction_dict[chain.id][inter_tple][1].id  # chain to be added
                                model_stech.setdefault(target_chain_id, 0)
                                model_number_chain = model_stech[target_chain_id]  # Get the number of repetitions
                                stech_dict.setdefault(target_chain_id, 0)
                                if stech_dict[target_chain_id] <= model_number_chain:  # If the number of this target
                                    # chain would surpass the stechemestry given, don't add the chain and
                                    if verbose:
                                        print("(S) Chain NOT added: interaction " + chain.id + ": " +
                                            str(inter_tple[:1]) + " ... " + str(inter_tple[-1]) + " to " + target_chain_id)
                                    continue # jump to the next interaction tuple
                            fix, to_move = interaction_dict[chain.id][inter_tple]  # Get the interaction chain instances
                            sup = Bio.PDB.Superimposer()  # Generates a superimposer instance
                            chain_atoms, fix_atoms = chain.get_common_atoms(fix) # Get common atoms between the
                            # macrocomplex chain and the one in the interaction dictionary
                            sup.set_atoms(chain_atoms, fix_atoms)  # Generate the superposition
                            move = to_move.copy()  # Make a copy of the chain to move
                            sup.apply(move)  # Apply superposition matrix
                            move_atoms = sorted(move.get_atoms())
                            # Now it checks if the target chain has clashes with the model
                            if not has_clashes(move_atoms, macrocomplex):  # If it hasn't
                                if verbose:
                                    print("Chain " + str(num_of_chains) + " added: interaction " + chain.id + ": " +
                                          str(inter_tple[0]) + " ... " + str(inter_tple[-1]) + " to " + move.id)
                                move.parent = None  # Sets the parent to none to evade biopython's strict id policy
                                macrocomplex.add(move)  # Adds the target chain to the model
                                model_stech.setdefault(move.id, 0)  # Updates stech dict
                                model_stech[move.id] += 1
                                num_of_chains += 1
                                if dirty:  # Generates a cif file for each step in the building of the model
                                    macrocomplex.save_to_mmCIF(output + str(i) + "_tmp_" + str(num_of_chains))
                            elif verbose:  # If it has don't add the target chain
                                print("Chain NOT added: interaction " + chain.id + ": " +
                                      str(inter_tple[:1]) + " ... " + str(inter_tple[-1]) + " to " + move.id)
                        chain.interactions = False  # Set the interaction attribute to 0, this chain now will be ignored
                    else:
                        if verbose:
                            print("Chain " + chain.id + " empty")
                        num_empty_chains += 1
                else:
                    run = False  # When the maximum chain treshold is reached stop running
                    break
            if num_empty_chains >= len(macrocomplex):  # If all chains are empty of interactions stop running
                run = False
        if verbose:
            stechometry_string = ""  # Print the model's stechometry
            for key in sorted(model_stech.keys()):
                stechometry_string += key + ":" + str(model_stech[key]) + ","
            stechometry_string = stechometry_string[:-1]
            print("Macrocomplex's"+str(i)+" Stoichiometry is: "+stechometry_string)
        print("Macrocomplex " + str(i) + " finished")
        out_models.append(macrocomplex)  # Add model to the models list
    return out_models

def get_template_stech_dict(template, seq_dict, verbose=False):
    """Generates a stechometry dictionary for a given pdb template"""
    template_stech_dict = {}  # Format: { "A": 2, "B": 3, ...}, where key is chain id and value
    # is the number of repetitions
    parser = PDBParser(PERMISSIVE=1, QUIET=True)
    template_object = parser.get_structure("template", template)[0]  # Generates pdb template object
    for chain in template_object:
        chain = CustomChain(chain)  # Transforms pdb chain object to CustomChain instance
        chain.parent = None  # Removes previous parent to evade biopython errors of id repetitions
        chain_seq = chain.get_sequence()
        if chain_seq in seq_dict:
            chain.id = seq_dict[chain_seq]  # Updates the template chain id to the corresponding by its sequence
            template_stech_dict.setdefault(chain.id, 0)
            template_stech_dict[chain.id] += 1  # Adds to chain id counter
    if verbose:  # Transforms the stech_dict to a string to be printed
        stechometry_string = ""
        for key in sorted(template_stech_dict.keys()):
            stechometry_string += key+":"+str(template_stech_dict[key])+","
        stechometry_string = stechometry_string[:-1]
        print("Template's Stoichiometry is: "+stechometry_string)
    return template_stech_dict

def get_string_stech_dict(stech_string):
    """Given a stechometry string it transsforms it to a dictionary"""
    stech_dict = {}
    try:
        stech_lst = stech_string.split(",")  # Generates a stech list: ["A:3", "B:2", ...]
        for stech in stech_lst:
            chain, number = stech.split(":")
            stech_dict[chain] = int(number)  # Chain id as key and number as value: { "A": 3, "B": 2, ...}
        return stech_dict
    except:
        sys.stderr.write("Stechometry string format is wrong, please follow this format: A:2,B:11,C:4, ...")
        sys.exit(1)


def build_macrocomplex(directory, output, max_chains=300, num_models=1, template=False, dirty=False, verbose=False, stech_string=False):
    """Main function that integrates all the important steps. First it reads the pdb models ans stores them in a list.
    Then it compares all the chains and unifies the chain ids of the pdb list updating them, it also generates a sequence
    key dictionary. Then it checks at each pdb model for chain interactions and stores them in a dictionary. After it
    updates the dictionary with information of itself. Next it generates the model/s using this interactions. Finally it
    saves the model/s in cif format."""
    print("Program is running, please wait...")
    # Reads and stores pdb objects in a list
    in_pdbmodels = read_pdbs(directory, verbose)
    # Unifies all ids by sequence, updates the pdb list with new chain ids and returns a sequence dictionary: {seq: id,}
    seq_dict = unify_ids(in_pdbmodels, verbose)
    # Checks each pdb object for chain interactions and stores it in a dictionary of dictionaries:
    # {
    #   Chain1_id : { residues_tuple_1_to_2 : chain1_object, chain2_object, residues_tuple_2_to_1}
    #   Chain2_id : {residues_tuple_2_to_1 : chain2_object, chain1_object, residues_tuple_1_to_2}
    #   ...
    # }
    interaction_dict = get_interaction_dict(in_pdbmodels, verbose=verbose)
    # Changes interaction_dict chain objects to CustomChain instances and adds the interactions to each instance
    update_interactions_dict(interaction_dict, verbose)
    stech_dict = {}
    # If a template or a string has been given to set Stoichometry, it generates a dictionary of it
    # { "A":5, "B":2, "C":6, .. }
    if template:
        stech_dict = get_template_stech_dict(template, seq_dict, verbose=verbose)
    elif stech_string:
        stech_dict = get_string_stech_dict(stech_string)
    # Starts iterating the interaction pair with more known interactions and generates the model/s
    out_pdbmodels = main_loop(num_models, output, interaction_dict, verbose, max_chains, dirty, stech_dict=stech_dict)
    # Saves the model/s to ciff format
    save_results(out_pdbmodels, output)
