#!/usr/bin/env python
# coding=utf-8
from MB.MacroB import build_macrocomplex
from argparse import ArgumentParser

parser = ArgumentParser(description='MacrocomplexBuilder is a python program designed to generate macrocomplex structures from simple pair inetractions pdb files.')

stech = parser.add_mutually_exclusive_group()
parser.add_argument("-i",
                    dest="directory",
                    action="store",
                    type=str,
                    help="Input directory where the pair interaction pdbs are located",
                    required=True)

parser.add_argument("-o",
                    dest="output",
                    action="store",
                    type=str,
                    default='macrocomplex.pdb',
                    help="Name of the output file, no extension is needed. The output will be saved on the current "
                         "working directory")

parser.add_argument('-c',
                    dest='max_chains',
                    action="store",
                    type=int,
                    default=300,
                    help="Maximum number of chains desired in the model")

parser.add_argument('-n',
                    dest='num_models',
                    action="store",
                    type=int,
                    default=1,
                    help="Number of models that the program will compute")

parser.add_argument('-d',
                    '--dirty',
                    dest='dirty',
                    action="store_true",
                    default=False,
                    help="Generates an output file for each added chain to track how the program builds the complex")

parser.add_argument('-v',
                    '--verbose',
                    dest='verbose',
                    action="store_true",
                    default=False,
                    help="Shows what the program is doing")

stech.add_argument('-t',
                    dest='template',
                    action="store",
                    type=str,
                    default=False,
                    help="Stechiometry given as a template pdb file")

stech.add_argument('-s',
                    dest='stech_string',
                    action='store',
                    type=str,
                    default=False,
                    help="Stechiometry given as a string: A:6,B:11,C:2 ...")

options = parser.parse_args()
build_macrocomplex(options.directory, options.output, options.max_chains, options.num_models, dirty=options.dirty,
                       verbose=options.verbose, template=options.template, stech_string=options.stech_string)
