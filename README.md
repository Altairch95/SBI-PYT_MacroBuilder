# MacrocomplexBuilder
## Constructing macromolecular complexes. 

*Pau Badia, Altaïr C.Hernández and Natàlia Segura*

## **Index**

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->
- [MacrocomplexBuilder](#macrocomplexbuilder)
- [Software Requirements](#software-requirements)
- [Download and Installation](#download-and-installation)
  - [Package tree](#package-tree)
- [Input Files](#input-files)
- [Tutorial](#tutorial)
  - [Command line arguments](#command-line-arguments) 
  - [GUI](#gui)
- [Analysis of examples](#analysis-of-examples)
  - [Enterovirus](#enterovirus)
  - [Proteosome](#proteosome)
  - [Nucleosome](#nucleosome) 
- [Next steps](#next-steps)
- [FAQS](#FAQS)
  
<!-- /TOC -->



### MacrocomplexBuilder

MacrocomplexBuilder is a stand-alone python3 program developed by **Pau Badia i Mompel**, **Altaïr C. Hernández** and **Natàlia Segura Alabart**. It builds protein macrocomplexes taking a set of protein-protein, protein-RNA, protein-DNA, RNA - DNA, RNA - RNA, and/or DNA - DNA interactions. This software could serve to study quaternary structures that are difficult to study *in vivo*.

Below is shown how to install and use this program as a stand-alone command line script (executing the script *MacroB.py*) or with the Graphical User Interface (*MB_GUI.py*).


### Software Requirements

These are the software and its versions required for the MacrocomplexBuilder functionality and execution:

  * [Python 3.6](https://www.python.org/downloads/)
 

For the GUI the following ones are also necessary:

  * [Tkinter (for the GUI interface)](https://wiki.python.org/moin/TkInter)
  * [Pymol](https://pymol.org/2/)


### Download and Installation

You can download our package using Git with the next command:
 
```bash
  $ git clone https://github.com/NataliaSegura/Macrocomplex-builder.git
  $ cd Macrocomplex-Builder
 ```
At this pont, the directory Macrocomplex-Builder should contain the files and directories described bellow:

#### Package tree

The package has the following structure:

    Macrocomplex-Builder/
      README.md
      README.pdf
      setup.py
      MacroBuilder/
          __init__.py
          MacroB.py
          CustomPDB.py
          MB_GUI.py
          MB_launcher.py
      Examples/
          enterovirus/
          hemo/
          microtubul/
          nucleosome/
          phosphate/
          proteasome/
      templates/
          1pma_proteosome.pdb.zip
          3kuy_nucleosome.pdb
          1a3n_hemoglobin.pdb
          3j23_enterovirus.pdb
          5syg_microtubul.pdb

      doc/
          report.md
          CustomPDB.m.html
          MB_GUI.m.html
          MBlauncher.m.html
          MacroB.m.html
          index.html
          

* README.md, README.pdf: the files containing the tutorial and information about our application.
* MB: a fold with the following scripts:
  - MBlauncher.py: the command-line script to launch the program.
  - MacroB.py: a module requiered by MBlauncher.py where are defined the classes of the program.
  - CustomPDB.py: a module required by MacroB.py where are defined the functions of the program.
  - MB_GUI.py: a module to launch the Graphical User Interface.
* Examples: a directory with several examples stored in sub-directories that serve as input to the program.
* Images: all the images used to generate the README and the Report.
* Templates: the raw PDB files from which we extracted the example pairwise interactions.
* setup.py script: to install the program in the python side-packages.
* doc: a folder with the report.md fie and report.pdf, as well as the documentation.

Check that all this information has been correctly downloaded and that there is the script called *setup.py*.

In order to be able to use all the scprits provided in MacrocomplexBuilder the user has to install the package in the python site-packages.

```bash
   $ sudo python3 setup.py install
```
Be sure to have the dependencies previously stated.


### Input Files

This program needs an input of PDB files holding the protein pairwise interactions needed to reconstruct the desired macrocomplex. The program can handle those scenarios: 

* The same sequence appearing in different PDB files has not to be identical, it can handle 95% of identity. 
* The same sequence appearing in different files with the same and/or different names. 
* Repeated chain interactions are not requiered as inputs (i.e. interaction A-A 10 times is treated as a single PDB). It solves infinite structures (i.e. *Microtuble*).
* Pairwise interactions wrongly given to the program. The program threshold for considering two chains as interacting together is 3.5 Amstrongs. If the user gives interactions with bigger distance, they are not considered as such.
* A template PDB file containing the structure of the model to use it as a guideline.


All the models that only have protein - protein interactions can be dispalyed with both USCF CHIMERA and Pymol but when the model has nucleic acid interactions it must be opened with Pymol.


### Tutorial

In this section we make a brief explanation of how to use MacrocomplexBuilder.

#### Command line arguments

```bash
    $ MBlauncher.py -h

    usage: MBlauncher.py [-h] -i DIRECTORY [-o OUTPUT] [-c MAX_CHAINS]
                         [-n NUM_MODELS] [-d] [-v]
                         [-t TEMPLATE | -s STOICH_STRING]

    MacrocomplexBuilder is a python program designed to generate macrocomplex
    structures from simple pair inetractions

    optional arguments:
      -h, --help        show this help message and exit
      -i DIRECTORY      Input directory where the pair interaction pdbs are
                        located
      -o OUTPUT         Name of the output file, no extension is needed
      -c MAX_CHAINS     Maximum number of chains that the user wants in the model
      -n NUM_MODELS     Number of models that the program will compute
      -d, --dirty       Generates an output file for each added chain to track how
                        the program builds the complex
      -v, --verbose     Shows what the program is doing
      -t TEMPLATE       To discriminate against different models, a template can
                        be given to calculate the RMSD
      -s STOICH_STRING  The user can also give the desired stechiometry in this
                        format: A:6,B:11,C:2 ...
```

#### GUI

Another way to use the program is using the the GUI. To do so run the following command:

```bash
$ MB_GUI.py
```
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/whole_GUI.png" alt="whole_GUI" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>GUI structure</b></h4>
          <p>Here is the main structure of the GUI. First there is the Options widged where the user can select what parammeters to use for the modeling. Then there is the console panel where both the stderr and stdout will be shown, if the user wants more information there is the option verbose which will print each action that the program does. Next there is the Sequence widget, where the model sequences and their id's will be shown. At its bottom, there's the Structure composition panel where the model's chain composition is shown. Finally, a n image of the resulting model is shown at the bottom. (For Sequence, Structure and Image, all of the information shown is from the first model). </p>
        </div>
      </div>
    </div>
  </div>

To use the program, first a directory must be selected. One can do this by going to File>Select Directory:
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/selecting_directory.png" alt="select_dir" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Selecting a directory</b></h4>
          <p>To correctly select an input directory, the user must enter INSIDE the directory in the navigation window and press OK. </p>
        </div>
      </div>
    </div>
  </div>

Then, the user can change any of the options in the option panel, let's make just one model of a Proteosome to show it really works and select the option verbose to see what the program is doing. We press RUN to tun the program:
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/options_before_run.png" alt="options_b4_run" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Modify options before run</b></h4>
          <p>Here the model will make just one model, with the name proteasoma_1.cif and without checking structure composition. </p>
        </div>
      </div>
    </div>
  </div>

Finally we can see how the program has build this proteasoma really fast. 
<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/final_result.png" alt="final_result" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Final result</b></h4>
          <p>The sequence, structure composition, structure and console have been updated</p>
        </div>
      </div>
    </div>
  </div>
  
### Analysis of examples

To get a better understanding of how to run the programme properly, we show different examples that represent different inputs that may be provided. The main aspects that may differ the inputs are: number of different chain interactions and number of atoms of the whole macrocomplex.

#### Enterovirus

The 3j23 PDB entry is the Enterovirus 71 empty capsids (https://www.rcsb.org/structure/3j23). EV71 is a single-stranded positive-sense RNA virus and a causative agent of hand, food, and mouth disease. 3j23 is a macrocomplex with three unique protein chains and a stoichiometry of hetero 180-mer-A60B60C60 with 4, 5 and 7 interaction sites in each chain respectively.
Giving a set of protein-protein interactions, MacrocomplexBuilder is able to construct the whole capsid macrocomplex. Comparing the structural composition of the model versus a template, we can observe any differences between both (*Figures 1*).   

To achieve this complex we can run:

```python
MBlauncher.py -i enterovirus/ -o enterovirus -v
```

Where enterovirus/ is the Directory containing all input files, enterovirus is the file where the output will be saved in the current directory, and -v means that the standard error will be printed. 
This is a clear example of one of the strong points of our program: given 8 interaction pairs is able to produce the whole capsid with 180 chains in correct position.


<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/ent.png" alt="enterovirus_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 1</b></h4>
          <p><i>Comparison between the 3j23 PDB entry (left) and the model created by MacrocomplexBuilder(right)</i></p>
        </div>
      </div>
    </div>
  </div>
  
#### Proteosome

The 1pma PDB entry is a proteosome from *Thermoplasma acidophilum* (https://www.rcsb.org/structure/1PMA). A proteosome is a protein macrocomplex which degrade proteins by proteolysis of the peptide bonds. 1pma is a macrocomplex with two unique protein chains and a stoichiometry of hetero 28-mer-A14B14 with 4 interactions in one chain and 7 interactions in the other.

To achieve this complex we can run:

```python
MBlauncher.py -i proteasome/ -o proteasoma -v -c 28
```

Where proteasome/ is the Directory containing all input files, proteasome is the file where the output will be saved in the current directory, -v means that the standard error will be printed, and 28 means that we are limiting the number of chains the model will have.

Even if in these case it is not necessary to limit the number of chains, we limited to 28 to show that the program will correctly construct the model. This will work with all the models this optional argument is given.

<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/proteosome.png" alt="proteosome_superimposed_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 2</b></h4>
          <p><i>In blue we see the 1pma PDB protein and in brown the model created by MacrocomplexBuilder</i></p>
        </div>
      </div>
    </div>
  </div>


#### Nucleosome

As an optional argument, MacrocomplexBuilder can accept global stechiometry. If the user desires, it can be given and the program will create the model according to the given stechiometry.To achieve this complex we can run:

```python
MBlauncher.py -i nucl/ -o nucl -v -s A6:B2
```
Where nucl/ is the Directory containing all input files, nucl is the file where the output will be saved in the current directory, -v means that the standard error will be printed, and A6:B2 means that the global stechiometry will be this one.

The 3kuy PDB entry is the DNA stretching in the nucleosome core of *Escherichia coli* (https://www.rcsb.org/structure/3kuy).The DNA stretching in the nucleosome core can cause dramatic structural distortions, which may influence compaction and factor recognition in chromatin. It has a Stoichiometry of hetero 8-mer-A6B2.
MacrocomplexBuilder is able to create this protein - nucleic acid macrocomplex with 4 protein chains and 2 nucleotic acid chains with a total of 28 different pairwise interactions.


<div class="row">
    <div class="col-md-12">
      <div class="thumbnail">
        <img src="/images/nucl.png" alt="model_nucleosome_image" style="width:500px;height:400px">
        <div class="caption">
          <h4><b>Figure 3</b></h4>
          <p><i>Comparison between the nucleosome created limiting it stechiometry to A:6,B:2 (right) and the one created without stechiometry limitations (left)</i></p>
        </div>
      </div>
    </div>
</div>


## FAQS

**Do I need PyMOL to launch and use the GUI?**
No, it is not necessary, but it is advisable to install [Pymol](https://pymol.org/2/) in order to see the whole macrocomplex/es once they have been done.

**What is the limit of chains that I can model for infinite macrocomplexes?**
In the GUI it's 1000, on the command line program there is no limit but more than 1000 chains could take a lot of time. By default is 100 in the GUI and 300 in the command line but in both programs it chan be changed. 

**Will heteroatoms also be included in the final macrocomplex?**
Yes indeed, if the inputs pdbs had heteroatoms in them they will be present in the final structure. 

**Can I open the resulting macrocomplex in chimera insted of pymol?**
Yes, but because the macrocomplex files are saved in mmCIF format they take 50 times more to open in chimera than normal pdb files so it's not recommended (at least for big structures like capsid). Also, macrocomplexes with DNA and RNA chains cannot be opened with chimera if the mmCIF format is used, so these structures will have to be opened with pymol. 

**How fast is the program?**
It's quite fast. It shows linear growth during the first 300 added chains and then it's exponential. For modeling the enterovirus capsid it takes around 2 - 3 minutes to compute. 

**In the GUI, why do I get only one image if I selected to create more than one model?**
Because the GUI will only use the information of the first model (ended with _1) to fill the GUI widgets.
