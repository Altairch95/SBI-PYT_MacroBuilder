# Setup.py python distribution script to install Macrocomplex-Builder in the computer

from distutils.core import setup

setup(
   name='MB',
   version='0.0.1',
   description='MacroBuilder is a python program designed to 	generate macrocomplex structures from simple pair interactions',
   long_description=open('README.md').read(),
   author='B. Pau, C. Altaïr, S. Natàlia',
   url='https://github.com/Altairch95/4SMacroBuilder',
   packages=['MB'],  
   install_requires=['biopython >= 1.73.0',
		                 'numpy >= 1.16.1'], 
   license='LICENSE.txt',
   classifiers=[
                "Programming Language :: Python :: 3",
	              "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"],
    scripts=["MB/MBlauncher.py" , "MB/MB_GUI.py"]
)
