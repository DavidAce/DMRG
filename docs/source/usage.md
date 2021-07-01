# Usage
The input configuration file `input/input.config` controls the properties of the simulation. `DMRG++` admits custom input files from command-line arguments, e.g. `./build/Release/DMRG++ -i ../CustomFolder/custom.config`.

If no configuration file is given as argument, the default is to look for a configuration file located in `input/input.config` relative to the project root folder.

## Configuration file
The default configuration file in `input/input.config` contains run-time instructions for DMRG++. Here you can choose the type of simulation, type of model, model parameters,
system size, precision as well as settings for profiling, data storage and console verbosity. Read the comments in the file to learn more.

## Models
Two models of quantum systems are included, the *(Random) Transverse-field Ising model* and the *Self-dual (random) transverse-field Ising model* located under `source/model`. To add another model, one currently
needs to hand-craft the *Matrix Product Operator* (MPO) of the model, and add it to the source code via a derived class just as the preexisting ones.
Once implemented, the model is selected using the input configuration file in `input/input.config`, with the option `model::model_type `.

## Output file
After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration file `input/input.config`.
By default this should be in `output/output.h5`. This file will contain values like the final energies, entanglement entropies, entanglement spectrums and
optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.


