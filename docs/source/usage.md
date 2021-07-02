# Usage

## Input Configuration File
The default input configuration file `input/input.cfg` sets simulation properties, such as algorithm, system size and precision.
`DMRG++` admits custom input files from command-line arguments, e.g. `./DMRG++ -c path/to/file.cfg`.  
The full list of configurable settings can be found under [Settings](Settings). 

## Output Data File
After execution the results are stored a binary file in HDF5 format. Its location is specified in the configuration file `input/input.cfg`.
By default this should be in `output/output.h5`. This file will contain values like the final energies, entanglement entropies, entanglement spectrums and
optionally the final state in MPS form.

To view the data you can use any hdf5-viewer, such as HDFCompass.

## Model Hamiltonians
These model Hamiltonians of 1D quantum systems are included:

- `ModelType::ising_sdual`: The Self-dual transverse-field Ising model.
- `ModelType::ising_tf_rf`: The Transverse-field Ising model with random on-site field. 
- `ModelType::lbit`: The l-bit Hamiltonian, used to describe a many-body localized phase (MBL) 
   in terms of its local integrals of motion (the l-bits).
 

The Hamiltonians are implemented as *Matrix Product Operators* (MPO), located under `source/tensors/model`.
The model type is selected using the input configuration file in `input/input.cfg`, with the option `settings::model::model_type`.
To add another model, one currently has to implement a new MPO and derive from `class_mpo_site` just like the existing models. 




