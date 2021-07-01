# Settings

Simulation settings are read from an input file. The command-line flag `-c <path>` specifies the path to a `.cfg` file.
By default, this path is `input/input.cfg` relative to the project root directory.
The fields in `ìnput.cfg` will modify values in the namespace `settings` documented below.

Each value is set with the following syntax:

```
<namespace>::<variable> = <value>   [//Optional comment] 
```
For example, the value of

```c++
namespace settings::precision{
    ...
    double svd_threshold = 1e-12;
    ...
}
```

can be modified from `ìnput.cfg` with the line

```
precision::svd_threshold = 1e-12    //Optional comment
```



## Input

```{eval-rst}
.. doxygennamespace:: settings::input
   :project: DMRG++
   :path: config/nmspc_settings.h
```

## Output

```{eval-rst}
.. doxygennamespace:: settings::output
   :project: DMRG++
   :path: config/nmspc_settings.h
```

## Model Hamiltonian


```{eval-rst}
.. doxygennamespace:: settings::model
   :project: DMRG++
   :path: config/nmspc_settings.h
```

## Algorithms

### iDMRG

```{eval-rst}
.. doxygennamespace:: settings::idmrg
   :project: DMRG++
   :path: config/nmspc_settings.h
```

### fDMRG

```{eval-rst}
.. doxygennamespace:: settings::fdmrg
   :project: DMRG++
   :path: config/nmspc_settings.h
```

### xDMRG

```{eval-rst}
.. doxygennamespace:: settings::xdmrg
   :project: DMRG++
   :path: config/nmspc_settings.h
```


### iTEBD

```{eval-rst}
.. doxygennamespace:: settings::itebd
   :project: DMRG++
   :path: config/nmspc_settings.h
```

### fLBIT

```{eval-rst}
.. doxygennamespace:: settings::flbit
   :project: DMRG++
   :path: config/nmspc_settings.h
```


## Algorithm Strategies

```{eval-rst}
.. doxygennamespace:: settings::strategy
   :project: DMRG++
   :path: config/nmspc_settings.h
```

## Algorithm Precision

```{eval-rst}
.. doxygennamespace:: settings::precision
   :project: DMRG++
   :path: config/nmspc_settings.h
```


## Threading

```{eval-rst}
.. doxygennamespace:: settings::threading
   :project: DMRG++
   :path: config/nmspc_settings.h
```


## Console

```{eval-rst}
.. doxygennamespace:: settings::console
   :project: DMRG++
   :path: config/nmspc_settings.h
```