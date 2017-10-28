#!/bin/bash 
#!/bin/bash

# trap ctrl-c and call ctrl_c()
trap ctrl_c INT
./build/Release/DMRG++
