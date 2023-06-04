# PIQM 2023 Advanced Tutorial: Getting Started

![N|Solid](https://ipi-code.org/images/ipi-logo-alpha.png)

![N|Solid](https://static.wixstatic.com/media/32ed9b_173c699a443e4ff18823ccfd12a9b54f~mv2.png/v1/crop/x_38,y_113,w_997,h_351/fill/w_187,h_66,al_c,q_85,usm_0.66_1.00_0.01,enc_auto/PIOM.png)

We will walk you though advanced path-integral quantum mechanics approaches implemented in the [i-PI](https://ipi-code.org/) code, such as 
- "Advanced thermostats and RP contraction" by  Michele Ceriotti
- "Free energy methods" by Venkat Kapil
- "Path integral approximations to real time correlations" by Mariana Rossi
- "RP Instanton Rate Theory" by Yair Litman
- "Bosonic and Fermionic PIMD" by Barak Hirshberg  

To run this tutorial you simply need to clone this [repository](https://github.com/i-pi/piqm2023-tutorial) and [the piqm2023 branch of the i-pi repository](https://github.com/i-pi/i-pi/tree/piqm2023). Find below instructions on setting up the tools, software, and notebooks to run this tutorial on your PC, laptop, or [cloud](https://courseware.epfl.ch/courses/course-v1:EPFL+X+2022/about).

## What do you need?

If you want to run the tutorial on your laptop or PC, you will need a unix-based operating system with a `python 3.6` or higher. If you have a windows OS you can use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) to setup linux. Alternatively, you can log into the MOOC on [Path Integral Methods in Atomistic Modelling]() and then access the cloud based `noto` interface through one of the chapter exercises or this  [link](https://noto-lti-1.epfl.ch/). Note that sometimes you may have to wait a couple of mins after you log in on the MOOC webpage to access the `noto` interface. 

# Setup

Open a terminal on your PC or laptop or on the `noto` interface and follow the instructions below

### Installing the piqm2023 branch of i-PI

Step 1: Clone the `piqm2023` branch of the i-PI code.
```sh
$ git clone -b piqm2023 https://github.com/i-pi/i-pi.git
```

Step 2: Add the i-PI folder to the default path by including the source command in the `.bashrc` file. 
```sh
$ echo "source IPI_PATH/env.sh" >> .bashrc
```
Here `IPI_PATH` is a placeholder for the correct path. Don't forget to `source IPI_PATH/env.sh` from your open terminal to set correctly the paths to the i-PI code.

Step 3: Compile the drivers that implement simple PES meant for examples and exercises.   
```sh
$ cd IPI_PATH/drivers/f90/
$ make
````
You will have access to `i-pi-driver` executable if you have sourced the `env.sh` file correctly!

Step 4: You can also install the following python packages
```sh
$ python -m pip install -U numpy matplotlib ase chemiscope
```
You can check for successful installation by running the following snippet on python.
```py
import numpy as np
import matplotlib.pyplot as plt
import ase, ase.io
import chemiscope
```

You can check for successful installation by running the following snippet on python.
```py
import numpy as np
import matplotlib.pyplot as plt
import ase, ase.io
import chemiscope
```

Step 5: For the free energy tutorial you will also need plumed2 (version 2.5 is compatible with i-PI). You can follow these instructions.

```sh
$ git clone -b v2.5 git@github.com:plumed/plumed2.git
$ cd plumed2
$ ./configure
$ make
```
This can take up to 15 mins on the noto. 

Remember to source the `sourceme.sh` file. It will add the path to plumed2.5 executables and libraties to releavnt environment variables.
```sh
source PLUMED_PATH/sourceme.sh
```

Now to install the python interface you neeed to 
```sh
$ python -m pip install -U plumed==2.5
```

To check, try the following snippet 
```
import plumed
plumed.Plumed()
```

Remember the above will work only if you have sourced the plumed `sourceme.sh` file.

### Cloning the PIQM 2023 tutorial. 

Cloning the PIQM 2023 tutorial is very easy. 

```sh
git clone https://github.com/i-pi/piqm2023-tutorial
```
