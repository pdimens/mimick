---
label: Installing Mimick
icon: desktop-download
order: 99
---

There are three ways you can install Mimick: conda, pixi, pip. Make note of the names of the code blocks, as they
describe installation and usage variations.

+++ conda
## Install with conda
To install Mimick using conda (or miniconda, mamba, etc.) into an existing environment, it's as simple as:
```you're already in the environment
$ conda install -c conda-forge bioconda::mimick
```
```you're not in the environment
conda install -n env_name -c conda-forge bioconda::mimic
```
where `env_name` is the name of the existing conda environment
you want to install Mimick into. Alternatively, you can create
a new conda environment and install Mimick into it using one command:
```global environment
conda create -n mimick -c bioconda -c conda-forge mimick
```
```local environment
conda create -p path/to/workdir/mimick -c bioconda -c conda-forge mimick
```

Then, run the install script provided by Mimick itself, which will install the internally-bundled MimickLinkedReads.jl
Julia package.
```bash
# activate environment with mimick...
install_mimick
```

### Using with conda
You'll need to activate the environment you installed Mimick into, then call it with:
```run Mimick
mimick options... args...
```

### Updating with conda
If you wish to update Mimick, just replace word `create`/`install` in the commands above with `update`, e.g.:
```you're already in the environment
conda update -c conda-forge bioconda::mimick
```
```you're not in the environment
conda update -n env_name -c conda-forge bioconda::mimick
```

Because of the Julia package backend, you will need to reinstall the Julia package:
```bash
# activate environment with mimick...
install_mimick
```

+++ pixi
### Install with pixi
!!!warning
Julia does not install correctly when using Pixi. **Do not use this method.**

See [here](https://discourse.julialang.org/t/managing-julia-versions-using-pixi/116165/12) for more information
!!!

+++ pip
### Install with pip
If neither conda nor pixi appeal to you, Mimick can be installed with pip too. To do that, you will need to first
download [the latest release](https://github.com/pdimens/mimick/releases) (recommended) or clone
[the git repository](https://github.com/pdimens/mimick):
```download the latest release
# (replace x.x.x with the actual version #)
wget -O mimick.tar.gz https://github.com/pdimens/mimick/releases/download/x.x.x/mimick.x.x.x.tar.gz
tar -xvzf mimick.tar.gz
```
```clone the repository
git clone https://github.com/pdimens/mimick.git
```

Then, you will need to enter the Mimick directory and do a local `pip` installation:
```local pip installation
cd mimick
pip install .
```

Then, run the install script provided by Mimick itself, which will install the internally-bundled MimickLinkedReads.jl
Julia package.
```bash
# activate environment with mimick...
install_mimick
```

### Using with pip
Just call `mimick` from the command line
```call Mimick
mimick options... args...
```

### Updating with pip
You'll need to repeat the [Install with pip](#install-with-pip) process above.

+++
