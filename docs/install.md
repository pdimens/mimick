# Installing Mimick

There are three ways you can install Mimick: conda, pixi, pip. Make note of the names of the code blocks, as they
describe installation and usage variations.

!> WIP: mimick is not yet available on conda.

<!-- tabs:start -->
#### **conda**

## Install with conda
To install Mimick using conda (or miniconda, mamba, etc.) into an existing environment, it's as simple as:
```you're already in the environment
$ conda install -c conda-forge bioconda::mimick
```
```you're not in the environment
conda install -n env_name bioconda::mimic
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


####  **pixi**

### Install with pixi
Pixi is a new ultra-fast Rust-based environment manager much like conda is. If you haven't tried it out, you ought to give it a shot.
With it, you can install Mimick to be accessible in your PATH, i.e. a "global" installation:
```global install
pixi global install -c conda-forge -c bioconda mimick
```
```local install
pixi init -c conda-forge -c bioconda projectname && cd projectname
pixi add mimick
```

!> Make sure `~/.pixi/bin` is in your PATH using `export PATH=~/.pixi/bin:$PATH`

### Using with pixi
If installed globally, this is very similar to the conda installation:
```install globally
mimick options... args...
```

If you installed Mimick locally using Pixi, you can do one of two things:
1. activate the pixi shell in the directory you installed it into, then run `mimick`:
```installed locally and activate environment
pixi shell
mimick options... args...
```
2. be in the directory you installed Mimick into and call it with `pixi run mimick` instead of just `mimick`:
```installed locally and not activating environment
pixi run mimick options... args...
```

### Updating with pixi
Once again, slightly different depending on whether it was installed globally or locally:
```installed globally
pixi global update mimick
```
```installed locally
cd path/to/projectdir
pixi update mimick
```

#### **pip**
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

### Using with pip
Just call `mimick` from the command line
```call Mimick
mimick options... args...
```

### Updating with pip
You'll need to repeat the [Install with pip](#install-with-pip) process above.

<!-- tabs:end -->
