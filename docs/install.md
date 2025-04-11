# Installing Mimick

There are three ways you can install Mimick:
1. using conda
2. using pixi
3. local pip installation


## Install with conda
!> WIP: mimick is not yet available on conda.


To install Mimick using conda (or miniconda, mamba, etc.) into an existing environment, it's as simple as:
```bash
# you're already in the environment
conda install -c conda-forge bioconda::mimick

# you're not in the environment
conda install -n env_name bioconda::mimic
```
where `env_name` is the name of the existing conda environment
you want to install Mimick into. Alternatively, you can create
a new conda environment and install Mimick into it using one command:
```bash
# global environment
conda create -n mimick -c bioconda -c conda-forge mimick

# local environment
conda create -p path/to/workdir/mimick -c bioconda -c conda-forge mimick
```

### Using with conda
You'll need to activate the environment you installed Mimick into, then call it with:
```bash
mimick options... args...

# start here
mimick --help
```

### Updating with conda
If you wish to update Mimick, just replace word `create`/`install` in the commands above with `update`, e.g.:
```bash
# you're already in the environment
conda update -c conda-forge bioconda::mimick

# you're not in the environment
conda update -n env_name -c conda-forge bioconda::mimick
```

## Install with pixi
Pixi is a new ultra-fast Rust-based environment manager much like conda is. If you haven't tried it out, you ought to give it a shot.
With it, you can install Mimick to be accessible in your PATH, i.e. a "global" installation:
```bash
pixi global install -c conda-forge -c bioconda mimick
```
!> Make sure ~/.pixi/bin is in your PATH
```bash
export PATH=~/.pixi/bin:$PATH
```

If you prefer a local environment installation into a project directory:
```bash
pixi init -c conda-forge -c bioconda projectname && cd projectname
pixi add mimick
```

### Using with pixi
If installed globally, this is very similar to the conda installation:
```bash
# installed globally
mimick options... args...
```

If you installed Mimick locally using Pixi, you can do one of two things:
1. activate the pixi shell in the directory you installed it into, then run `mimick`:
```bash
pixi shell
mimick options... args...
```
2. be in the directory you installed Mimick into and call it with `pixi run mimick` instead of just `mimick`:
```bash
pixi run mimick options... args...
```

### Updating with pixi
Once again, slightly different depending on whether it was installed globally or locally:
```bash
# installed globally
pixi global update mimick

# installed locally
cd path/to/projectdir
pixi update mimick
```

## Install with pip
If neither conda nor pixi appeal to you, Mimick can be installed with pip too. To do that, you will need to first
download [the latest release](https://github.com/pdimens/mimick/releases) (recommended) or clone [the git repository](https://github.com/pdimens/mimick):
```bash
# download the latest release (replace x.x.x with the actual version #)
wget -O mimick.tar.gz https://github.com/pdimens/mimick/releases/download/x.x.x/mimick.x.x.x.tar.gz
tar -xvzf mimick.tar.gz

# clone the repository
git clone https://github.com/pdimens/mimick.git
```

Then, you will need to enter the Mimick directory and do a local `pip` installation:
```bash
cd mimick
pip install .
```

### Using with pip
Just call `mimick` from the command line
```bash
mimick options... args...
```

### Updating with pip
You'll need to repeat the [Install with pip](#install-with-pip) process above.