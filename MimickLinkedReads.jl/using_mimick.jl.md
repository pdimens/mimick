# Using Mimick's Julia port
These instructions are _for now_.

## 1. Install Julia
Install Julia by whichever way is recommended by your system.

## 2. Install MimickLinkedReads.jl
Until this package matures a bit and I figure out how I want the API to look, the
most sensible way to install it would be:

1. Clone the Mimick repository
```bash
git clone https://github.com/pdimens/mimick.git
```
2. Change to the `mimick` directory and swap to the `mimick.jl` branch
```bash
cd mimick

git checkout mimick.jl
```

3. Start a Julia REPL session. You can do this however feels natural: using the 
standard REPL (which is great), Jupyter, VScode + Julia addon, etc. Julia needs
to be invoked with a thread number, so use `-t X` where `X` is the number of 
threads you would like it to use.

```bash
julia -t 15
```

4. Within Julia, use the `Pkg` shell to install the local version of `MimickLinkedReads.jl`. This can be done one of two ways:

    - **option 1**: using the `Pkg` shell by pressing `]` on an empty line and using the `dev` command (idiomatic method)
    ```julia
    julia>]
    (@v1.11) pkg> dev path/to/MimickLinkedReads.jl
    ```
    - **option 2**: using `Pkg` explicitly (Jupyter Notebook method)
    ```julia
    julia> using Pkg
    julia> Pkg.develop("path/to/MimickLinkedReads.jl")
    ```

5. Load `MimickLinkedreads.jl` into the userspace. Don't include the `.jl` part.
```julia
julia> using MimickLinkedReads
```

Now you should have access to the two exported `mimick()` methods. At any point, you can press `?` on an empty line to switch the REPL into help mode:
```julia
julia>?
help?> mimick
```

You can also press a semicolon `;` on an empty line to switch to shell mode, where you can now quickly use shell commands in a quick-and-dirty bash terminal:
```julia
julia>;
shell> ls
```

Note: A backspace on an empty line reverts the REPL back into regular mode.

## The `mimick` command
Here is the gist of the `mimick` command (all of it's described in the `?` docstring). It's "unzipped" here for visual clarity.
```julia
mimick(
    "fasta.fa",
    "file.vcf",
    "standard:haplotagging",
    prefix = "simulated/",
    coverage = 10,
    n_molecules = 2,
    mol_cov = 0.2,
    mol_len = 80000,
    insert_length = 500,
    insert_stdev = 50,
    read_len = [150,150],
    singletons = 0.35,
    circular = false,
    attempts = 25,
    seed = 0
)
```