# Chimera Evolve

This repository contains the source code for Chimera Evolve, an evolutionary algorithm that designs a single CDS for multiple organisms using the [Chimera ARS score](https://academic.oup.com/bioinformatics/article/31/8/1161/212401).


<hr>

## Building

The algorithm is written in Rust. After downloading and installing Rust ([instructions here](https://www.rust-lang.org/tools/install)), you can generate an executable to run the algorithm as follows:

    git clone git@github.com:intbio-ncl/chimera_evolve.git
    cd chimera_evolve
    cargo build --release
    
This creates an executable appropriate for your system in `chimera_evolve/target/release/`. For example, on Linux, from the `chimera_evolve` directory, you can then run `./target/release/chimera-evolve --help`:


    David Skelton <d.j.skelton@newcastle.ac.uk>
    Optimises a single coding sequence for multiple organisms using the Chimera ARS score

    USAGE:
        chimera-evolve [OPTIONS] <protein> <cds>... --outfile <outfile>

    FLAGS:
        -h, --help       Prints help information
        -V, --version    Prints version information

    OPTIONS:
        -c, --crossovers <crossovers>         Sets the number of crossover events to carry out per generation [default: 100]
        -s, --gen_start <generation_start>    Sets the generation start size. [default: 200]
        -g, --generations <generations>       Sets the number of generations to run the algorithm for [default: 1000]
        -q, --method <method>                 Method to use to score solutions (weighted, min) [default: min]
        -m, --mutations <mutations>           Sets the number of mutation events to carry out per generation [default: 300]
        -o, --outfile <outfile>               Name of the file to which result will be written
        -w, --weights <weights>               Comma separated list of weights if weighted method is used

    ARGS:
        <protein>    Sets the protein to optimise a CDS for
        <cds>...     Coding sequences from organisms to optimise for
