
# Recommended reading / tutorials

* You should be comfortable using the unix shell. There is a good introduction [here](https://swcarpentry.github.io/shell-novice/).
* The pipeline is driven by a Makefile. [Here](http://byronjsmith.com/make-bml/) is a tutorial for how to use `make` to manage a data analysis pipeline.
* I do all my version control with `git`. You won't need advanced knowledge of `git` to run the pipeline but if you want to change it, it would be a good idea to track your changes with `git`. There is a good tutorial [here](https://swcarpentry.github.io/git-novice/).

# System requirements
You can run the pipeline on macOS or linux. Windows is not supported.

* Python 3.6+
* GNU make
* git

You can check whether a program is installed with the command `which`:
```bash
which python3
which make
which git
```

It should print something like:
```bash
$ which python3
/Users/Aaron/miniconda3/bin/python3
$ which make
/usr/bin/make
$ which git
/usr/local/bin/git
```

If you don't see a the path to the executable for all of the programs, you will need to install them.

On macOS you can install `make` and `git` with Command Line Developer Tools:
```bash
xcode-select --install
```

You can install Python3 with Homebrew following the instructions [here](https://docs.python-guide.org/starting/install3/osx/).

You can check your Python version with:
```bash
python3 --version
```

# Install the pipeline code
```bash
git clone --recurse-submodules https://github.com/troyanr/VAMP_analysis
```

# Running the pipeline
1\. First move the uncompressed `.fastq` files to `VAMP_analysis/data`

2\. Setup the project environment
```bash
cd h1_tiled_library_analysis # navigate to the project directory
make code # setup project virtual python environment
```

3\. Run the pipeline!
```bash
make -j data
```

# How it works
The pipeline is driven by GNU Make. The `Makefile` manages the dependencies between input, intermediate, and final data files. It also automatically sets up a python virtual environment and installs all the required python packages.

## Editable parameters
Depending on the analysis, the end user will need to modify parameters in 2 places:

* On top of the `Makefile` are settings for the filtering Q score and mismatch filters
* The `params.yaml` file has settings for aggregating the filtered sequence and generating sequence logos.
