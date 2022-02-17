
# Protein Stability Prediction with a Gaussian Network Model <br> (PSP-GNM)


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
Protein stability prediction upon point mutation is a key biological problem. Single point mutations can alter protein function resulting in disease incidence. A very significant example is that of sickle cell anemia, wherein a single genomic mutation results in a single amino acid change and impairs the function of hemoglobin. It is therefore essential to develop methods that can predict the impact of point mutations on protein stability. More specifically, such methods should enable the estimation of free energy change (ddG) upon point mutation and be able to tell us whether a mutation reduces, increases or doesn't change the thermodynamic stability of proteins. 

### Built With

* [Python3.8](https://www.python.org/downloads/release/python-380/)
* [scikit-learn=1.0.1](https://scikit-learn.org/stable/whats_new/v1.0.html#version-1-0-1)
* [scipy=1.7.2](https://scipy.github.io/devdocs/release.1.7.2.html)
* [statsmodels=0.13.1](https://www.statsmodels.org/stable/release/version0.13.1.html)
* [biopython=1.79](https://biopython.org/docs/1.79/api/Bio.html)
* [click=8.0.3](https://click.palletsprojects.com/en/8.0.x/)
* [numpy=1.21.4](https://numpy.org/devdocs/release/1.21.4-notes.html)
* [seaborn=0.11.2](https://seaborn.pydata.org/installing.html)



<!-- GETTING STARTED -->
## Getting Started
You will need access to a Linux machine (or a MacOS or Windows 10 with Windows Subsystem for Linux enabled). The specific instructions below work best in a Linux (Ubuntu 18.02/Ubuntu20.04) platform

To get a local copy up and running follow these simple steps:
- Install miniconda: 

```
wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.11.0-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

- Clone the repository with `git clone https://github.com/sambitmishra0628/PSP-GNM`
- On a terminal window, go to the PSP-GNM directory
- Create the conda environment using the provided yaml file: `conda env create -f psp_gnm_env.yaml` 
- Activate the environment: conda activate psp_gnm
- Checkout the command line arguments needed to run psp_gnm: `python scripts/psp_gnm.py --help`

The following help manual should be displayed

```
Usage: psp_gnm.py [OPTIONS]

Options:
  --data_file TEXT    Name of the .csv file containing the information on ddG
                      for the mutants  [required]
  --outfile TEXT      Name of the file to     which the PSP-GNM-calculated
                      energies and experimental energies will be written
                      [required]
  --outdir TEXT       Name of the directory to     which the intermittent
                      result files will be written to  [required]
  --wt_pdb_dir TEXT   Directory containing the wild type atomic pdb files
                      [required]
  --num_jobs TEXT     Maximum number of jobs to be run in parallel  [required]
  --dist_cutoff TEXT  Distance cutoff for interactions in GNM  [default: 9;
                      required]
  --num_modes TEXT    Number of modes to be used  [default: 10; required]
  --help              Show this message and exit.
```

### Prerequisites
You will need access to a Linux machine (or a MacOS or Windows 10 with Windows Subsystem for Linux enabled). You will need to have conda installed.


<!-- USAGE EXAMPLES -->
## Usage

There are primarily 2 scripts included in the `scripts/` directory:
- `scripts/psp_gnm_benchmark_data.py`
- `scripts/psp_gnm.py`

`psp_gnm_benchmark_data.py` is written to be run on the benchmark data (`datasets/`)
`psp_gnm.py` is the generic version that is written to be run on any generic data

An example run using the benchmark data is given below. Make sure that your current working directory is the `PSP_GNM` folder having the `scripts` directory.

```
python psp_gnm_benchmark_data.py --datafile test_data/S350_test_benchmark_run --outdir S350_test_run_output --outfile S350_test_benchmark_run_out.csv --wt_pdb_dir test_data/pdb_test --num_jobs 4 --dist_cutoff 9 --num_modes 10
```
In the above:
  `num_jobs` should be the number of cores in your machine minus 1

This run will create an output directory `S350_test_run_output` and store all the intermediate files containing information on the contacts broken during partial unfolding.

It will then create S350_test_benchmark_run_out.csv containing the calculated ddG values.



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Sambit Mishra - sambit.mishra@nih.gov

Project Link: [PSP-GNM](https://github.com/sambitmishra0628/PSP-GNM)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





