
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
        <li><a href="#installation">Installation</a></li>
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
- Clone the repository with `git clone https://github.com/sambitmishra0628/PSP-GNM`
- On a terminal window, go to the PSP-GNM directory
- Create the conda environment using the provided yaml file: `conda env create -f psp_gnm_env.yaml` 
- Activate the environment: conda activate psp_gnm
- Checkout the command line arguments needed to run psp_gnm: python scripts/psp_gnm.py --help

### Prerequisites
You will need access to a Linux machine (or a MacOS or Windows 10 with Windows Subsystem for Linux enabled). The specific instructions below work best in a Linux (Ubuntu 18.02/Ubuntu20.04) platform

### Installation

1. Clone the repo
   ```sh
   git clone https://github.com/sambitmishra0628/repo_name.git
   ```
2. Install NPM packages
   ```sh
   npm install
   ```



<!-- USAGE EXAMPLES -->
## Usage

Use this space to show useful examples of how a project can be used. Additional screenshots, code examples and demos work well in this space. You may also link to more resources.

_For more examples, please refer to the [Documentation](https://example.com)_




<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.



<!-- CONTACT -->
## Contact

Your Name - [@twitter_handle](https://twitter.com/twitter_handle) - email

Project Link: [https://github.com/sambitmishra0628/repo_name](https://github.com/sambitmishra0628/repo_name)



<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* []()
* []()
* []()





