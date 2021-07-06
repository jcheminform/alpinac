# ALPINAC

ALgorithmic Process for Identification of Non-targeted Atmospheric Compounds

## Pre-print
**Automated fragment identification for electron ionisation mass spectrometry:
application to atmospheric measurements of halocarbons**,
Myriam Guillevic, Aurore Guillevic, Martin Vollmer, Paul Schlauri, Matthias
Hill, Lukas Emmenegger, Stefan Reimann.
[hal-03176025](https://hal.inria.fr/hal-03176025),
[arxiv 2103.13807](https://arxiv.org/abs/2103.13807).

## First release on June 29, 2021 with LGPL licence, Copyright EMPA and Inria, 2019 - 2021.

## How to use git on both Windows and Linux platforms
On Windows, type:
```
git config --global core.autocrlf true
```
On Linux, type:
```
git config --global core.autocrlf input
```

## Installation requirements
On any platform: install Anaconda.

### Anaconda installation on Windows
(`rdkit` is not strictly necessary).

1.   Install Anaconda.
2.   Install the `rtdkit env` in Anaconda:
     http://www.rdkit.org/docs/Install.html
3.   Then activate this environment: `conda activate my-rdkit-env`.
4.   Then, with this environment activated in the Conda prompt, install Spyder at this environment (even
     if spyder is installed in the base):
     `conda install spyder` (**not** `pip install spyder`! this would not work).
5.   Launch Spyder.
     Now you can import the `rdkit` in Spyder.

You may have to install via the conda prompt:
```shell
conda install matplotlib
conda install h5py
conda install -c conda-forge lmfit
pip install pysmiles
```
  
Install the _Networkx_ package in the Conda prompt:
```sell
conda install -c anaconda networkx
```

### Linux

```shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install python3 python3-numpy python3-scipy python3-matplotlib ipython3 python3-notebook python3-pandas python3-sympy python3-nose python3-h5py python3-networkx
sudo apt-get install python3-lmfit
```
The alternative installation tool is `pip3`, installing it with:
```shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install pip3
```
Then installing additional packages with `pip3`:
```shell
pip3 install pysmiles
```

### Pysmiles package
For converting SMILES chemical formulas, we use the package **pysmiles**:
[Github link](https://github.com/pckroon/pysmiles)
Installation:
```shell
pip3 install pysmiles
```
Usage:
```python
import pysmiles
from pysmiles.read_smiles import read_smiles
from pysmiles.write_smiles import write_smiles
nx_graph_molecule = read_smiles(smiles, explicit_hydrogen=True, zero_order_bonds=True, reinterpret_aromatic=True)
write_smiles(nx_graph_molecule, default_element='*', start=None)
```

## Where to start

The main file is `mode_identification.py`.
You can open it with e.g. Spyder3 to have a look inside.

It can be called in two different ways.

### Command line
[`mode_identification.py`](./mode_identification.py) can be executed through command line, providing arguments. Command-line
examples are provided in
[`mode_identification_script_train.sh`](./mode_identification_script_train.sh)
and [`mode_identification_script_val.sh`](mode_identification_script_val.sh).

The first compulsory argument is the filename, preceeded by the directory.

Additional arguments can be provided with:
- `--target-elements`: chemical elements potentially present, e.g. CHFClBr. The order has no importance.

- `--target-formula`: the chemical formula of the molecular ion. ALPINAC will use only this formula or smaller fragments. For example with CCl4 given, ALPINAC will use only C, Cl, CCl, CCl2, CCl3, CCl4 and Cl2 (and all possible isotopologues).

### Pre-defined Python script
The file 
[`mode_identification_script.py`](./mode_identification_script.py) is configured to call `mode_identification.py` for each file of the training set, using default parameters. You can either type its name in an anaconda prompt, making sure you are in its folder, or open and run it using e.g. Spyder3.

## Input test data

Data to test ALPINAC are provided in `data\nontarget_screening\fragments`. Input data are text files with the following columns:
1. `RT`: Retention time in seconds (unused)
2. `mass`: Measured mass, m/z
3. `mass_u_ppm`: Uncertainty of the measured mass (measurement noise), ppm, 1 sigma
4. `mass_cal_u_ppm`: Uncertainty of the mass calibration, ppm, 1 sigma
5. `area`: Measured signal intensity
6. `area_u`: Uncertainty of the measured intensity (unsued)
7. `peak_width`: Width of the mass peak, 1 sigma
8. `peak_alpha`: Fraction of Lorenzian peak shape, use zero if your peak is Gaussian.
9. `compound_bin`: Index of time bin, all co-eluting masses should have the same value.
10. `LOD`: Limit of detection for the measured mass (unsued)
11. `SN`: Ratio signal/noise (unused)
12. `Ionisation`: Ionisation type (unused), note that ALPINAC is developed for EI only.

A few columns are not used yet but may be taken into account in future developments.

## Output data

ALPINAC produces a text file with results for the identification of the masses provided as input. The text file is saved at `data\nontarget_screening\formulas`. The text file has the following format:

1. Percentage of signal for the chemical formula, compared to the sum of measured signal
2. Retention time
3. Percentage of signal compared to the largest peak - the largest peak has a value of 100.
4. Measured mass, m/z
5. Reconstructed exact mass, m/z
6. Combined uncertainty of the measured mass as computed from Eq. (1) in the paper, ppm
7. Mass difference between measured and exact masses, ppm
8. Reconstructed chemical formula
9. DBE (double bound equivalent)
10. Absolute intensity assigned to this chemical formula
11. Fraction of the measured peak assigned to this chemical formula
12. Ionisation type


## Organisation

We describe where each box of Figure 1 of the paper is coded.
First the input data and chemical data.

-   **Input: experimental data**. A file with uncertainties is required either as
command-line argument when running `python3 mode_identification.py <FILE.txt>`,
or with a wrapper script (see `.sh` files, or 
[`mode_identification_script.py`](./mode_identification_script.py).
-   **Exact mass and valence of atoms**: this is provided in the file
    [`periodic_table.py`](./periodic_table.py).
-   **Environmental abundance of isotopes**: this is also provided in the file
    [`periodic_table.py`](./periodic_table.py).
-   **Output**: a graph of mass spectrum with labeled peaks (for the identified ones).

Measured data are read and saved in the class structure `Compound` defined in the file `compound.py`.

We now describe each step of Figure 1 and provide the road map of the code.

-   **Step 1: Knapsack**. The functions are in `utils_identification.py`.
-   **Step 2: Initialise directed graph**. The functions are methods of the
    classes in `utils_graph`.py. 
-   **Step 3: Initialise isotopocule sets**. Functions are methods of the class
    `Isotopologues` defined in the file `isotopologue.py`.
-   **Step 4: Compute maximum contribution of each set**.
-   **Step 5: Rank Fragments according to likelihood estimator**.
-   **Step 6: Select most likely fragments**.
-   **Step 7: Optimise multiple isotopocule sets**. Calls `lmfit` Python package.
-   **Step 8: Update directed graph, remove singletons**.
-   **Step 9: Eliminate not optimised sets**.
-   **Step 10: Reconstruct most likely molecular ion(s)**.


## How to read the plotted mass spectra

-  Orange bar means identified fragment
-  Pink bar means non-identified fragment
-  Brown bar: fragment expected because of isotopes, but not found in measurements.




## Some measured masses have no solutions. What can I do?

1. Make sure the mass in question is a real signal, not e.g. electronic noise. If the peak is ionised twice, it will not be identified. Doulbe ionisation is not supported by ALPINAC.

2. Using default constant values in `const_identification.py`, run ALPINAC once. Based on the results, make your own best estimation of which chemical elements are present.

3. Run ALPINAC a second time using, as additional constrain, your chosen list of present atoms. You can also either increase the value of `max_opt_signal` to e.g. 99%, or decrease the value of `max_ratio_not_opt_masses` to e.g. 0.01.

4. If this does not work, do as in 3. but at the same time, increase the mass uncertainty `m_u_k` to e.g. 3.0.

## How can I be sure that the results are correct?

- We have observed that in case only the chemical elements C, H, N, and O are present, with a mass resolution of approx. 3500, results are sometimes inaccurate. We strongly advice to run ALPINAC several times while giving different possible combinations of chemical elements, e.g., "CH", "CHO", "CHN", "CHNO". 

- If ALPINAC suggests one chemical element is present in some formulae but you think this is wrong: this can happen, especially if no rare-isotope fragment was detected. Run ALPINAC one more time with a given list of chemical elements that do not contain the one you think is wrong, and check if ALPINAC can still reconstruct the same percentage of measured signal.

- With a mass resolution of approx. 3500, it can also be difficult to correctly reconstruct the number of hydrogen atoms. This number can be overestimated. We advice to run ALPINAC with a given chemical formula as constrain (in `mode_identification_script.py` or through the command line). Try multiple times with formulae that contain a different number of hydrogen atoms.

## ALPINAC is still running after 3 minutes. What can I do?

This can happen if the mass uncertainty is high, or if there are many measured masses, causing a number of knapsack formulae above 1000.

1. Use a lower value for `max_opt_signal` in `const_identification.py`, e.g. 80%.

2. Run ALPINAC using only a limited list of present chemical elements.

## My data are not measured with EI but with a softer ionisation. Can I use ALPINAC?

The main difference between EI and softer ionisation (such as chemical ionisation) is that no adduct formation occurs in an EI source. In practice, this means that we can assume that the formula of any detected mass has a non-negative DBE (double bound equivalent). This assumption is used to speed up the enumeration of the knapsack.

If your measured data contains adducts with formulae of negative DBE, this is not supported by ALPINAC. The correct formulae will not be reconstructed.

	






