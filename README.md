# vaccine-speed-vs-prioritization

This repository contains the code accompanying Agarwal, Komo, Patel, Pathak, and Ünver (2021). It is a VSEIR model implementation to simulate epidemiological implications of various scenarios of vaccine rollout policy, rollout speed, and mitigation policy.

### Repository Organization

This repository connects two pipelines to produce exhibits from the raw data. The process can be thought of as: ``./raw_data/`` -- via ``./code/R/`` -- > ``./inputs/`` -- via ``./sir_code/`` -- > ``./ouput`` and ``./exhibits``

To run from start to finish, run the shell script ``make.sh`` (in the root directory, type ``bash make.sh``). Please see the replication section to download the raw data and check the dependencies section to ensure the appropriate installations are in place.

### Replication

1. Clone this repository into your local root directory
2. Download the raw data as a compressed zip file from [this](https://www.dropbox.com/s/nk2menrpx7k159b/raw_data.zip?dl=0) link
3. Save the downloaded unzipped raw data folder in the local root directory as ``./raw_data``
4. Check the dependencies section to ensure the appropriate languages and packages are installed
5. Run the shell script ``make.sh``. The estimated total run time is approximately two days.

### Dependencies

This project was written and tested in R (4.0.0), Matlab 2020b and Artleys' KNITRO 12.2.0.

The R required packages (version used in codebase testing in parentheses) are:
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) (1.0.2)
* [stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0) (1.4.0)
* [feather](https://cran.r-project.org/web/packages/feather/index.html) (0.3.5)
* [assertr](https://cran.r-project.org/web/packages/assertr/index.html) (2.7)
* [tidyr](https://cran.r-project.org/web/packages/tidyr/index.html) (1.1.2)
* [data.table](https://cran.r-project.org/web/packages/data.table/index.html) (1.13.4)
* [zoo](https://cran.r-project.org/web/packages/zoo/index.html) (1.8.8)
* [openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html) (4.2.2)
* [survey](https://cran.r-project.org/web/packages/survey/index.html) (4.0)
* [SASxport](https://cran.r-project.org/web/packages/SASxport/index.html) (1.7.0)

The MATLAB installation requires the interface with the KNITRO solver. 

### Maintainers
* Nikhil Agarwal (agarwaln@mit.edu)
* Andrew Komo (komo@mit.edu)
* Chetan Patel (cpatel@mit.edu)
