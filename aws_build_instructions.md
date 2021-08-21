# Instructions for building reports

## Setting up AMI

This procedure was tested on a `t2.large` AWS instance with Ubuntu 20.04.

### Setting up `R`

To install the most recent version of `R` from CRAN, these dependencies can be installed. 

```bash
# libraries needed to install over https
sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common

# add CRAN repository to your system sources’ list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# install R
sudo apt-get update 
# will take a couple minutes
sudo apt-get install -y r-base r-base-dev 

# install recent version of pandoc + citeproc
sudo apt-get install -y pandoc
sudo apt-get install -y pandoc-citeproc

# install texlive + extras
sudo apt-get install -y texlive
sudo apt-get install -y texlive-latex-extra
sudo apt-get install -y texlive-fonts-extra libfontconfig1-dev
```

### Setting up GitHub credentials

We will use the `R` package [`credentials`](https://github.com/r-lib/credentials) to set your GitHub credentials to access downloads from GitHub over https. First we need two additional libraries installed via `apt-get`.

```bash
sudo apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev libfontconfig1-dev
```

Log in to GitHub and click on your profile picture thumbnail in the upper right-hand corner. Click on Settings and in the left menu select "Developer Settings". In the new menu, select Personal access tokens. Generate a new token given it an arbitrary name. __Do not leave this page__. The access token is available only this once. Copy the access token to your clipboard.

Instead, head back to your running instance and start `R`. Install the credentials package.

```r
install.packages("credentials")
```

If you are prompted, type "yes" to create your local `R` package library. Or you can log into the session as `root`. 

Next, set your GitHub access token.

```r
credentials::set_github_pat()
```

At the prompt that says `Password for 'https://PersonalAccessToken@github.com':` paste your access token (it may not be visible after being pasted) and hit return. Check your credentials: 

```r
credentials::set_github_pat()
```

The output should read something like `Using GITHUB_PAT from ...` You can now clone the GitHub repository using `git`. Exit `R` and at the command line, `clone` the repository.

```bash 
git clone https://github.com/covpn/correlates_reporting
```

### Setting up the `R` environment using `renv`

Start `R` and install the `renv` and `here` packages.

```r
install.packages("renv")
install.packages("here")
```

Exit out of `R` and from the command line, move into the `correlates_reporting` directory and open an `R` session.

```bash
cd correlates_reporting && R
```

Once `R` has started, it should load the `renv` package automatically. When you get an `R` command prompt, execute the following command to download and configure all necessary `R` packages.

```R
renv::restore()
```

This will take about 45 minutes. During this time you will see messages like this.

````
Installing (some package) [version] ...
  OK [built from source]
````

### Installations for machine learning pipeline

To run the machine learning pipeline, `python3` must be available and installed.
On the Hutch servers, we are using version 3.8. If needed, this version of
python can be installed as follows.

```bash
# add repository to install older versions of python
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt-get update
# install 3.8
sudo apt-get install -y python3.8
# after installation check python version
python --version
```

Next, we use `pip` to install the `pytorch` library. We are using version 1.7.1,
without GPU acceleration. This version can be installed as follows.

```bash
# install pip if needed 
sudo apt-get install -y python3-pip
# install pytorch
pip install torch==1.7.1+cpu torchvision==0.8.2+cpu torchaudio==0.7.2 -f https://download.pytorch.org/whl/torch_stable.html
```

Finally, we need to create a virtual environment in the package repository and
install requirements.

```bash
# need to be in correlates_reporting/cor_surrogates folder
cd cor_surrogates
# install virtualenv if needed
pip install virtualenv
# create virtual environment
virtualenv ./pytorch
# install requirements
pip install -r requirements.txt
```

## Building the report

### Placing data file

The pipeline expects data to be formatted according to documentation [...]
(https://github.com/covpn/correlates_reporting/). A data set with this format
should be placed in `correlates_reporting/data_raw/[folder]`, where `[folder]`
is `data_raw_dir` specified in the `config.yml` file (see below).

### Workflow for building reports

#### Updating `renv`

Once an image has been built that includes the GitHub repository, the general work flow will be to
- `git pull` the updated `correlates_reporting` repository;
- open `R` in the `correlates_reporting` directory;
- run `renv::restore()` to update `R` package list;
- exit `R` and proceed to report building below.


#### Setting `config`

To change certain inputs to the pipeline, the `config.yml` file can be edited.
The head of this file is displayed below.

````
# the default analysis is moderna-like analysis
default: &default
  data_raw_dir: moderna
  subset_variable: None
  subset_value: All
  assays: [bindSpike, bindRBD, pseudoneutid50, pseudoneutid80]
  times: [B, Day29, Day57, Delta29overB, Delta57overB, Delta57over29]

moderna_real:
  <<: *default
  data_in_file: real_moderna_data.csv # Weiping can set
  study_name: COVE
  study_name_code: COVE
 ````

Here, we see several options that are set by default:
- `data_raw_dir`: the subfolder or `data_raw` where the raw data will be placed
- `subset_variable`: variable to subset the entire data set (e.g., by region or
baseline serostatus). Leave as `None` to run the pipeline on the entire data
set.
- `subset_value`: The value of `subset_variable` to subset the data to.
- `assays`: Names of assays to be included in the analysis pipeline.
- `times`: Times at which assays are measured

Additionally, these default options are inherited by the remaining
configurations. For example, the `moderna_real` configuration includes
additional variables:
- `data_in_file`: The name of the raw data file placed in `data_raw_dir`
- `study_name`: The name of the study to be printed in the report
- `study_name`: The name of the study that is accessed within the code.

Before running the pipeline, an environment variable needs to be set indicating
which configuration is active. This can be done at the command line using
`export`. For example,

```bash
export TRIAL=moderna_real
```

This makes the `moderna_real` configuration "active".

#### Immunogenecity report

The general workflow that we have developed relies on the following steps:
1. pre-processing the raw data;
2. creating immunogenecity figures and tables; and
3. compiling the immunogenecity report.

There are `make` commands for each of these three steps.

```bash
make immuno_report
```

The compiled `pdf` report will appear in `_report_immuno/covpn_correlates_immuno.pdf`.

The immunogenecity report is outcome blinded and may be used both to validate the veracity of the assay data and to inform other aspects of the correlates analysis. 

#### Baseline risk score report

The baseline risk score report can be generated using a similar `make` command. 

```bash
make risk_report
```

This compiled `pdf` report will appear in `_report_riskscore/covpn_correlates_report_riskscore.pdf`.

This command also saves an additional data object to the `data_clean` folder that is subsequently used in the correlates of risk reporting detailed below. Thus, __this report must be compiled prior to building the correlates of risk report__.


#### Tier 1 correlates report

After the immunogenecity report has been satisfactorily examined, a correlates report can be generated using `make` commands as follows.

```bash
make cor_report
```

The compiled `pdf` report will appear in `_report_cor/covpn_correlates_cor.pdf`.

For test builds, we recommend turning town the number of bootstrap and permutation resamples in `cor_coxph/code/params.R` by setting the values of the variables [`B`](https://github.com/CoVPN/correlates_reporting/blob/fb1e0c976e6ffb8ed939325dbd20a6c59f44f82b/cor_coxph/code/params.R#L2) and [`numPerm`](https://github.com/CoVPN/correlates_reporting/blob/fb1e0c976e6ffb8ed939325dbd20a6c59f44f82b/cor_coxph/code/params.R#L3) both to `5`.