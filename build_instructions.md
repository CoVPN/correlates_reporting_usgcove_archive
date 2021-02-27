# Building the correlates report on AWS

## AMI

This procedure was tested on a `t2.large` instance with Ubuntu 20.04.

## Setting up `R`

To install the most recent version of `R` from CRAN, these dependencies can be installed. 

```bash
# libraries needed to install over https
sudo apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common

# add CRAN repository to your system sources’ list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# install R
sudo apt-get update 
sudo apt-get install r-base r-base-dev 
```

## Setting up GitHub credentials

We will use the `R` package [`credentials`](https://github.com/r-lib/credentials) to set your GitHub credentials to access downloads from GitHub over https. First we need two additional libraries installed via `apt-get`.

```bash
sudo apt-get install libssl-dev libcurl4-openssl-dev libxml2-dev
```

Log in to GitHub and click on your profile picture thumbnail in the upper right-hand corner. Click on Settings and in the left menu select "Developer Settings". In the new menu, select Personal access tokens. Generate a new token given it an arbitrary name. __Do not leave this page__. The access token is available only this once. Copy the access token to your clipboard.

Instead, head back to your running instance and start `R`. Install the credentials package.

```r
install.packages("credentials")
```

Next, set your GitHub access token.

```r
credentials::set_github_pat()
```

At the prompt that says `Password for 'https://PersonalAccessToken@github.com':` paste your access token (it may not be visible after being pasted) and hit return. Check your credentials: 

```r
credentials::set_github_pat()
```

The output should read something like `Using GITHUB_PAT from ...` You can now clone the GitHub repository using `git`.

```bash 
git clone https://github.com/covpn/correlates_reporting
```

## Setting up the `R` environment using `renv`

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

You now have all of the software needed to build the reports.

## Placing data file

The pipeline expects data to be formatted according to documentation [...](https://github.com/covpn/correlates_reporting/). A data set with this format should be placed in `correlates_reporting/raw_data`.

## Making the immunogenecity report

To build the immunogenecity report you can execute the command

```bash
make immuno_report
```




