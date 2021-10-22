# CoVPN/USG Correlates Analysis Reporting [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5593129.svg)](https://doi.org/10.5281/zenodo.5593129)

__Note:__ As of 16 October 2021, this repository has been archived. Generalized
versions of its data processing and analysis modules have been split in two and
migrated to a [correlates data processing
workflow](https://github.com/CoVPN/correlates_processing/) and a [correlates of
risk and protection analysis
workflow](https://github.com/CoVPN/correlates_reporting2/). This archive serves
as a reference of the workflows used in developing reports for the immune
correlates analyses of the Moderna and Janssen (ENSEMBLE) COVID-19 vaccine
efficacy trials.

The _Statistical Analysis Plan_ is available
  * [on Figshare](https://doi.org/10.6084/m9.figshare.13198595)
  * [on Overleaf](https://www.overleaf.com/project/5ecd5bcc18e1d30001c913ec)

## Collaboration Guide

Please consult [this blog
post](https://davidbphd.com/project-organization-for-reproducible-data-science/),
which outlines most aspects of our project organization recommendations.

* [Code style guide](https://style.tidyverse.org/), with some modifications;
  this will largely be enforcd with [`styler`](https://styler.r-lib.org/).
* Project organization: _mostly_ independent subdirectories, each incorporating
  [`here`](https://here.r-lib.org/) for path resolution.
* Package version control and virtual environments using
  [`renv`](https://rstudio.github.io/renv/).
* Code review procedure: see our [contribution
   guidelines](https://github.com/CoVPN/correlates_reporting_usgcove_archive/blob/master/CONTRIBUTING.md).

---

## CITATION

When citing the analysis workflow or analytic results produced by its use,
please cite the following

        @software{gilbert2021usgcove,
          author = {Gilbert, Peter B and Fong, Youyi and Benkeser, David and
            Hejazi, Nima S and Hughes, Ellis and Borate, Bhavesh and Yu,
            Chenchen and Lu, Yiwen and Li, Kendrick Q and {van der Laan}, Lars
            WP and Simpkins, Brian},
          title = {{COVID-19 Prevention Network Immune Correlates Analyses},
          year  = {2021},
          doi = {10.5281/zenodo.5593129},
          url = {https://github.com/CoVPN/correlates_reporting_usgcove_archive}
        }

---

## License

The contents of this repository are distributed under the GPL-3 license. See
file [`LICENSE.md`](https://github.com/CoVPN/correlates_reporting_usgcove_archive/blob/master/LICENSE.md)
for details.
