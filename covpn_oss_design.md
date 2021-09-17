---
title: Open Source Design of the CoVPN Immune Correlates Reporting Repository
author: Nima Hejazi and David Benkeser
---

The immune correlates analysis of the joint COVID-19 Vaccine Prevention Network
(CoVPN) and US Government (USG) statistics team is a large-scale data science
project, with a distributed team of collaborators, implementing the publicly
available USG/CoVPN immune correlates statistical analysis plan (SAP), available
at https://doi.org/10.6084/m9.figshare.13198595. From its inception, the project
has placed a strong emphasis on transparency, computational reproducibility,
manual source code verification, and portability across modern computing
infrastructures. In line with these objectives, the data analysis and reporting
architecture adheres closely to best practices developed in the free and open
source software community; the unrestricted source code is publicly available as
a single version-controlled repository on GitHub, a widely used collaborative
programming platform, at https://github.com/CoVPN/correlates_reporting.

The repository includes distinct modules, each mapped closely to a section of
the SAP, allowing all members of the implementation team to develop code for
components of the analysis closely aligned with their respective areas of
expertise. This modular structure ensures that multiple, distinct analyses can
be developed simultaneously, while minimizing the frequency with which analysts
face impasses arising from development work outside of their analysis module.
This modular structure has important benefits for manual code review, an
external robustness check, as well, by allowing for detailed feedback and
suggested edits to be provided (via GitHub's pull requests feature) throughout
the development cycle. Within each module, analysis code is developed in the `R`
language and environment for statistical computing (https://www.r-project.org/),
with the `RMarkdown` format (https://rmarkdown.rstudio.com/) used to write
module-specific reports and the `bookdown` package (https://bookdown.org/) used
to compile these reports into the larger, distributed reports, thematically
organized around immunogenicity and correlates of risk/protection analyses. The
dependency structure across analysis modules is encoded and controlled by GNU
Make (https://www.gnu.org/software/make/), which simplifies the report building
process into single commands. The GitHub repository is tightly coupled with
Travis-CI (https://travis-ci.org/), a continuous integration service that
automatically builds the reports in a small cloud computing environment upon
proposed changes to the repository's stable source code. By integrating these
tools into the data science team's workflow, the repository design ensures that
the analyses are -- both individually and collectively -- free of and robust to
programming errors, transparently developed, and readily produce analysis
reports containing independently verified results.

The CoVPN immune correlates analysis and reporting repository represents, to our
knowledge, the first entirely publicly developed and available data science
framework for the analysis of any vaccine efficacy trial. Incorporating modern
statistical computing, software engineering, and open source software tooling
and best practices allows the repository to embody a public, persistent, and
continuous peer review of the statistical methodology and scientific results of
the CoVPN/USG COVID-19 immune correlates analyses. Beyond serving as a reference
implementation of the CoVPN/USG COVID-19 immune correlates SAP, the success of
the repository in facilitating large-scale, distributed collaboration on the
evolving analysis of vaccines critical to curbing the current global pandemic
-- while closely adhering to strict standards of modularity, reproducibility,
and transparency -- demonstrates the great scientific value that can be added by
the careful integration of such tooling. While many improvements can be made in
future iterations, the repository serves as a model for how modern data science
tools and open science standards can guide statistical analysis development and
reporting in future vaccine trials.
