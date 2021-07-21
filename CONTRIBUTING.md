# CoVPN Guide to working with the repository

This set of instructions describes a suggested workflow and best practices for
working smoothly across the CoVPN Correlates Report GitHub repository. The goal
is to ensure that we all follow the same process when merging new code into the
`master` branch intended for use by collaborators (including our future selves).
Ultimately, the aim is to make sure that our code works across as many instances
and compute setups as possible; thus, for all contributed code, the original
programmer is responsible for resolving merge conflicts prior to requesting code
review from other members.

## Git, GitHub, Branching, Pull Requests

### Branching

Git is a powerful collaborative tool: its major strength is that multiple people
can work across the same source repository at the same time, with the ability to
seamlessly combine their work. This is done through several steps: branching,
committing changes, and generating pull requests.

_Branching_ is the process of creating a distinct line of work (i.e., a sequence
of commits) in a local copy of the repository. Your work on a branch does _not_
impact the state of the `master` branch. This has the advantage of allowing you
to develop and test new code, without any disruptions to the stable version of
the code stored on the shared `master` branch. This way, the programmer is able
to make as many changes as desired, then commit and push the updated code to
make it available to other programmers; however, the new code will not impact
the `master` branch until it is manually requested to be added in a pull
request.

### Pull Requests

A [_pull request_](https://help.github.com/articles/about-pull-requests/) (PR)
is a feature specific to GitHub that allows a programmer to request that changes
in a given branch be _pulled_ into the `master` branch (or any other target
branch). Upon _opening_ a new PR, a few things can happen. First, GitHub will
check that there are no conflicts with the `master` branch (see later),
automatically flagging any cases of such conflicts. Next, an opportunity to
request a review is made available. In a review of a PR, the programmer is able
to tag others to check their programming logic, a chance to ensure that the
underlying logic is clear and that the code included in the PR runs on more than
the machine on which it was written. The interface GitHub provides is very nice:
reviewer can point to specific areas in the code with questions and comments, or
even suggest changes to individual lines of code. As you work, make sure you've
created a new branch, commit often with clear commit messages, and generate a PR
once the code is ready to be reviewed and merged.

### Conflicts

An important downside of branching and pulling branches into `master` is that
two programmers can unwittingly work on the same file across multiple branches.
Should a PR updating a particular file (or files) also changed on your branch be
approved and merged into `master` before the merging of your PR, Git will be
unsure of which version of the updated file ought to be used. This results in
a _merge conflict_: Git cannot be sure whether to accept your changes or those
made by another programmer (on a different branch). This can happen for many
reasons. The best method for avoiding such issues is simply for each programmer
to make changes only to those files that they need to edit and that such changes
be explicit and contained.

If, despite your best efforts, a merge conflict does occur when making a PR, the
conflict(s) will need to be evaluated and resolved manually prior to merging the
PR into the `master` branch. _Merge conflicts should be resolved by the original
programmer_, as they are most familiar with the code being contributed in the PR
in question. Once all changes to be included in your PR have been committed on
your branch, `checkout` the `master` branch (i.e., `git checkout master`) and
pull any/all changes from the remote (i.e., `git pull REMOTE master`). Now,
switch back to your branch (`git checkout MYBRANCH`) and attempt to merge in the
up-to-date `master` branch (i.e., `git merge master`). If there are any
conflicts, a warning message will appear in the command line, indicating the
presence of merge conflicts requiring manual resolution. Use `git status` to
identify the files with the conflicts, or if you are using RStudio, look at the
Git pane (orange boxes with a "U" in them will appear before the file name).

For each of the files identified to have a merge conflict, read through the code
-- the differences should be highlighted. Sometimes, these will be large chunks
of code, other times just single lines. Differences start with `<<<<<<< BRANCH
NAME`, followed by all code from your branch, and then then several equal signs
(`=======`) followed by the code from the `master` branch, which itself ends
with `>>>>>>> master`. Select which version of the code, or blend of them, to
keep and delete the remainder of the code. There may be multiple conflicts in
a file, so make sure that all are found and resolved. After correcting each
conflict across all affected files, ensure the code still runs as intended.
Commit the edits, push, and voila...the PR should no longer indicate any merge
conflicts.

### Rules for pull requests

Please adhere to the following rules before submitting a pull request.
- As noted above, do not submit a PR that has merge conflicts.
- Do not add compiled graphics and PDFs to your pull requests (with the
possible exception of verification documentation). Pull requests should only be
code. If you need to store compiled files locally, add them to `.gitignore`.
- Update the `Makefile` to include relevant documentation for how to
execute the code and dependencies between scripts.
- Include a rule `all` at the top of the `Makefile` that compiles all of the
output that is needed to build the report for your subdirectory.
- Include a rule `clean` that removes all compiled graphs and data sets from
your subdirectory.
- Confirm that after executing your `make clean` command (to remove any
compiled results from your local directory), your `make all` command executes on
your system.

### Subdirectories

Git has an interesting approach to tracking how folder structures behave: it
does not look at the directory structure of the repository but instead tracks
files and their paths. This means that folders that do not contain files are not
added to git for tracking; therefore, such files will not be added or shared to
the remote repository. To get around this, add a simple empty file `.gitkeep` to
the folder. Now, git has a file to track and will add the folder to version
control.

This is important because if your code expects a folder to exist and does not
check first if it needs to create it, the code will stop running. Since the goal
of this project is to run seamlessly, adding the `.gitkeep` file will make this
possible with minimal effort.

`.gitkeep` can be created a few different ways. Through the command line, use
`touch /path/to/folder/.gitkeep`; or through R, use
`file.create("/path/to/folder/.gitkeep")`.

### Git Resources

Here are a few useful notes and readings that may be useful in remedying any
problems that may arise when working with `git`. These range the spectrum from
applied to fairly technical:

* ["Version Control with `git`" (Software
    Carpentry)](https://swcarpentry.github.io/git-novice/): a comprehensive
    introduction to the uses of `git` and social coding with GitHub. This is
    aimed towards students and research professionals.
* ["Happy `git` and GitHub for the useR" (Jenny Bryan,
    RStudio)](http://happygitwithr.com/): a comprehensive introduction to both
    the inner workings of `git` and best practices in using GitHub, with a focus
    on integrating these tools into your R workflow. Aimed towards
    masters-level university students.
* ["Tools for Reproducible Research" (Karl
    Broman)](http://kbroman.org/Tools4RR/): a semester-long course on best
    practices for modern reproducible research, including a couple of lectures
    on `git` and GitHub.
    * ["Version Control with `git` and
        GitHub"](http://kbroman.org/Tools4RR/assets/lectures/04_git.pdf)
* ["Introduction to `git`" (Berkeley's Stat 159/259, Fall 2015, KJ
    Millman)](http://www.jarrodmillman.com/rcsds/standard/git-intro.html): a
    thorough walkthrough of some common `git` commands as well as explanations
    of what these commands do when invoked.

Some basic commands:

To clone the repo, 

```bash
git clone https://usr@github.com/CoVPN/correlates_reporting.git
```




## Tracking the Package Environment

In an effort to promote computational reporoducibility, we are making an effort
to _all use the same version of R and R packages_, using {renv}. When opening
the project in R, {renv}'s tracking system should start automatically and
highlight steps to follow. Things to note include

- Make sure that all R packages being used are the version specified in the
  {renv} lockfile (`renv.lock`), ensuring consistency across all programmers'
  work by calling `renv::restore()`. _Do this every time you manually start
  working,_ frequently updating from the `master` branch to ensure that the
  lockfile is itself up-to-date.
- To add a new package, there are a few steps. First, ensure there are no
  differences in your installed packages and the ones specified in the lockfile
  by calling `renv::restore()`. This will enforce the correct versions. Next,
  install any package(s) as usual (e.g., `install.packages()` or
  `remotes::install_github()`) and integrate use of the package into your code
  as necessary. Finally, call `renv::snapshot()`. This will record use of the
  new package(s) into the package library, adding to the lockfile.

Also take note of [the {renv} collaboration
guide](https://rstudio.github.io/renv/articles/collaborating.html).

## Configuring the analysis

To accommodate code developed to span multiple analyses (e.g., for Moderna and
for Janssen), we make use of the [`config`]
(https://cran.r-project.org/web/packages/config/vignettes/introduction.html)
package. The file `config.yml` in the base repository directory contains
specifications pertaining to each analysis. This includes specifications for 
the mock and real Moderna data analysis, as well as several
specifications for Janssen's pooled and region-specific analyses. More
configurations may be added in the future.

To load an analysis configuration, all the user must do is set an environment
variable in a shell session. This can be achieved via the command line as
follows.

```bash
export TRIAL=moderna_mock
```

```tcsh
setenv TRIAL moderna_mock
```

On Windows, replace export with set in command line prompt. This can also be achieved in an interactive `R` session as follows.

```r
Sys.setenv(TRIAL = "moderna_mock")
```

Note that if the `TRIAL` variable __has not been set__, the data processing
script __will fail__. Users should get in the habit of setting this variable
before executing any R code. Once `TRIAL` *has* been set, configurations can be 
loaded into an R session using `config::get(config = Sys.getenv("TRIAL"))`. 
However, in general a user will not need to explicitly call `config::get()` in
their code, __so long as `_common.R` is
`source`'d into the code__. 

Examining `_common.R`, we see that `config::get` is called in that script and
arguments in the `config` list are read into the global environment of that `R`
session (so users do not need to refer to e.g., `config$study_name` and instead
can simply write `study_name`).

## Automated Testing and Continuous Integration

In order to facilitate _reproducible reporting_, this project makes use of
automated testing and checking via the continuous integration service
[Travis-CI](https://travis-ci.com/github/CoVPN/correlates_reporting/builds).
For any contributed or changed code, the Travis-CI service will automatically
attempt to build the reports corresponding to each of the statistical analyses
outline in the statistical analysis plan. This automatic check serves to ensure
that all dependencies and analyses, as well as the reports themselves, can be
executed/generated from a bare machine instance (Linux Ubuntu, by default). Any
contributed code that fails to pass these checks will be deemed unsuitable for
merging into the project; thus, it is imperative that contributors not only
verify their updated code by building reports locally (via the `Make` recipes
provided) but also check that their contributions pass on the machine instances
spawned by the Travis-CI service. For any given build, logs are provided to help
elucidate any potential sources of error. Upon passing a continuous integration
check, Travis-CI will automatically post the resultant reports to the `gh-pages`
branch at https://github.com/CoVPN/correlates_reporting/tree/gh-pages.

## Optimal Flow

With all this being said, if we all work together, we can try to minimize cases
in which we need to resolve merge conflicts. _Work in your own folder, make
changes specific and manageable, and commit/push regularly._

Follow this flow closely:
1. Check out `master` if not there already. If using RStudio, select the `master`
   branch; if on the command line: `git checkout master` will do.
2. Pull any updates from the GitHub repo. If in RStudio, click the "pull"
   button; if on the command line, simply `git pull`.
3. Make a new branch, check it out, and push to GitHub. If using RStudio, make
   a new branch via the GUI. If on the command line: `git branch NAME`, then
   `git checkout NAME`, then `git push -u origin NAME` (n.b., the first two
   of these may be combined into `git checkout -b NAME`).
4. Make sure that the packages you are using match those in the shared library
    by calling `renv::restore()`. Resolve any issues as necessary.
5. Write your code and commit regularly. Install packages as needed. Alternate
   between `renv::restore()` and `renv::snapshot()` to capture changes to the
   active library. To add folders that do not have content in them, or whose
   content are ignored through `.gitignore`, add a `.gitkeep` to the folder and
   add the folder to be tracked. If using RStudio, use the git pane and add the
   folder; if on the command line, use `git add /path/to/folder/.gitkeep`.
   Commit the file to complete the addition to tracking the folder through git.
6. At the beginning of starting work on a branch for the day or reverting back
   to prior work, make sure as much of the code is up-to-date with the remote.
   Make sure all your edits are committed, `checkout master`, pull updates from
   the remote `master` branch, `checkout` your branch and merge `master` into
   your branch. Resolve conflicts as necessary.
7. When preparing to open a PR, make sure all your edits are committed. Check
   out the `master` branch, pull in changes to the remote, check out your
   branch and merge `master` into your branch. Resolve conflicts. Commit
   and push any updates.
8. Go to GitHub and open a PR. Identify the updates clearly, and tag reviewers,
   who will provide feedback on how to improve/streamline the contributed code.
   Note that the PR will be subject to automated continuous integration checks
   by the [Travis-CI
   service](https://travis-ci.com/github/CoVPN/correlates_reporting/builds); a
   PR that fails to pass these automated checks _will not be merged_.
9. Once the PR has been accepted, merge normally (i.e., do _not_ use a special
   merge), then delete the feature branch.
10. On your machine, clean things up by checking out `master` and pulling.
    Delete the local copy of that working branch (`git branch -d NAME`), and
    prune your remote copies (`git remote prune origin`).
