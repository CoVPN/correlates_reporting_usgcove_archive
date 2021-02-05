## CoVPN Guide to working with Git and {renv}

These are a set of instructions that are suggested to be followed to smooth the process of working on the CoVPN Correlates Report github repository. 
The goal of this is to make sure we all follow the same process when we go to merge our new code into the master branch that is intended for use.
Additionally, we want to make sure that our code works across as many instances as possible, and that means that the original programmer should be responsible for ensuring there are no merge conflicts before requesting other members perform code review.

### Git, Github, Branching, and PRs

#### Branching

The strength of using git as a collaborative tool is that multiple people can we working across the same repository at the same time, and seamlessly combine their work. This is done through a process of branching and pull requests. 

Branching is creating a copy of the repository for you to work in, without impacting the master branch. This allows the programmer to make as many changes as they want, commit and push to make it available to other programmers, but will not impact the master branch until it is requested to be added in a Pull Request (PR).

#### Pull Requests

A Pull Request (PR) is a feature in github where the programmer of the branch _requests_ the changes made on their branch be _pulled_ into the master branch. This gives an opportunity for a few things to happen. First, github will check that there are no conflicts (covered later), and will flag cases where there are. Next, it gives an opportunity to request a review. This is where we tag another programmer to check our programming logic, and ensure it runs on more than your machine. The interface github gives too is very nice in that the reviewer can point to specific areas in the code with questions, etc.

#### Conflicts

However, the downside of branching and pulling branches into master is that someone can unwittingly work on the same file in another branch, have their pull request approved before yours, and all of a sudden git is not sure which version to use. This is a merge conflict: Git is not sure whether to accept your changes or the changes another person made. This can happen for any number of reasons. The best method to avoid this issue is to make sure everyone only makes changes to the files they need to, and that the changes are explicit and contained. 

If despite your best efforts, a merge conflict occurs when making a pull request, the conflicts will need to be resolved manually before merging. This should be done by the original programmer, since they know the expectations. Once all the changes on your branch have been committed, checkout the master branch (`git checkout master`) and pull(`git pull`) all the changes. Checkout your branch again(`git checkout MYBRANCH`) and attempt to merge in the master branch (`git merge master`). If there are conflicts, a warning message in the command line will appear that there were merge conflicts that need to be resolved. Use `git status` to identify the files with the conflicts, or if you are using RStudio, look at the git pane and orange boxes with a "U" in them will appear before the file name. 

For each of the files identified to have a merge conflict, read through the code, and it should highlight the differences. Sometimes they will be large chunks of code, and other times single lines. Differences start with `<<<<<<< BRANCH NAME`, which all code following is the code in your branche, then several equal signs (`=======`) then the code from the master branch, which ends with `>>>>>>> master`. Select which version, or a blend of both, to keep, and delete the other code. There may be multiple conflicts in a file, so make sure to have found and resolved each one. After doing this for each file, ensure the code still runs as intended. Commit the edits and voila, the pull request should no longer show merge conflicts.

### Tracking Environment

For this project, in an effort to ensure we are all using the same packages, we use {renv}. When opening the project, it should automatically start {renv}'s tracking system and highlight steps to follow. Things to note include:

	- Make sure that all packages being used are the version specified from {renv} lockfile to ensure consistency across all programmers work by calling `renv::restore()`. Do this every time you manually start working.
	
	- To add a new package, there are a few steps. First, ensure there are no differences in your installed packages and the ones specified in the lockfile by calling `renv::restore()`. This will enforce the correct versions. Next, install the package and apply it into your code as necessary. Call `renv::snapshot()` now. This will identify that the new package has been added and used in the code and add it to the lock file.

### Optimal Flow

With all this being said, if we all work together, we can try to minimize cases that we need to resolve merge conflicts. work in your own folder, make changes specific and manageable, and commit/push regularly.

Follow this flow:
1. checkout master if not there already (If in RStudio, Select master branch. If in command line: `git checkout master`)
2. pull any updates from github repo(If in RStudio, click "pull" button. If in command line: `git pull`)
3. make a new branch, check it out, and push to github. (If in RStudio, just make new branch. If in command line: `git branch NAME`,`git checkout NAME`,`git push -u origin NAME`).
4. Make sure that the packages you are using are the same as other by calling `renv::restore()`. Resolve any issues from that as necessary. 
5. Write your code and commit regularly. Install packages as needed. Alternate between `renv::restore()` and `renv::snapshot()` to capture active library use. 
6. At the beginning of starting working on a branch for the day or reverting back to prior work, make sure as much of the code is up to date with remote. Make sure all your edits are committed, check out master, pulling updates, check out your branch and merging master into your branch. Resolve conflicts as necessary. 
7. When preparing to open a pull request, make sure all your edits are committed. Check out the master branch, pull in the remote changes, check out your branch and merge master into your branch. Resolve conflicts. Commit and push any updates.
8. Go to github and open a PR. Identify the updates clearly, and tag reviewers. 
9. Once PR is accepted, merge normally (Do not use a special merge), and delete the branch.
10. On your machine, clean things up by checking out master and pulling. delete the local copy of that working branch (`git branch -d NAME`), and prune your remote copies (`git remote prune origin`).