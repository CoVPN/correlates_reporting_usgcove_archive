#! /bin/bash

# configure your name and email if you have not done so
git config --global user.email "benkeser@emory.edu"
git config --global user.name "David Benkeser"
git config --global http.postBuffer 100000000

# clone the repository
git clone -b gh-pages \
  https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git \
  correlates_reporting

# remove contents from existing gh-pages branch
cd correlates_reporting
git rm -rf *
echo "All files in /correlates_reporting after git rm"
ls -l 
# replace with contents from master branch /website
cp -r ../_report_cor/* ./
# move tmp_lectures and tmp_homeworks in and rename
cp -r ../_report_immuno/* ./

echo "All files in /correlates_reporting after copies"
ls -l 

COMMIT_MESSAGE="update the test build reports."
git add --all *
git commit -m "${COMMIT_MESSAGE}"
git push -q origin gh-pages

