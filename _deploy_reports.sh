#!/bin/bash

# configure your name and email if you have not done so
git config --global user.email "benkeser@emory.edu"
git config --global user.name "David Benkeser"
git config --global http.postBuffer 100000000

# clone the repository's gh-pages branch
git clone -b gh-pages \
  https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git \
  correlates_reporting

# overwrite contents from existing gh-pages branch
cd correlates_reporting
echo "Files in correlates_reporting/ _before_ copying:"
ls -l

# replace with reports and note R version
if [ "$REPORT_TYPE" == "IMMUNO" ]
then
  echo "copying Immuno report"
  ls -s $TRAVIS_BUILD_DIR/_report_immuno/*
  cp -rf $TRAVIS_BUILD_DIR/_report_immuno/* ./
elif [ "$REPORT_TYPE" == "COR" ]
then
  echo "copying COR report and Risk Score"
  ls -s $TRAVIS_BUILD_DIR/_report_riskscore
  ls -s $TRAVIS_BUILD_DIR/_report_cor
  cp -rf $TRAVIS_BUILD_DIR/_report_riskscore/* ./
  cp -rf $TRAVIS_BUILD_DIR/_report_cor/* ./
elif [ "$REPORT_TYPE" == "COP" ]
then
  echo "copying COP report and Risk Score"
  ls -s $TRAVIS_BUILD_DIR/_report_riskscore
  ls -s $TRAVIS_BUILD_DIR/_report_cop
  cp -rf $TRAVIS_BUILD_DIR/_report_riskscore/* ./
  cp -rf $TRAVIS_BUILD_DIR/_report_cop/* ./
elif [ "$REPORT_TYPE" == "RISK" ]
then
  echo "copying RISK report"
  ls -s $TRAVIS_BUILD_DIR/_report_riskscore
  cp -rf $TRAVIS_BUILD_DIR/_report_riskscore/* ./
fi
echo "Reports built with R version $TRAVIS_R_VERSION"

# check what files have been copied to branch gh-pages
echo "All files in correlates_reporting/ _after_ copying:"
ls -l

# stage, commit, push copied files to branch gh-pages
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]
then
  COMMIT_MESSAGE="Update reports via ${TRAVIS_COMMIT}."
else
  COMMIT_MESSAGE="Update reports: PR #${TRAVIS_PULL_REQUEST} ($TRAVIS_COMMIT)."
fi
git add --all *
git commit -m "${COMMIT_MESSAGE}"
git push -q origin gh-pages
