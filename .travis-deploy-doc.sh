#!/bin/bash -x
# Copying and distribution of this file, with or without modification,
# are permitted in any medium without royalty provided the copyright
# notice and this notice are preserved.  This file is offered as-is,
# without any warranty.
set -e
ABS_DEPLOY_KEY="`pwd`/.travis_ci_gh_pages_deploy_key"
if [[ -r "$ABS_DEPLOY_KEY" ]]; then
    echo "Deployment key exists, attempting to upload"
    if [[ -z "${DEPLOY_DOC_TO_REPOSITORY}" ]]; then
        DEPLOY_DOC_TO_REPOSITORY="${TRAVIS_REPO_SLUG}"
    fi
    chmod 600 "$ABS_DEPLOY_KEY"
    export GIT_SSH_COMMAND="ssh -v -i $ABS_DEPLOY_KEY -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no"
    rm -Rf gh-pages
    git clone --depth 1 git@github.com:${DEPLOY_DOC_TO_REPOSITORY}.git --depth 1 --branch=gh-pages gh-pages
    BUILT_DOCS_DIR=`cd docs/build/html && pwd`
    cd gh-pages
    rm -Rf ./${DEPLOY_DOC_TO_DIRECTORY}/*
    mkdir -p ./${DEPLOY_DOC_TO_DIRECTORY}
    cp -R $BUILT_DOCS_DIR/* ./${DEPLOY_DOC_TO_DIRECTORY}/
    git add --all .
    git config user.name "Travis CI"
    git config user.email "nobody@example.org"
    if git commit -m "Automatic upload of documentation built from ${TRAVIS_COMMIT}"; then
	git push origin gh-pages
    fi
fi
/usr/bin/killall python2 || true
