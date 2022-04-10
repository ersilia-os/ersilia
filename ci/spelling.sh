#!/usr/bin/env bash

npm install -g cspell
git fetch origin master:master
git diff --name-only --diff-filter=AM master $HEAD | grep '\.md$' | xargs --no-run-if-empty -L1 npx cspell -c ./ci/spelling-config.json