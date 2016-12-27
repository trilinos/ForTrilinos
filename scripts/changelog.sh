#!/bin/bash

[[ ! -f LICENSE ]] && echo "Call the script from ForTrilinos root" && exit 1

ARGS=(
    # --token XXX
    --no-pull-requests
    --include-labels 'bug,enhancement,New Feature'
    --enhancement-labels 'enhancement,New Feature'
)
github_changelog_generator Trilinos/ForTrilinos "${ARGS[@]}"
