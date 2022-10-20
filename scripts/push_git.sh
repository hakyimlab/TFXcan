#!/bin/bash

git pull
git add {*.r,*.R,*.ipynb,*.py,*.qmd,*.sh,*.rmd,*.Rmd,*.py~,*.html,*.tex,*.pdf}
#git reset -- "$( ls -l */ )"
git commit -m 'Updated scripts'
git push

git filter-branch -f --index-filter 'git rm --cached --ignore-unmatch scripts/*_cache/*/*.rdb'