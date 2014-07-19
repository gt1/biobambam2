#! /bin/bash

# update branches
git checkout experimental
git pull
git checkout experimental-debian
git pull
git checkout debian
git pull
git checkout master
git pull

# merge
git checkout master
git merge experimental
git push

git checkout debian
git merge experimental-debian
git push

# add release tag/branch
git checkout master
VERSION=`grep <configure.ac "AC_INIT" | perl -p -e "s/.*AC_INIT\(//" | awk -F ',' '{print $2}'`
DATE=`date +"%Y%m%d%H%M%S"`
RELEASE=${VERSION}-release-${DATE}
git checkout -b ${RELEASE}-branch master
PATH=/software/hpag/autotools/bin:${PATH} autoreconf -i -f
ADDFILES="INSTALL Makefile.in aclocal.m4 autom4te.cache compile config.guess config.h.in config.sub configure depcomp install-sh ltmain.sh m4/libtool.m4 m4/ltoptions.m4 m4/ltsugar.m4 m4/ltversion.m4 m4/lt~obsolete.m4 missing src/Makefile.in"
mv .gitignore .gitignore_
git add ${ADDFILES}
git commit -m "Release ${RELEASE}"
mv .gitignore_ .gitignore
git tag ${RELEASE}
git push origin ${RELEASE}
git checkout master
git branch -D ${RELEASE}-branch
git checkout experimental
