#! /bin/bash
VERSION=`grep AC_INIT < configure.ac | awk -F',' '{print $2}'`
FIRST=`echo $VERSION | awk -F'.' '{print $1}'`
SECOND=`echo $VERSION | awk -F'.' '{print $2}'`
THIRD=`echo $VERSION | awk -F'.' '{print $3}'`
NEXTTHIRD=`expr ${THIRD} + 1`

awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^AC_INIT/ {gsub(first"."second"."third,first"."second"."third+1);print} ; !/^AC_INIT/{print}' < configure.ac | \
	awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^LIBRARY_VERSION=/ {gsub("="first"."third"."second,"="first":"third+1":"second);print} ; !/^LIBRARY_VERSION=/{print}' \
	> configure.ac.tmp
mv configure.ac.tmp configure.ac

pushd debian
export DEBEMAIL=gt1@sanger.ac.uk
export DEBFULLNAME="German Tischler"
# dch --distribution unstable -v ${FIRST}.${SECOND}.${NEXTTHIRD}-0
dch --distribution unstable -v ${FIRST}.${SECOND}.${NEXTTHIRD}
dch --release
# dch --release -v ${FIRST}.${SECOND}.${NEXTTHIRD}-1
popd

git add debian/changelog
git add configure.ac
git commit
git push

TAG=biobambam_experimental_${FIRST}_${SECOND}_${NEXTTHIRD}
git tag -a ${TAG} -m "biobambam experimental version ${FIRST}_${SECOND}_${NEXTTHIRD}"
git push origin ${TAG}
