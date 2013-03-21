#! /bin/bash
VERSION=`grep AC_INIT < configure.in | awk -F',' '{print $2}'`
FIRST=`echo $VERSION | awk -F'.' '{print $1}'`
SECOND=`echo $VERSION | awk -F'.' '{print $2}'`
THIRD=`echo $VERSION | awk -F'.' '{print $3}'`
NEXTTHIRD=`expr ${THIRD} + 1`

awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^AC_INIT/ {gsub(first"."second"."third,first"."second"."third+1);print} ; !/^AC_INIT/{print}' < configure.in | \
	awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^LIBRARY_VERSION=/ {gsub("="first"."third"."second,"="first":"third+1":"second);print} ; !/^LIBRARY_VERSION=/{print}' \
	> configure.in.tmp
mv configure.in.tmp configure.in

pushd debian
dch --distribution UNRELEASED -v ${FIRST}.${SECOND}.${NEXTTHIRD}-1
popd

git add debian/changelog
git add configure.in
git commit
git push

git tag -a biobambam_${FIRST}_${SECOND}_${NEXTTHIRD} -m "biobambam version ${FIRST}_${SECOND}_${NEXTTHIRD}"
git push origin biobambam_${FIRST}_${SECOND}_${NEXTTHIRD}
