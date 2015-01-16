#! /bin/bash
LIBMAUSVERSION=$1

if [ -z ${LIBMAUSVERSION} ] ; then
	exit 1
fi

sed -i -e "s/libmaus >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmaus >= ${LIBMAUSVERSION}/" configure.ac
sed -i -e "s/libmausdigests >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmausdigests >= ${LIBMAUSVERSION}/" configure.ac
git add configure.ac
git commit -m "bump libmaus version"
git push

git checkout experimental-debian
sed -i -e "s/libmaus-dev (>= [0-9][0-9]*.[0-9][0-9]*.[0-9][0-9]*)/libmaus-dev (>= ${LIBMAUSVERSION})/" debian/control
git add debian/control
git commit -m "bump libmaus version"
git push

git checkout experimental
