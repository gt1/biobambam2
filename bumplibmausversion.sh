#! /bin/bash
LIBMAUS2VERSION=$1

if [ -z ${LIBMAUS2VERSION} ] ; then
	exit 1
fi

sed -i -e "s/libmaus2 >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmaus2 >= ${LIBMAUS2VERSION}/" configure.ac
sed -i -e "s/libmaus2digests >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmaus2digests >= ${LIBMAUS2VERSION}/" configure.ac
sed -i -e "s/libmaus2irods >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmaus2irods >= ${LIBMAUS2VERSION}/" configure.ac
sed -i -e "s/libmaus2seqchksumsfactory >= [0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*/libmaus2seqchksumsfactory >= ${LIBMAUS2VERSION}/" configure.ac
git add configure.ac
git commit -m "bump libmaus2 version"
git push

git checkout experimental-debian
sed -i -e "s/libmaus2-dev (>= [0-9][0-9]*.[0-9][0-9]*.[0-9][0-9]*)/libmaus2-dev (>= ${LIBMAUS2VERSION})/" debian/control
git add debian/control
git commit -m "bump libmaus2 version"
git push

git checkout experimental
