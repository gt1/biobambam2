#! /bin/bash
VERSION=`grep AC_INIT < configure.ac | awk -F',' '{print $2}'`
FIRST=`echo $VERSION | awk -F'.' '{print $1}'`
SECOND=`echo $VERSION | awk -F'.' '{print $2}'`
THIRD=`echo $VERSION | awk -F'.' '{print $3}'`
NEXTTHIRD=`expr ${THIRD} + 1`

function cleanup
{
	if [ ! -z "${COMMITFILE}" ] ; then
		if [ -f "${COMMITFILE}" ] ; then
			rm -f "${COMMITFILE}"
		fi
	fi
}

COMMITFILE=commit_msg_$$.txt

trap cleanup EXIT SIGINT SIGTERM

# make sure we have the latest version
git checkout experimental
git pull

# create commit log message
joe "${COMMITFILE}"

if [ ! -s "${COMMITFILE}" ] ; then
	echo "Empty commit log, aborting"
	exit 1
fi

# update to next minor version
awk -v first=${FIRST} -v second=${SECOND} -v third=${THIRD} '/^AC_INIT/ {gsub(first"."second"."third,first"."second"."third+1);print} ; !/^AC_INIT/{print}' < configure.ac > configure.ac.tmp
mv configure.ac.tmp configure.ac

# commit file
git add configure.ac
git commit -F "${COMMITFILE}"
git push

# switch to experimental debian branch
git checkout experimental-debian
git pull
git merge experimental

# create change log message
pushd debian
export DEBEMAIL=gt1@sanger.ac.uk
export DEBFULLNAME="German Tischler"
dch --distribution unstable -v ${FIRST}.${SECOND}.${NEXTTHIRD}
dch --release
popd
git add debian/changelog
git commit -F "${COMMITFILE}"
git push

# back to experimental branch
git checkout experimental

TAG=biobambam2_experimental_${FIRST}_${SECOND}_${NEXTTHIRD}
git tag -a ${TAG} -m "biobambam2 experimental version ${FIRST}_${SECOND}_${NEXTTHIRD}"
git push origin ${TAG}

exit 0
