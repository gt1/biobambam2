#! /bin/bash
SCRIPTDIR=`dirname "${BASH_SOURCE[0]}"`
pushd ${SCRIPTDIR}
SCRIPTDIR=`pwd`
popd

source ${SCRIPTDIR}/dupsingle.sh
source ${SCRIPTDIR}/dupsinglemarked.sh

function runmark
{
	dupsingle | ../src/bammarkduplicates
}

function runmark2
{
	dupsingle | ../src/bammarkduplicates
}

./bamcmp <(runmark) <(dupsinglemarked)

if [ $? -ne 0 ] ; then
	exit 1
fi

./bamcmp <(runmark2) <(dupsinglemarked)

if [ $? -ne 0 ] ; then
	exit 1
fi
