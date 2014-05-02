#! /bin/bash
source dupsingle.sh
source dupsinglemarked.sh

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
