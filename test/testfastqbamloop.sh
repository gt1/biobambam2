#! /bin/bash
../src/fastqtobam precollated.fq | ../src/bamtofastq | cmp precollated.fq

# copy pipe return status array
PIPESTAT=( ${PIPESTATUS[*]} )

if [ ${PIPESTAT[0]} -ne 0 ] ; then
	echo 'fastqtobam failed'
	exit 1
elif [ ${PIPESTAT[1]} -ne 0 ] ; then
	echo 'bamtofastq failed'
	exit 1
elif [ ${PIPESTAT[2]} -ne 0 ] ; then
	echo 'cmp failed'
	exit 1
fi

exit 0
