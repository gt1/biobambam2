#! /bin/bash
function precollated_fq_1
{
cat <<EOF
@A/1
ACGTACGT
+
HGFEDCBA
EOF
}

function precollated_fq_2
{
cat <<EOF
@A/2
GTCAGTCA
+
GHGFEDCB
EOF
}
function precollated_fq
{
	precollated_fq_1
	precollated_fq_2
}

function runopts
{
	../src/fastqtobam $* <(precollated_fq) | ../src/bamtofastq | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}

function runopts2
{
	../src/fastqtobam $* <(precollated_fq_1) <(precollated_fq_2) | ../src/bamtofastq | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}

function runoptsfile
{
	../src/fastqtobam $* <(precollated_fq) | ../src/bamtofastq filename=/dev/stdin | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}

function runoptsfile2
{
	../src/fastqtobam $* <(precollated_fq_1) <(precollated_fq_2) | ../src/bamtofastq filename=/dev/stdin | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}

function runoptsgz
{
	../src/fastqtobam $* <(precollated_fq) | ../src/bamtofastq gz=1 | zcat | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}

function runoptsgz2
{
	../src/fastqtobam $* <(precollated_fq_1) <(precollated_fq_2) | ../src/bamtofastq gz=1 | zcat | cmp <(precollated_fq)

	# copy pipe return status array
	PIPESTAT=( ${PIPESTATUS[*]} )

	if [ ${PIPESTAT[0]} -ne 0 ] ; then
		echo 'fastqtobam failed'
		return 1
	elif [ ${PIPESTAT[1]} -ne 0 ] ; then
		echo 'bamtofastq failed'
		return 1
	elif [ ${PIPESTAT[2]} -ne 0 ] ; then
		echo 'cmp failed'
		return 1
	else
		return 0
	fi
}


#
runopts ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runopts namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runopts namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runopts namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runopts namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

runopts2 ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runopts2 namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runopts2 namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runopts2 namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runopts2 namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

#
runoptsfile ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runoptsfile namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runoptsfile namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runoptsfile namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runoptsfile namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

runoptsfile2 ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runoptsfile2 namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runoptsfile2 namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runoptsfile2 namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runoptsfile2 namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

# gzip
runoptsgz ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runoptsgz namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runoptsgz namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runoptsgz namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runoptsgz namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

runoptsgz2 ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runoptsgz2 namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runoptsgz2 namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runoptsgz2 namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runoptsgz2 namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi


exit ${R}
