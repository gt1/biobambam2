#! /bin/bash
function precollated_fq
{
cat <<EOF
@A/1
ACGTACGT
+
HGFEDCBA
@A/2
GTCAGTCA
+
GHGFEDCB
EOF
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

runopts ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
runopts namescheme=generic ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi
echo "This is expected to fail:"
runopts namescheme=c18s ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
echo "This is expected to fail:"
runopts namescheme=c18pe ; R=$? ; if [ ${R} -eq 0 ] ; then exit 1 ; fi
runopts namescheme=pairedfiles ; R=$? ; if [ ${R} -ne 0 ] ; then exit ${R} ; fi

exit ${R}
