#! /bin/bash
SCRIPTDIR=`dirname "${BASH_SOURCE[0]}"`
pushd ${SCRIPTDIR}
SCRIPTDIR=`pwd`
popd

source ${SCRIPTDIR}/dupsinglemarked.sh
source ${SCRIPTDIR}/dupsinglemarkedsortedqreset.sh

function reset
{
	dupsinglemarked | ../src/bamsort SO=queryname | ../src/bamreset resetsortorder=0
}

./bamcmp <(reset) <(dupsinglemarkedsortedqreset)

if [ $? -ne 0 ] ; then
	exit 1
fi

SORTORDER=unknown
for i in `reset | ./bamtosam | head -n 1 | egrep "^@HD"` ; do
	if [ ! -z "`echo $i | egrep \"^SO:\"`" ] ; then
		SORTORDER=`echo $i | perl -p -e "s/^SO://"`
	fi
done

if [ "${SORTORDER}" != "queryname" ] ; then
	exit 1
fi

# check reset by column
COL2=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $2}' | sort -u`
COL2CNT=`echo $COL2 | wc -l`

if [ "${COL2CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL2}" != "4" ] ; then
	exit 1
fi

COL3=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $3}' | sort -u`
COL3CNT=`echo "${COL3}" | wc -l`

if [ "${COL3CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL3}" != "*" ] ; then
	exit 1
fi

COL4=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $4}' | sort -u`
COL4CNT=`echo "${COL4}" | wc -l`

if [ "${COL4CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL4}" != "0" ] ; then
	exit 1
fi

COL5=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $5}' | sort -u`
COL5CNT=`echo "${COL5}" | wc -l`

if [ "${COL5CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL5}" != "0" ] ; then
	exit 1
fi

COL6=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $6}' | sort -u`
COL6CNT=`echo "${COL6}" | wc -l`

if [ "${COL6CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL6}" != "*" ] ; then
	exit 1
fi

COL7=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $7}' | sort -u`
COL7CNT=`echo "${COL7}" | wc -l`

if [ "${COL7CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL7}" != "*" ] ; then
	exit 1
fi

COL8=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $8}' | sort -u`
COL8CNT=`echo "${COL8}" | wc -l`

if [ "${COL8CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL8}" != "0" ] ; then
	exit 1
fi

COL9=`reset | ./bamtosam | egrep -v "^@" | awk -F '\t' '{print $9}' | sort -u`
COL9CNT=`echo "${COL9}" | wc -l`

if [ "${COL9CNT}" -ne 1 ] ; then
	exit 1
fi

if [ "${COL9}" != "0" ] ; then
	exit 1
fi

exit 0
