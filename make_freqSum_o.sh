
options=$(getopt -l "config:withName,packed,out:" -o "c:npo:" -- "$@")
eval set -- "${options[@]}"

CONFIG_FILE=""
WITH_NAME=false
PACKED=false
OUT_FILE=""

while [ -n "$1" ]
do
    case "$1" in
        -c|--config)
            shift
            CONFIG_FILE=$1
        ;;
        -n|--withName)
            WITH_NAME=true
        ;;
        -p|--packed)
            PACKED=true
        ;;
        -o|--out)
            shift
            OUT_FILE=$1
        ;;
        --)
            shift
            break
        ;;
        *)
            echo "$1 is not an option"
            exit 1
        ;;
    esac
    shift
done

if [[ ${CONFIG_FILE} == "" ]]
then
    echo "Missing configuration file!"
    exit 2
fi

PASTE_CMD="paste"

N=`wc -l < ${CONFIG_FILE}`
for i in `seq 1 ${N}`
do
    if ${WITH_NAME}
    then
        read FILE_DIR FIELD NAME <<< `sed -n "${i}p" ${CONFIG_FILE}`
    else
        read FILE_DIR FIELD <<< `sed -n "${i}p" ${CONFIG_FILE}` # read the n-th row in CONFIG_FILE
    fi
    if [[ ${FILE_DIR##*.} == "geno" ]]
    then
        CMD="cut -c ${FIELD} ${FILE_DIR} | sed 's/9/-1/g'"
    else
        CMD="cut -f ${FIELD} ${FILE_DIR}"
    fi
    #echo ${CMD}
    PASTE_CMD+=" <(${CMD})"
done

if ${PACKED}
then
    PASTE_CMD+=" | awk -v OFS='\t' '{a[\$0]++} END {for (i in a) print i,a[i]}'"
fi

if [[ ${OUT_FILE} != "" ]]
then
    PASTE_CMD+=" > ${OUT_FILE}"
fi

#echo ${PASTE_CMD}
eval ${PASTE_CMD}
