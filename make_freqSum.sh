
options=$(getopt -l "config:,packed,out:" -o "c:po:" -- "$@")
eval set -- "${options[@]}"

CONFIG_FILE=""
PACKED=false
OUT_FILE=""

while [ -n "$1" ]
do
    case "$1" in
        -c|--config)
            shift
            CONFIG_FILE=$1
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
    ROW=`sed -n "${i}p" ${CONFIG_FILE}` # read the n-th row in CONFIG_FILE
    FILE_DIR=`cut -f 1 <<< ${ROW}`
    FIELD=`cut -f 2 <<< ${ROW}`
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
