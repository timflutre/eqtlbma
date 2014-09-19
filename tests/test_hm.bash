#!/usr/bin/env bash

# Aim: launch a basic functional test for eqtlbma_hm
# Author: Timothée Flutre
# Not copyrighted -- provided to the public domain

progVersion="1.0"

#------------------------------------------------------------------------------

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
    msg="\`${0##*/}' launches a basic functional test for eqtlbma_hm.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="  -b, --p2b\tabsolute path to the 'eqtlbma_bf' binary\n"
    msg+="  -H, --p2h\tabsolute path to the 'eqtlbma_hm' binary\n"
    msg+="  -R, --p2R\tabsolute path to the 'functional_tests.R' script\n"
    msg+="  -n, --noclean\tkeep temporary directory with all files\n"
    msg+="\n"
    msg+="Report bugs to <eqtlbma-users@googlegroups.com>."
    echo -e "$msg"
}

# Display version and license information on stdout.
function version () {
    msg="${0##*/} ${progVersion}\n"
    msg+="\n"
    msg+="Written by Timothée Flutre.\n"
    msg+="\n"
    msg+="Not copyrighted -- provided to the public domain\n"
    echo -e "$msg"
}

# http://www.linuxjournal.com/content/use-date-command-measure-elapsed-time
function timer () {
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local startRawTime=$1
        endRawTime=$(date '+%s')
        if [[ -z "$startRawTime" ]]; then startRawTime=$endRawTime; fi
        elapsed=$((endRawTime - startRawTime)) # in sec
        nbDays=$((elapsed / 86400))
        nbHours=$(((elapsed / 3600) % 24))
        nbMins=$(((elapsed / 60) % 60))
        nbSecs=$((elapsed % 60))
        printf "%01dd %01dh %01dm %01ds" $nbDays $nbHours $nbMins $nbSecs
    fi
}

# Parse the command-line arguments.
# http://stackoverflow.com/a/4300224/597069
function parseCmdLine () {
    getopt -T > /dev/null # portability check (say, Linux or Mac OS?)
    if [ $? -eq 4 ]; then # GNU enhanced getopt is available
	TEMP=`getopt -o hVv:b:H:R:n -l help,version,verbose:,p2b:,p2h:,p2R:,noclean \
        -n "$0" -- "$@"`
    else # original getopt is available (no long options, whitespace, sorting)
	TEMP=`getopt hVv:b:H:R:n "$@"`
    fi
    if [ $? -ne 0 ]; then
	echo "ERROR: "$(which getopt)" failed"
	getopt -T > /dev/null
	if [ $? -ne 4 ]; then
	    echo "did you use long options? they are not handled \
on your system, use -h for help"
	fi
	exit 2
    fi
    eval set -- "$TEMP"
    while [ $# -gt 0 ]; do
        case "$1" in
            -h | --help) help; exit 0; shift;;
            -V | --version) version; exit 0; shift;;
            -v | --verbose) verbose=$2; shift 2;;
            -b | --p2b) pathToBf=$2; shift 2;;
	    -H | --p2h) pathToHm=$2; shift 2;;
	    -R | --p2R) pathToRscript=$2; shift 2;;
	    -n | --noclean) clean=false; shift;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed, use -h for help"; exit 1;;
        esac
    done
    if [ ! -f "${pathToBf}" ]; then
	echo "ERROR: can't find path to 'eqtlbma_bf' -> '${pathToBf}'\n"
	help
	exit 1
    fi
    if [ ! -f "${pathToHm}" ]; then
	echo -e "ERROR: can't find path to 'eqtlbma_hm' -> '${pathToHm}'\n" 1>&2
	help
	exit 1
    fi
    if [ ! -f "${pathToRscript}" ]; then
	echo -e "ERROR: can't find path to 'functional_tests.R' -> '${pathToRscript}'\n" 1>&2
	help
	exit 1
    fi
}

#------------------------------------------------------------------------------

function simul_data_and_calc_exp_res () {
    if [ $verbose -gt "0" ]; then
	echo "simulate data and calculate expected results ..."
    fi
    ${pathToRscript} --verbose 1 --dir $(pwd) >& stdout_simul_exp
}

function calc_obs_res () {
    if [ $verbose -gt "0" ]; then
	echo "analyze data to get observed results ..."
    fi
    $pathToBf --geno list_genotypes.txt --scoord snp_coords.bed.gz \
	--exp list_phenotypes.txt --gcoord gene_coords.bed.gz --cis 5 \
	--out obs_bf --outss --outw --analys join --bfs all \
	--gridL grid_phi2_oma2_general.txt.gz \
	--gridS grid_phi2_oma2_with-configs.txt.gz \
	-v 1 >& stdout_bf
    $pathToHm --data obs_bf_l10abfs_raw.txt.gz --nsubgrp 3 --dim 7 --ngrid 10 \
    	--out obs_hm.txt.gz --getbf --getci >& stdout_hm
    zcat obs_hm.txt.gz | grep "#grid" | cut -f2 > obs_grid_probas.txt
    zcat obs_hm.txt.gz | grep "#config" | awk '{split($1,a,"."); print a[2]"\t"$2}' > obs_config_probas.txt
    zcat obs_hm.txt.gz | grep "#" | sed 's/#//g' > obs_file_for_ci.txt
    # eqtlbma_avg_bfs --in obs_bf_l10abfs_raw.txt.gz --gwts obs_grid_probas.txt --model configs --nsubgrp 3 --dim 7 --cwts obs_config_probas.txt --save bf+post --pi0 2.2517e-01 --post a+b+c+d --bestsnp 2 --alldim --out obs_avg_bfs.txt.gz -v 2
    
    # ~/src/eqtlbma/src/hm/eqtlbma_hm --data obs_bf_l10abfs_raw.txt.gz --nsubgrp 3 --model types --dim 4 --ngrid 10 --out obs_hm_types.txt.gz --configs "1|2|3" --getbf >& stdout_hm_types &
    # zcat obs_hm_types.txt.gz | grep "#type\|#subgroup" | sed 's/#//g' | cut -f1,2 > obs_type_subgroup_probas.txt
    # ~/src/eqtlbma/src/eqtlbma_avg_bfs --in obs_bf_l10abfs_raw.txt.gz --gwts obs_grid_probas.txt --model types --nsubgrp 3 --dim 4 --tswts obs_type_subgroup_probas.txt --save bf+post --pi0 2.1345e-01 --post a+b --alldim --out obs_avg_bfs_types.txt.gz
}

function comp_obs_vs_exp () {
    if [ $verbose -gt "0" ]; then
	echo "compare obs vs exp results ..."
    fi
    
    for i in {1..3}; do
    # nbDiffs=$(diff <(zcat obs_bf_sumstats_s${i}.txt.gz) <(zcat exp_bf_sumstats_s${i}.txt.gz) | wc -l)
    # if [ ! $nbDiffs -eq 0 ]; then
	if ! zcmp -s obs_bf_sumstats_s${i}.txt.gz exp_bf_sumstats_s${i}.txt.gz; then
	    echo "file 'obs_bf_sumstats_s${i}.txt.gz' has differences with exp"
	    exit 1
	fi
    done
    
    if ! zcmp -s obs_bf_l10abfs_raw.txt.gz exp_bf_l10abfs_raw.txt.gz; then
    	echo "file 'obs_bf_l10abfs_raw.txt.gz' has differences with exp"
	exit 1
    fi
    
    if ! zcmp -s obs_bf_l10abfs_avg-grids.txt.gz exp_bf_l10abfs_avg-grids.txt.gz; then
    	echo "file 'obs_bf_l10abfs_avg-grids.txt.gz' has differences with exp"
	exit 1
    fi
    
    if [ ! $(zcat obs_hm.txt.gz | grep "pi0" | cut -f2) == "2.2517e-01" ]; then
	echo "file 'obs_hm.txt.gz' has differences with exp"
	exit 1
    fi
    
    if [ $verbose -gt "0" ]; then
	echo "all tests passed successfully!"
    fi
}

#------------------------------------------------------------------------------

verbose=1
pathToBf=$bf_abspath
pathToHm=$hm_abspath
pathToRscript=$Rscript_abspath
clean=true
parseCmdLine "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@
    echo -e $msg
fi

cwd=$(pwd)

uniqId=$$ # process ID
testDir=tmp_test_${uniqId}
rm -rf ${testDir}
mkdir ${testDir}
cd ${testDir}
if [ $verbose -gt "0" ]; then echo "temp dir: "$(pwd); fi

simul_data_and_calc_exp_res

calc_obs_res

comp_obs_vs_exp

cd ${cwd}
if $clean; then rm -rf ${testDir}; fi

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
