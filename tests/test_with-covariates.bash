#!/usr/bin/env bash

set -o errexit -o pipefail

# Aim: launch a functional test for eqtlbma with covariates
# Author: Timothee Flutre
# Not copyrighted -- provided to the public domain

#------------------------------------------------------------------------------

function help () {
    msg="\`${0##*/}' launches a functional test for eqtlbma with covariates.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="      --p2e\tabsolute path to the 'eqtlbma' binary\n"
    msg+="      --p2R\tabsolute path to the 'functional_tests.R' script\n"
    msg+="      --noclean\tkeep temporary directory with all files\n"
    echo -e "$msg"
}

function version () {
    msg="${0##*/} 1.0\n"
    msg+="\n"
    msg+="Written by Timothee Flutre.\n"
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

function parseArgs () {
    TEMP=`getopt -o hVv: -l help,version,verbose:,p2e:,p2R:,noclean \
        -n "$0" -- "$@"`
    if [ $? != 0 ] ; then echo "ERROR: getopt failed" >&2 ; exit 1 ; fi
    eval set -- "$TEMP"
    while true; do
        case "$1" in
            -h|--help) help; exit 0; shift;;
            -V|--version) version; exit 0; shift;;
            -v|--verbose) verbose=$2; shift 2;;
            --p2e) pathToEqtlBma=$2; shift 2;;
	    --p2R) pathToRscript=$2; shift 2;;
	    --noclean) clean=false; shift;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed"; exit 1;;
        esac
    done
    if [[ ! -f $pathToEqtlBma ]]; then
	echo "ERROR: can't find path to 'eqtlbma' -> '${pathToEqtlBma}'"
	exit 1
    fi
    if [[ ! -f $pathToRscript ]]; then
	echo "ERROR: can't find path to 'functional_tests.R' -> '${pathToRscript}'"
	exit 1
    fi
}

#------------------------------------------------------------------------------

function simul_data_and_calc_exp_res () {
    if [ $verbose -gt "0" ]; then
	echo "simulate data and calculate expected results ..."
    fi
    R --no-restore --no-save --slave --vanilla --file=${pathToRscript} \
	--args --verbose 1 --dir $(pwd) --cvrt >& stdout_simul_exp
}

function calc_obs_res () {
    if [ $verbose -gt "0" ]; then
	echo "analyze data to get observed results ..."
    fi
    $pathToEqtlBma -g list_genotypes.txt --scoord snp_coords.bed.gz \
	-p list_phenotypes.txt --fcoord gene_coords.bed.gz --cis 5 \
	-o obs_eqtlbma --outss --outraw --step 3 --bfs all \
	--gridL grid_phi2_oma2_general.txt.gz \
	--gridS grid_phi2_oma2_with-configs.txt.gz \
	--covar list_covariates.txt --outcv \
	-v 1 >& stdout_eqtlbma
}

function comp_obs_vs_exp () {
    if [ $verbose -gt "0" ]; then
	echo "compare obs vs exp results ..."
    fi
    
    for i in {1..3}; do
    # nbDiffs=$(diff <(zcat obs_eqtlbma_sumstats_s${i}.txt.gz) <(zcat exp_eqtlbma_sumstats_s${i}.txt.gz) | wc -l)
    # if [ ! $nbDiffs -eq 0 ]; then
	if ! zcmp -s obs_eqtlbma_sumstats_s${i}.txt.gz exp_eqtlbma_sumstats_s${i}.txt.gz; then
	    echo "file 'obs_eqtlbma_sumstats_s${i}.txt.gz' has differences with exp"
		exit 1
	fi
    done
    
    if ! zcmp -s obs_eqtlbma_l10abfs_raw.txt.gz exp_eqtlbma_l10abfs_raw.txt.gz; then
    	echo "file 'obs_eqtlbma_l10abfs_raw.txt.gz' has differences with exp"
		exit 1
    fi
    
    if ! zcmp -s obs_eqtlbma_l10abfs_avg-grids.txt.gz exp_eqtlbma_l10abfs_avg-grids.txt.gz; then
    	echo "file 'obs_eqtlbma_l10abfs_avg-grids.txt.gz' has differences with exp"
		exit 1
    fi
    
    if [ $verbose -gt "0" ]; then
	echo "all tests passed successfully!"
    fi
}

#------------------------------------------------------------------------------

verbose=1
pathToEqtlBma=$eqtlbma_abspath
pathToRscript=$Rscript_abspath
clean=true
parseArgs "$@"

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
