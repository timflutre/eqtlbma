#!/usr/bin/env bash

# Aim: used to launch eqtlbma in parallel
# Author: Timothee Flutre
# Not copyrighted -- provided to the public domain

function help () {
    msg="\`${0##*/}' is used to launch eqtlbma in parallel.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="      --p2b\tpath to the binary 'eqtlbma'\n"
    msg+="      --ftrD\tdirectory with lists of features to analyze (BED files)\n"
    msg+="\t\tfile names have to be like '<anything>_<batchId>.<anything>'\n"
    msg+="      --snpD\tdirectory with lists of SNPs to analyze\n"
    msg+="\t\tfile names have to be like '<anything>_<batchId>.<anything>'\n"
    msg+="\t\teach SNP file should correspond to a feature file, in the same order\n"
    msg+="      --seedF\tfile with seeds (as many as files in --ftrD)\n"
    msg+="\t\toptional, default=list_seeds.txt.gz (should be gzipped)\n"
    msg+="      --task\ttask identifier (not for SGE, for SLURM only)\n"
    msg+="\n"
    msg+="Options from 'eqtlbma' (run 'eqtlbma -h | less' for more details):\n"
    msg+="      --geno\tfile with absolute paths to genotype files\n"
    msg+="\t\tdefault=list_genotypes.txt\n"
    msg+="      --scoord\tfile with the SNP coordinates\n"
    msg+="      --pheno\tfile with absolute paths to phenotype files\n"
    msg+="\t\tdefault=list_phenotypes.txt\n"
    msg+="      --anchor\tfeature boundary(ies) for the cis region\n"
    msg+="\t\tdefault=FSS\n"
    msg+="      --cis\tlength of half of the cis region (in bp)\n"
    msg+="\t\tdefault=100000\n"
    msg+="      --out\tprefix for the output files\n"
    msg+="\t\tdefault=out_eqtlbma\n"
    msg+="      --step\tstep of the analysis to perform (1/2/3/4/5)\n"
    msg+="\t\tdefault=1\n"
    msg+="      --outss\twrite the output file with all summary statistics\n"
    msg+="      --outcv\twrite the results for the other covariates when --outss\n"
    msg+="      --outraw\twrite the output file with all raw ABFs\n"
    msg+="      --qnorm\tquantile-normalize the phenotypes to a N(0,1)\n"
    msg+="      --maf\tminimum minor allele frequency\n"
    msg+="\t\tdefault=0\n"
    msg+="      --covar\tfile with absolute paths to covariate files\n"
    msg+="      --gridL\tfile with a 'large' grid for prior variances\n"
    msg+="\t\tdefault=grid_phi2_oma2_general.txt.gz\n"
    msg+="      --gridS\tfile with a 'small' grid for prior variances\n"
    msg+="\t\tdefault=grid_phi2_oma2_with-configs.txt.gz\n"
    msg+="      --bfs\twhich Bayes Factors to compute for the joint analysis\n"
    msg+="\t\tdefault=all\n"
    msg+="      --mvlr\tuse the multivariate version of the ABF\n"
    msg+="      --fitsig\tparam used when estimating the variance of the errors with --mvlr\n"
    msg+="\t\tdefault=0.0\n"
    msg+="      --nperm\tnumber of permutations\n"
    msg+="\t\tdefault=10000\n"
    msg+="      --trick\tapply trick to speed-up permutations\n"
    msg+="\t\tdefault=2\n"
    msg+="      --tricut\tcutoff for the trick\n"
    msg+="\t\tdefault=10\n"
    msg+="      --permsep\twhich permutation procedure for the separate analysis\n"
    msg+="\t\tdefault=1\n"
    msg+="      --pbf\twhich BF to use as the test statistic for the joint-analysis permutations\n"
    msg+="\t\tdefault=all\n"
    msg+="      --maxbf\tuse the maximum ABF over SNPs as test statistic for permutations\n"
    msg+="      --sbgrp\tidentifier of the subgroup to analyze\n"
    msg+="\n"
    msg+="Examples:\n"
    msg+="  ${0##*/} --p2b ~/bin/eqtlbma --ftrD ./lists_genes/\n"
    echo -e "$msg"
}

function version () {
    msg="${0##*/} 1.1\n"
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
    TEMP=`getopt -o hVv: -l help,version,verbose:,p2b:,ftrD:,snpD:,seedF:,task:,geno:,scoord:,pheno:,fcoord:,anchor:,cis:,out:,step:,outss,outcv,outraw,qnorm,maf:,covar:,gridL:,gridS:,bfs:,mvlr,fitsig:,nperm:,trick:,tricut:,permsep:,pbf:,maxbf,sbgrp: \
	-n "$0" -- "$@"`
    if [ $? != 0 ] ; then echo "ERROR: getopt failed" >&2 ; exit 1 ; fi
    eval set -- "$TEMP"
    while true; do
	case "$1" in
            -h|--help) help; exit 0; shift;;
            -V|--version) version; exit 0; shift;;
            -v|--verbose) verbose=$2; shift 2;;
	    --p2b) pathToBin=$2; shift 2;;
	    --ftrD) ftrDir=$2; shift 2;;
	    --snpD) snpDir=$2; shift 2;;
	    --seedF) seedFile=$2; shift 2;;
	    --task) task=$2; shift 2;;
	    --geno) geno=$2; shift 2;;
	    --scoord) scoord=$2; shift 2;;
	    --pheno) pheno=$2; shift 2;;
	    --fcoord) fcoord=$2; shift 2;;
	    --anchor) anchor=$2; shift 2;;
	    --cis) cis=$2; shift 2;;
	    --out) out=$2; shift 2;;
	    --step) step=$2; shift 2;;
	    --outss) outss=true; shift;;
	    --outcv) outss=true; shift;;
	    --outraw) outraw=true; shift;;
	    --qnorm) qnorm=true; shift;;
	    --maf) maf=$2; shift 2;;
	    --covar) covar=$2; shift 2;;
	    --gridL) gridL=$2; shift 2;;
	    --gridS) gridS=$2; shift 2;;
	    --bfs) bfs=$2; shift 2;;
	    --mvlr) mvlr=true; shift;;
	    --fitsig) fitsig=$2; shift 2;;
	    --nperm) nperm=$2; shift 2;;
	    --trick) trick=$2; shift 2;;
	    --tricut) trickCutoff=$2; shift 2;;
	    --pbf) pbf=$2; shift 2;;
	    --permsep) permsep=$2; shift 2;;
	    --maxbf) maxbf=true; shift;;
	    --sbgrp) sbgrp=$2; shift 2;;
            --) shift; break;;
            *) echo "ERROR: options parsing failed"; exit 1;;
	esac
    done
    if [ -z "${pathToBin}" ]; then
	echo -e "ERROR: missing compulsory option --p2b\n"
	help
	exit 1
    fi
    if [ ! -f "${pathToBin}" ]; then
	echo -e "ERROR: can't find binary '${pathToBin}'\n"
	help
	exit 1
    fi
    if [ ! -x "${pathToBin}" ]; then
	echo -e "ERROR: can't execute '${pathToBin}'\n"
	help
	exit 1
    fi
    if [ -z "${ftrDir}" ]; then
	echo -e "ERROR: missing compulsory option --ftrDir\n"
	help
	exit 1
    fi
    if [ ! -d "${ftrDir}" ]; then
	echo -e "ERROR: can't find feature directory '${ftrDir}'\n"
	help
	exit 1
    fi
    if [ ! -f "${seedFile}" ]; then
	echo -e "ERROR: can't find seed file '${seedFile}'\n"
	help
	exit 1
    fi
    if [ -z "${task}" -a ! -z "${SGE_TASK_ID}" ]; then
	task=${SGE_TASK_ID}
    fi
}

verbose=1
pathToBin=""
ftrDir=""
snpDir=""
seedFile=""
task=""
geno="list_genotypes.txt"
scoord=""
pheno="list_phenotypes.txt"
anchor="FSS"
cis=100000
out="out_eqtlbma"
step=1
outss=false
outcv=false
outraw=false
qnorm=false
maf=0
covar=""
gridL="grid_phi2_oma2_general.txt.gz"
gridS="grid_phi2_oma2_with-configs.txt.gz"
bfs="all"
mvlr=false
fitsig=0
nperm=10000
trick=2
trickCutoff=10
permsep=1
pbf="all"
maxbf=false
sbgrp=""
parseArgs "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@
    echo -e $msg
    uname -a
fi

# prepare the feature list + misc
fcoord=$(ls ${ftrDir}/* | awk -v i=${task} 'NR==i{print;exit}');
split=$(echo ${fcoord} | awk '{n=split($0,a,"_"); split(a[n],b,"."); print b[1]}');
snp=""
if [ ! -z "${snpDir}" ]; then
    snp=$(ls ${snpDir}/* | awk -v i=${task} 'NR==i{print;exit}');
fi
if [ ! -z "${seedFile}" ]; then
    seed=$(zcat ${seedFile} | sed -n ${task}p)
fi

# build the command-line
cmd="${pathToBin}"
cmd+=" -g ${geno}"
if [ ! -z "${scoord}" ]; then
    cmd+=" --scoord ${scoord}"
fi
cmd+=" -p ${pheno}"
cmd+=" --fcoord ${fcoord}"
cmd+=" --anchor ${anchor}"
cmd+=" --cis ${cis}"
cmd+=" --out ${out}_${split}"
cmd+=" --step ${step}"
if $outss; then
    cmd+=" --outss"
fi
if $outcv; then
    cmd+=" --outcv"
fi
if $outraw; then
    cmd+=" --outraw"
fi
if $qnorm; then
    cmd+=" --qnorm"
fi
cmd+=" --maf ${maf}"
if [ ! -z "${covar}" ]; then
    cmd+=" --covar ${covar}"
fi
if [ "x${step}" != "x1" -a "x${step}" != "x2" ]; then
    cmd+=" --gridL ${gridL}"
    cmd+=" --gridS ${gridS}"
    cmd+=" --bfs ${bfs}"
fi
if $mvlr; then
    cmd+=" --mvlr"
    cmd+=" --fitsig ${fitsig}"
fi
if [ "x${step}" != "x1" -a "x${step}" != "x3" ]; then
    cmd+=" --nperm ${nperm}"
    if [ ! -z "${seed}$" ]; then
	cmd+=" --seed ${seed}"
    fi
    cmd+=" --trick ${trick}"
    cmd+=" --tricut ${trickCutoff}"
    cmd+=" --permsep ${permsep}"
    if [ "x${step}" == "x4" -o "x${step}" == "x5" ]; then
	cmd+=" --pbf ${pbf}"
	if $maxbf; then
	    cmd+=" --maxbf"
	fi
    fi
fi
if [ ! -z "${snpDir}" ]; then
    cmd+=" --snp ${snp}"
fi
if [ ! -z "${sbgrp}" ]; then
    cmd+=" --sbgrp ${sbgrp}"
fi
cmd+=" -v ${verbose}"
#echo $cmd; exit

# run the program
eval $cmd
if [ $? != 0 ]; then
    echo "ERROR: eqtlbma didn't finished successfully" >&2
    exit 1
fi

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
