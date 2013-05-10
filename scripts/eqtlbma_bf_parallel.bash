#!/usr/bin/env bash

# Aim: used to launch eqtlbma_bf in parallel
# Author: Timothee Flutre
# Not copyrighted -- provided to the public domain

function help () {
    msg="\`${0##*/}' is used to launch eqtlbma_bf in parallel.\n"
    msg+="\n"
    msg+="Usage: ${0##*/} [OPTIONS] ...\n"
    msg+="\n"
    msg+="Options:\n"
    msg+="  -h, --help\tdisplay the help and exit\n"
    msg+="  -V, --version\toutput version information and exit\n"
    msg+="  -v, --verbose\tverbosity level (0/default=1/2/3)\n"
    msg+="      --p2b\tpath to the binary 'eqtlbma_bf'\n"
    msg+="      --geneD\tdirectory with lists of features to analyze (BED files)\n"
    msg+="\t\tfile names have to be like '<anything>_<batchId>.<anything>'\n"
    msg+="      --snpD\tdirectory with lists of SNPs to analyze (optional)\n"
    msg+="\t\tfile names have to be like '<anything>_<batchId>.<anything>'\n"
    msg+="\t\teach SNP file should correspond to a feature file, in the same order\n"
    msg+="      --seedF\tfile with seeds (as many as files in --geneD)\n"
    msg+="\t\toptional, default=list_seeds.txt.gz (should be gzipped)\n"
    msg+="      --task\ttask identifier (not for SGE, for SLURM only)\n"
    msg+="\n"
    msg+="Options from 'eqtlbma_bf' (run 'eqtlbma_bf -h | less' for more details):\n"
    msg+="      --geno\tfile with absolute paths to genotype files\n"
    msg+="\t\tdefault=list_genotypes.txt\n"
    msg+="      --scoord\tfile with the SNP coordinates\n"
    msg+="      --exp\tfile with absolute paths to expression level files\n"
    msg+="\t\tdefault=list_expressions.txt\n"
    msg+="      --anchor\tfeature boundary(ies) for the cis region\n"
    msg+="\t\tdefault=TSS\n"
    msg+="      --cis\tlength of half of the cis region (in bp)\n"
    msg+="\t\tdefault=100000\n"
    msg+="      --out\tprefix for the output files\n"
    msg+="\t\tdefault=out_eqtlbma\n"
    msg+="      --type\ttype of analysis to perform (sep/join)\n"
    msg+="\t\tdefault=sep\n"
    msg+="      --outss\twrite the output file with all summary statistics\n"
    msg+="      --outraw\twrite the output file with all raw ABFs\n"
    msg+="      --qnorm\tquantile-normalize the expression levels to a N(0,1)\n"
    msg+="      --maf\tminimum minor allele frequency\n"
    msg+="\t\tdefault=0\n"
    msg+="      --covar\tfile with absolute paths to covariate files\n"
    msg+="      --gridL\tfile with a 'large' grid for prior variances\n"
    msg+="\t\tdefault=grid_phi2_oma2_general.txt.gz\n"
    msg+="      --gridS\tfile with a 'small' grid for prior variances\n"
    msg+="\t\tdefault=grid_phi2_oma2_with-configs.txt.gz\n"
    msg+="      --bfs\twhich Bayes Factors to compute for the joint analysis\n"
    msg+="\t\tdefault=gen\n"
    msg+="      --error\tmodel for the errors (uvlr/mvlr/hybrid)\n"
    msg+="\t\tdefault=uvlr\n"
    msg+="      --fiterr\tparam used when estimating the variance of the errors\n"
    msg+="\t\tdefault=0.5\n"
    msg+="      --nperm\tnumber of permutations\n"
    msg+="\t\tdefault=0\n"
    msg+="      --trick\tapply trick to speed-up permutations\n"
    msg+="\t\tdefault=2\n"
    msg+="      --tricut\tcutoff for the trick\n"
    msg+="\t\tdefault=10\n"
    msg+="      --permsep\twhich permutation procedure for the separate analysis\n"
    msg+="\t\tdefault=0\n"
    msg+="      --pbf\twhich BF to use as the test statistic for the joint-analysis permutations\n"
    msg+="\t\tdefault=none\n"
    msg+="      --maxbf\tuse the maximum ABF over SNPs as test statistic for permutations\n"
    msg+="      --thread\tnumber of threads for the permutations\n"
    msg+="\t\tdefault=1\n"
    msg+="      --sbgrp\tidentifier of the subgroup to analyze\n"
    msg+="\n"
    msg+="Examples:\n"
    msg+="  ${0##*/} --p2b ~/bin/eqtlbma_bf --geneD ./lists_genes/\n"
    echo -e "$msg"
}

function version () {
    msg="${0##*/} 1.2\n"
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
    TEMP=`getopt -o hVv: -l help,version,verbose:,p2b:,geneD:,snpD:,seedF:,task:,geno:,scoord:,exp:,fcoord:,anchor:,cis:,out:,type:,outss,outraw,qnorm,maf:,covar:,gridL:,gridS:,bfs:,error:,fiterr:,nperm:,trick:,tricut:,permsep:,pbf:,maxbf,thread:,sbgrp: \
	-n "$0" -- "$@"`
    if [ $? != 0 ] ; then echo "ERROR: getopt failed" >&2 ; exit 1 ; fi
    eval set -- "$TEMP"
    while true; do
	case "$1" in
            -h|--help) help; exit 0; shift;;
            -V|--version) version; exit 0; shift;;
            -v|--verbose) verbose=$2; shift 2;;
	    --p2b) pathToBin=$2; shift 2;;
	    --geneD) geneDir=$2; shift 2;;
	    --snpD) snpDir=$2; shift 2;;
	    --seedF) seedFile=$2; shift 2;;
	    --task) task=$2; shift 2;;
	    --geno) geno=$2; shift 2;;
	    --scoord) scoord=$2; shift 2;;
	    --exp) exp=$2; shift 2;;
	    --fcoord) fcoord=$2; shift 2;;
	    --anchor) anchor=$2; shift 2;;
	    --cis) cis=$2; shift 2;;
	    --out) out=$2; shift 2;;
	    --type) type=$2; shift 2;;
	    --outss) outss=true; shift;;
	    --outraw) outraw=true; shift;;
	    --qnorm) qnorm=true; shift;;
	    --maf) maf=$2; shift 2;;
	    --covar) covar=$2; shift 2;;
	    --gridL) gridL=$2; shift 2;;
	    --gridS) gridS=$2; shift 2;;
	    --bfs) bfs=$2; shift 2;;
	    --error) error=$2; shift 2;;
	    --fiterr) fiterr=$2; shift 2;;
	    --nperm) nperm=$2; shift 2;;
	    --trick) trick=$2; shift 2;;
	    --tricut) trickCutoff=$2; shift 2;;
	    --pbf) pbf=$2; shift 2;;
	    --permsep) permsep=$2; shift 2;;
	    --maxbf) maxbf=true; shift;;
	    --thread) thread=$2; shift 2;;
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
    if [ -z "${geneDir}" ]; then
	echo -e "ERROR: missing compulsory option --geneDir\n"
	help
	exit 1
    fi
    if [ ! -d "${geneDir}" ]; then
	echo -e "ERROR: can't find feature directory '${geneDir}'\n"
	help
	exit 1
    fi
    if [ ! -z "${snpDir}" -a ! -d "${snpDir}" ]; then
	echo -e "ERROR: can't find SNP directory '${snpDir}'\n"
	help
	exit 1
    fi
    if [ ! -z "${seedFile}" -a ! -f "${seedFile}" ]; then
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
geneDir=""
snpDir=""
seedFile=""
task=""
geno="list_genotypes.txt"
scoord=""
exp="list_expressions.txt"
anchor="TSS"
cis=100000
out="out_eqtlbma"
type="sep"
outss=false
outraw=false
qnorm=false
maf=0
covar=""
gridL="grid_phi2_oma2_general.txt.gz"
gridS="grid_phi2_oma2_with-configs.txt.gz"
bfs="gen"
error="uvlr"
fiterr=0.5
nperm=0
trick=2
trickCutoff=10
permsep=0
pbf="none"
maxbf=false
thread=1
sbgrp=""
parseArgs "$@"

if [ $verbose -gt "0" ]; then
    startTime=$(timer)
    msg="START ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+="\ncmd-line: $0 "$@
    echo -e $msg
    uname -a
fi

# prepare the gene list + misc
gcoord=$(ls ${geneDir}/* | awk -v i=${task} 'NR==i{print;exit}');
split=$(echo ${gcoord} | awk '{n=split($0,a,"_"); split(a[n],b,"."); print b[1]}');
snp=""
if [ ! -z "${snpDir}" ]; then
    snp=$(ls ${snpDir}/* | awk -v i=${task} 'NR==i{print;exit}');
fi
if [ ! -z "${seedFile}" ]; then
    seed=$(zcat ${seedFile} | sed -n ${task}p)
fi

# build the command-line
cmd="${pathToBin}"
cmd+=" --geno ${geno}"
if [ ! -z "${scoord}" ]; then
    cmd+=" --scoord ${scoord}"
fi
cmd+=" --exp ${exp}"
cmd+=" --gcoord ${gcoord}"
cmd+=" --anchor ${anchor}"
cmd+=" --cis ${cis}"
cmd+=" --out ${out}_${split}"
cmd+=" --type ${type}"
if $outss; then
    cmd+=" --outss"
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
if [ "x${type}" == "xjoin" ]; then
    cmd+=" --gridL ${gridL}"
    cmd+=" --gridS ${gridS}"
    cmd+=" --bfs ${bfs}"
    cmd+=" --error ${error}"
    cmd+=" --fiterr ${fiterr}"
fi
if [ $nperm -gt "0" ]; then
    cmd+=" --nperm ${nperm}"
    if [ ! -z "${seed}$" ]; then
	cmd+=" --seed ${seed}"
    fi
    cmd+=" --trick ${trick}"
    cmd+=" --tricut ${trickCutoff}"
    if [ $permsep != "0" ]; then
	cmd+=" --permsep ${permsep}"
    fi
    if [ "x${pbf}" != "xnone" ]; then
	cmd+=" --pbf ${pbf}"
	if $maxbf; then
	    cmd+=" --maxbf"
	fi
    fi
    cmd+=" --thread ${thread}"
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
    echo "ERROR: eqtlbma_bf didn't finished successfully" >&2
    exit 1
fi

if [ $verbose -gt "0" ]; then
    msg="END ${0##*/} $(date +"%Y-%m-%d") $(date +"%H:%M:%S")"
    msg+=" ($(timer startTime))"
    echo $msg
fi
