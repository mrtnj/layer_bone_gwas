## Helper function for running GWAS with GEMMA

function run_gwas () {
    NAME=$1
    PLINK_BASE=$2
    TRAIT_FILE=$3
    COVAR_FILE=$4

    if [ ! -d gwas/$NAME ]; then
	     mkdir gwas/$NAME
    fi

    cd gwas/$NAME

    if [ ! -f $NAME.bed ]; then
	    ln -s ../$PLINK_BASE.bim $NAME.bim
	    ln -s ../$PLINK_BASE.bed $NAME.bed
	    ln -s ../$PLINK_BASE.map $NAME.map
	    ln -s ../$TRAIT_FILE $NAME.fam
	    if [ COVAR_FILE != "NULL" ]; then
	        ln -s ../$COVAR_FILE covariates.txt
	    fi
    fi

    if [ $COVAR_FILE != "NULL" ]; then
	    $GEMMA_PATH/gemma -bfile $NAME \
			  -k ../output/${PLINK_BASE}_grm.sXX.txt \
			  -c covariates.txt \
			  -lmm 4 \
			  -o $NAME
    fi
    if [ $COVAR_FILE == "NULL" ]; then
	    $GEMMA_PATH/gemma -bfile $NAME \
			  -k ../output/${PLINK_BASE}_grm.sXX.txt \
			  -lmm 4 \
			  -o $NAME
    fi
    
    cd ../..
}