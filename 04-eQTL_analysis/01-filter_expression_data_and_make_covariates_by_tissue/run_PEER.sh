#!/bin/bash

# module load PEER

# expressiondatafile=$1
covariatedatafile=$1
tissue=$2
inputdir=$3
outputdir=$3

# expressiondatafilename=`echo $expressiondatafile | sed 's/.txt//g'`
covariatedatafilename=`echo $covariatedatafile | sed 's/.txt//g'`
# cat $inputdir/$tissue/$expressiondatafile | sed 's/\t/,/g' > $inputdir/$tissue/${expressiondatafilename}.csv
# peertool -f $inputdir/$tissue/${expressiondatafilename}.csv -n 15 --has_header -o $outputdir/$tissue/

cat $inputdir/$tissue/W.csv | sed 's/,/\t/g' > $inputdir/$tissue/temp1.txt
cat $inputdir/header_PEER_factor.txt $inputdir/$tissue/temp1.txt > $inputdir/$tissue/temp2.txt

awk '
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' $inputdir/$tissue/temp2.txt > $inputdir/$tissue/temp3.txt

cat $inputdir/$tissue/$covariatedatafile $inputdir/$tissue/temp3.txt > $inputdir/$tissue/temp4.txt
mv -f $inputdir/$tissue/temp4.txt $inputdir/$tissue/${covariatedatafilename}_v2.txt
rm -f $inputdir/$tissue/temp1.txt
rm -f $inputdir/$tissue/temp2.txt
rm -f $inputdir/$tissue/temp3.txt

echo -e "

USAGE: 

   peertool  [--sigma_off <float>] [--var_tol <float>] [--bound_tol
             <float>] [--e_pb <float>] [--e_pa <float>] [--a_pb <float>]
             [--a_pa <float>] [-i <int>] [-n <int>] [--prior <string>] [-c
             <string>] [--var_file <string>] -f <string> [-o <string>]
             [--has_header] [--add_mean] [--no_a_out] [--no_z_out]
             [--no_w_out] [--no_x_out] [--no_res_out] [--] [--version]
             [-h]


Where: 

   --sigma_off <float>
     Variance inactive component

   --var_tol <float>
     Variation tolerance

   --bound_tol <float>
     Bound tolerance

   --e_pb <float>
     Eps node prior parameter b

   --e_pa <float>
     Eps node prior parameter a

   --a_pb <float>
     Alpha node prior parameter b

   --a_pa <float>
     Alpha node prior parameter a

   -i <int>,  --n_iter <int>
     Number of iterations

   -n <int>,  --n_factors <int>
     Number of hidden factors

   --prior <string>
     Factor prior file

   -c <string>,  --cov_file <string>
     Covariate data file

   --var_file <string>
     Expression uncertainty (variance) data file

   -f <string>,  --file <string>
     (required)  Expression data file

   -o <string>,  --out_dir <string>
     Output directory

   --has_header
     Expression and covariates files have a header

   --add_mean
     Add a covariate to model mean effect

   --no_a_out
     No output of weight precision

   --no_z_out
     No output of posterior sparsity prior

   --no_w_out
     No output of estimated factor weights

   --no_x_out
     No output of estimated factors

   --no_res_out
     No output of residual values

   --,  --ignore_rest
     Ignores the rest of the labeled arguments following this flag.

   --version
     Displays version information and exits.

   -h,  --help
     Displays usage information and exits.


   Probabilistic estimation of expression residuals (PEER)

" > /dev/null

