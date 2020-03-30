# Simple local test suite for pfasta

test_check(){
    if [ "$1" != "$2" ]
    then
	echo "Test $3 failed"
	echo $2
    else
	echo "Test $3 passed"
    fi
}

get_random_string(){
    x=$(cat /dev/urandom | base64 | tr -dc '0-9a-zA-Z' | head -c10)
    echo "output_files/${x}"
}
    

####
echo "Running simple testing framework"
date
echo "Test setup..."
if [ -e output.fasta ]
then
    rm output.fasta
fi
rm output_files/*


outname='output_files/test.fa'

#test 0 defaut
filename='100_seqs_with_allvalid.fasta'
pfasta $filename > /dev/null 2>&1
v=$(cat output.fasta | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 0
rm output.fasta

# ................................................................................................
## test 1
filename='100_seqs_with_allvalid.fasta'

outname=$(get_random_string)
pfasta $filename --invalid-sequence ignore --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.1


outname=$(get_random_string)
pfasta $filename --invalid-sequence convert-all --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.2


outname=$(get_random_string)
pfasta $filename --invalid-sequence convert-res --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.3


outname=$(get_random_string)
pfasta $filename --invalid-sequence convert-all-ignore --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.4


outname=$(get_random_string)
pfasta $filename --invalid-sequence convert-res-ignore --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.5


outname=$(get_random_string)
pfasta $filename --invalid-sequence fail --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 1.6




# ................................................................................................
## test 2



## test 2
outname=$(get_random_string)
pfasta $filename --print-statistics --invalid-sequence ignore --random-subsample 100 -o ${outname}  > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 100 $v 2


## test 2 (
pfasta $filename --random-sub 50 -o $outname > /dev/null 2>&1
v=$(cat ${outname} | grep ">" |wc | awk {'print $1'}) 
test_check 50 $v 3

## test 3 (silent mode)
tmp=$(pfasta input_test.fasta --random-sub 100 -o $outname --silent)
if [ ! -z $tmp ]
then
    echo "Test 3 failed"
else
    echo "Test 3 passed"
fi
    
outname='output_files/test_general.fa'
pfasta input_test.fasta --random-sub 20 -o $outname > /dev/null 2>&1
v=$(cat $outname | grep ">" |wc | awk {'print $1'}) 
test_check 20 $v 4
