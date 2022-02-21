#!/bin/bash

# Test GRIN function with test multiplex network and gene set with restart = 0.7
echo "Testing GRIN installation."
echo "Running GRIN:"
GRIN_DIR=$PWD
Rscript $GRIN_DIR/R/GRIN.R -d $GRIN_DIR/test/suicide_weighted_Multiplex_0.5Delta.RData -g $GRIN_DIR/test/TestGenes.txt -r 0.7 -m Test_Install_User --tau 1,1,1,1,1,1,1,1,1,1 -o $GRIN_DIR/test/test_output

echo "Analyzing output:"
diff $GRIN_DIR/test/test_output/GRIN__Test_Install_User__Retained_Genes.txt $GRIN_DIR/test/test_output/GRIN__Test_Install_Reference__Retained_Genes.txt > output1.txt
diff $GRIN_DIR/test/test_output/GRIN__Test_Install_User__Removed_Genes.txt $GRIN_DIR/test/test_output/GRIN__Test_Install_Reference__Removed_Genes.txt > output2.txt
cat output1.txt output2.txt > output.txt

# If output from diff is not blank, then installation is not correct
correct=$(awk 'END{print NR}' output.txt)

if [ "$correct" == "0" ]; then
	echo "GRIN Installation successful!"
	echo "Removing test files..."
	rm output1.txt output2.txt output.txt test/test_output/GRIN__Test_Install_User__Retained_Genes.txt test/test_output/GRIN__Test_Install_User__Removed_Genes.txt
	echo "Installation test complete."
else
	echo "ERROR: improper output. Ensure that no files were moved following download and that conda installation is correct before reattempting GRIN."
	rm output1.txt output2.txt output.txt test/test_output/GRIN__Test_Install_User__Retained_Genes.txt test/test_output/GRIN__Test_Install_User__Removed_Genes.txt
fi