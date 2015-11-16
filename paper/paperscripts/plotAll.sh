
rm -rf *~ *.pdf

for script in Hist_mHH_bst_noPU  MVAdat_noPU Network_bst_noPU
do

    echo "script being executed = " plot$script.py

    python2.7 ./plot$script.py

done
