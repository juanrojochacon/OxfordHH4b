
rm -rf *~ *.pdf

for script in MVAdat_noPU Network_bst_noPU
do

    echo "script being executed = " plot$$script.py

    python2.7 ./plot$script.py

done
