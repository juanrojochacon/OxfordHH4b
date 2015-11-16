
rm -rf *~ *.pdf

for script in Hist_mHH_bst_noPU \
		  Hist_m_H0_bst_comp  \
		  Hist_m_H0_res_comp  \
		  Network_bst_noPU  \
		  Network_bst_PU \
		  MVAdat_noPU MVAdat_PU \
		  MVAdat_comp_bst \
		  MVAdat_comp_int \
		  MVAdat_comp_res
do

    echo "script being executed = " plot$script.py

    python2.7 ./plot$script.py

done
