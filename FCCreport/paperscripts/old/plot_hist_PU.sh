
# First of all some cleaning
rm -rf *.pdf

# Fig 9
# m_htot_res_signal_PUnoSK.pdf
./plotHist_m_htot_res_PUnoSK.py

# Fig 10
# m_H0_bst_comp
# m_H1_bst_comp
./plotHist_m_H0_bst_comp.py
./plotHist_m_H1_bst_comp.py
# Additional comparison plots for the mH distribution
./plotHist_m_H0_res_comp.py
./plotHist_m_H1_res_comp.py
./plotHist_m_H0_int_comp.py
./plotHist_m_H1_int_comp.py


# Fig 11
# pt_H0_C2_res_comp
# pt_H0_C2_bst_comp
# m_HH_C2_res_comp
# m_HH_C2_bst_comp
./plotHist_ptH0_res_comp.py
./plotHist_ptH0_bst_comp.py
./plotHist_mHH_res_comp.py
./plotHist_mHH_bst_comp.py
# Additional plots for comparison of kinematic variables in the noPU and PU80SK cases
./plotHist_pt_leadSJ_fj1_bst_comp.py
./plotHist_pt_subleadSJ_fj1_bst_comp.py

# Fig 12
# D2_h0_bst_comp
# tau21_h0_bst_comp
./plotHist_D2_h0_bst_comp.py
./plotHist_tau21_h0_bst_comp.py
# Additional comparison plots for other substructure variables
./plotHist_split12_h0_bst_comp.py
./plotHist_C2_h0_bst_comp.py


# Fig 13
# pt_h0_bst_comp_back
# tau21_h1_bst_comp_back
# split12_h0_bst_comp_back
# m_h0_bst_comp_back
./plotHist_pt_h0_bst_comp_back.py
./plotHist_tau12_h1_bst_comp_back.py
./plotHist_split12_h0_bst_comp_back.py
./plotHist_m_hh_bst_comp_back.py
# Additional comparison plots for signal vs background in the case of PU80SK
./plotHist_C2_h0_bst_comp_back.py
./plotHist_D2_h0_bst_comp_back.py
./plotHist_pt_leadSJ_fj1_bst_comp_back.py
./plotHist_pt_leadSJ_fj2_bst_comp_back.py
./plotHist_m_h0_bst_comp_back.py
./plotHist_m_h1_bst_comp_back.py
./plotHist_m_hh_res_comp_back.py


# Fig 14
# pt_h0_res_comp_back.pdf
# m_h0_res_comp_back.pdf
./plotHist_pt_h0_res_comp_back.py
./plotHist_m_h0_res_comp_back.py


# Now the plot the ANN weights
./plotNetwork_bst_PU.py
./plotNetwork_res_PU.py
./plotNetwork_int_PU.py

