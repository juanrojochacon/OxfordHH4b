
echo " "
echo " *********************** "
echo " "

for plots in  Boosted_disc_noPU.pdf           m_HH_C2_res_noPU.pdf            pt_leadSJ_fj2_noPU.pdf D2_h0_bst_comp.pdf              m_h0_res_comp_back.pdf          pt_smallRjets_res_noPU.pdf EEC_C2_h0_C1_boost.pdf          m_hh_C2_bst_back_noPU.pdf       pt_subleadSJ_fj2_noPU.pdf Intermediate_disc_noPU.pdf      m_hh_C2_res_back_noPU.pdf     Resolved_disc_noPU.pdf          m_htot_bst_signal_PUnoTrim.pdf  res_wgthist_SKPU80.pdf           m_htot_res_signal_PUnoSK.pdf    res_wgthist_noPU.pdf bst_nnarch_noPU.pdf             nev2_SKPU80.pdf                 roc_SKPU80.pdf bst_wgthist_SKPU80.pdf          nev2_noPU.pdf                   roc_noPU.pdf bst_wgthist_noPU.pdf            pt_H0_C2_bst_comp.pdf           sb_SKPU80.pdf eta_smallRjets_res_noPU.pdf     pt_H0_C2_res_comp.pdf           sb_noPU.pdf hhFeyn.pdf                      pt_H0_bst_C1d_noPU.pdf          split12_h0_bst_comp_back.pdf        pt_H1_bst_C1d_noPU.pdf          split12_h1_C1_boost.pdf int_wgthist_SKPU80.pdf          pt_HH_C2_bst_noPU.pdf           ssb_SKPU80.pdf int_wgthist_noPU.pdf            pt_HH_C2_res_noPU.pdf           ssb_noPU.pdf m_H0_bst_C1d_noPU.pdf           pt_h0_C2_bst_back_noPU.pdf      tau21_h0_C1_boost.pdf m_H0_res_C1d_noPU.pdf           pt_h0_C2_res_back_noPU.pdf      tau21_h0_bst_comp.pdf m_HH_C2_bst_comp.pdf            pt_h0_bst_comp_back.pdf         tau21_h1_C1_boost.pdf m_HH_C2_bst_noPU.pdf            pt_h0_res_comp_back.pdf         tau21_h1_bst_comp_back.pdf m_HH_C2_res_comp.pdf            pt_leadSJ_fj1_bst_comp_back.pdf

do
    echo "---------------------------"
    echo "Grepping for plot = "$plots
    grep -i $plots ../*.tex
    

done
