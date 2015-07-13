
rm -rf *.pdf

python2.7 ./plotHist_boosted_tau21_fj1.py
python2.7 ./plotHist_boosted_tau21_fj2.py
python2.7 ./plotHist_boosted_d12_fj2.py
python2.7 ./plotHist_boosted_C2_fj1.py

cp *.pdf ../plots/
