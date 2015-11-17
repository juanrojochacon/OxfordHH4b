default:
	pdflatex HH4b.tex
	pdflatex HH4b.tex	
	bibtex HH4b
	pdflatex HH4b.tex	


clean:
	rm -rf *~ *.log *.pdf *.aux *.bbl *.toc *.dvi *.blg *.out
	rm -rf plotscripts/*~ plotscripts/*.pdf ../analysis/res/*~ ../analysis/res/*.dat paperscripts/*.pdf paperscripts/*~
