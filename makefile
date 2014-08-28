paper.pdf: concept.tex conclusions.tex correction.tex discussion.tex functions.tex introduction.tex paper.tex paper.bib
	@rm -f paper.{aux,toc,lof,lot}
	@echo '---------------------------'
	@echo 'First pass pdflatex paper'
	@echo '---------------------------'
	pdflatex paper
	@echo '---------------------------'
	@echo 'BiBTeX pass bibtex paper'
	@echo '---------------------------'
	bibtex paper
	@echo '---------------------------'
	@echo 'Second pass pdflatex paper'
	@echo '---------------------------'
	pdflatex paper
	@echo '---------------------------'
	@echo 'Third pass pdflatex paper'
	@echo '---------------------------'
	pdflatex paper

clean:
	@rm -f paper.pdf paper.log paper.aux
	@rm -f *.bbl *.blg *.lof *.cut
	@rm -f *.lot *.out *.toc
