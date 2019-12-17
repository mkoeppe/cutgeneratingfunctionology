LATEXMK = latexmk
PDFLATEX = pdflatex

LATEX_OPTIONS = -interaction=nonstopmode

LATEXMK_OPTIONS = -pdf -g -f --quiet
# latexmk 1:4.41-1 on GitHub runner (2019-12-17) does not understand -Werror

check: check-bibtex check-biblatex

check-bibtex:
	$(LATEXMK) $(LATEXMK_OPTIONS) -pdflatex="$(PDFLATEX) $(LATEX_OPTIONS)" bibliography-bibtex 2>&1 | tee $@.log
	@if grep -i 'Warning\|Error' $@.log; then echo "BibTeX produces warnings"; exit 1; fi

check-biblatex:
	$(LATEXMK) $(LATEXMK_OPTIONS) -pdflatex="$(PDFLATEX) $(LATEX_OPTIONS)" bibliography-biblatex 2>&1 | tee $@.log
	@if grep -i 'Warning\|Error' $@.log; then echo "Biber produces warnings"; exit 1; fi
