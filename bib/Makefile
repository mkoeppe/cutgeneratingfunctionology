LATEXMK = latexmk
PDFLATEX = pdflatex
BIBER = biber

LATEX_OPTIONS = -interaction=nonstopmode

LATEXMK_OPTIONS = -pdf -g -f
# Omit -Werror because latexmk 1:4.41-1 on GitHub runner (2019-12-17) does not understand it
# Omit --quiet for better debuggability

check: check-encoding check-bibtex check-biblatex

BIBFILES = cutgeneratingfunctionology.bib

check-encoding:
	@if LC_ALL=C grep -H -v -n '^[ -~]*$$' $(BIBFILES) ; then echo "Offending characters found."; exit 1; else echo "bib files are ASCII and have no tabs and no ^Ms. Good."; exit 0; fi

fix-encoding:
	LC_ALL=en_US.UTF-8 sed -i.bak "s/–/--/g;s/’/\'/g;s/“/\`\`/g;s/”/\'\'/g;s/≡/\\\\equiv{}/g;s/≥/\\\\geq{}/g;s/{\\\\textparagraph}/ /g;" $(BIBFILES)

check-bibtex:
	$(LATEXMK) $(LATEXMK_OPTIONS) -pdflatex="$(PDFLATEX) $(LATEX_OPTIONS)" bibliography-bibtex 2>&1 | tee $@.log
	@if grep -i 'Warning\|Error' $@.log; then echo "BibTeX produces warnings"; exit 1; fi

check-biblatex:
	$(LATEXMK) $(LATEXMK_OPTIONS) -pdflatex="$(PDFLATEX) $(LATEX_OPTIONS)" bibliography-biblatex
	$(BIBER) bibliography-biblatex 2>&1 | tee $@.log
	@if grep -i 'Warning\|Error' $@.log; then echo "Biber produces warnings"; exit 1; fi
