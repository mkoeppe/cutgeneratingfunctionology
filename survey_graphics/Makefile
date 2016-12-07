all:

check-graphics: all_graphics.pdf
	@echo "#######################################################"
	@echo "Now compare all_graphics.pdf with all_graphics_good.pdf"
	@echo "#######################################################"

all_graphics.pdf: \
		survey-graphics \
		extreme_notes_graphics/extreme_notes_graphics.pdf \
		mip_notes_graphics/mip_notes_graphics.pdf
	pdflatex -synctex=1 -src-specials -interaction=nonstopmode all_graphics

#survey_graphics/survey_graphics.pdf survey_graphics/compendium_graphics.pdf:
survey-graphics:
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_survey.sage")')

# this one is dependent on survey_graphics/sage-commmands.tex
extreme_notes_graphics/extreme_notes_graphics.pdf: survey-graphics
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_extreme_notes.sage")')

# this one is dependent on survey_graphics/sage-commmands.tex
mip_notes_graphics/mip_notes_graphics.pdf: survey-graphics
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_computer_based_search.sage")')

