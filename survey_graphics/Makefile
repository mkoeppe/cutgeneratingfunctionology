all:

check-graphics: all_graphics.pdf
	@echo "#######################################################"
	@echo "Now compare all_graphics.pdf with all_graphics_good.pdf"
	@echo "#######################################################"

all_graphics.pdf: \
		survey-graphics \
		extreme_notes_graphics/extreme_notes_graphics.pdf \
		mip_notes_graphics/mip_notes_graphics.pdf \
		software_paper_graphics/software_paper_graphics.pdf \
		crazy_perturbation_graphics/crazy_perturbation_graphics.pdf
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

software_paper_graphics/software_paper_graphics.pdf:
	rm -f software_paper_graphics/animation_2d_diagram_disc*.pdf software_paper_graphics/gj2slope*.pdf
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_software_paper.sage")')

.PHONY: software_paper_graphics/software_paper_graphics.pdf

crazy_perturbation_graphics/crazy_perturbation_graphics.pdf:
	(cd crazy_perturbation_graphics && rm -f crazy_perturbation_graphics.pdf has_crazy_perturbation.png proof_uniform_V_case12.png proof_uniform_V_case3.png proof_uniform_W_case12.png proof_uniform_W_case3a.png proof_uniform_W_case3b.png proof_uniform_cont_case_1.png proof_uniform_cont_case_2.png proof_uniform_cont_case_3.png zhou_two_sided_discontinuous_cannot_assume_any_continuity-covered_intervals.png zhou_two_sided_discontinuous_cannot_assume_any_continuity-perturbation-1.png)
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_crazy_perturbation.sage")')

.PHONY: crazy_perturbation_graphics/crazy_perturbation_graphics.pdf

algo-paper-graphics:
	(cd .. && sage -c 'sys.path = [""] + sys.path; load("survey_graphics/graphics_for_algo_paper.sage")')
