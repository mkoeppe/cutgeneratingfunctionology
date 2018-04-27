SAGE=sage

SAGEFILES =					\
	bug_examples.sage			\
	compendium_procedures.sage		\
	continuous_case.sage			\
	discontinuous_case.sage			\
	extreme_functions_in_literature.sage	\
	extreme_functions_sporadic.sage		\
	intervals.sage				\
	real_number_field.sage			\
	fast_linear.sage			\
	functions.sage				\
	parametric.sage				\
	simple_extremality_test.sage		\
	survey_examples.sage			\
	extreme_functions_mlr_cpl3.sage		\
	quasi_periodic.sage			\
	crazy_perturbation_examples.sage	\
	crazy_perturbation.sage			\
	kslope_ppl_mip.py			\
	vertex_enumeration.py			\
	kslope_pattern.sage			\
	2q_mip.sage				\
	kslope_mip.sage

# Dual feasible functions
SAGEFILES +=					\
	dff_functions.sage			\
	dff_test_plot.sage			\
	discontinuous_dff.sage			\
	computer_based_search_naive_dff.sage

all:
	@echo "No need to 'make' anything. Just run it in Sage; see README.rst"

install:
	@echo "No need to install anything. Just run it in Sage; see README.rst"

check: check-encoding
	$(SAGE) -tp 4 $(SAGEFILES)

check-long: check-encoding
	cp .check-long-timings.json .tmp_check-long-timings.json
	$(SAGE) -tp 4 --long --stats-path .tmp_check-long-timings.json $(SAGE_CHECK_FLAGS) $(SAGEFILES)
	rm .tmp_check-long-timings.json

check-encoding:
	@if LC_ALL=C grep -v -n '^[ -~]*$$' $(SAGEFILES) ; then echo "Offending characters found."; exit 1; else echo "All Sage files are ASCII and have no tabs. Good."; exit 0; fi

## Checking graphics takes long and requires manual inspection, so it's not part of 'make check'.
check-graphics:
	(cd survey_graphics && $(MAKE) check-graphics)

tags: $(SAGEFILES)
	etags $(SAGEFILES)

doc:
	cd docs && $(SAGE) -sh -c "make html"

doc-pdf:
	cd docs && $(SAGE) -sh -c "make latexpdf"

clean: clean-doc
	rm -f .tmp_check-long-timings.json
	rm -rf build dist *.egg-info
	rm -rf $(PACKAGE)/*.c

clean-doc:
	cd docs && $(SAGE) -sh -c "make clean"

.PHONY: all build install test coverage sdist pip-install pip-uninstall pip-develop clean clean-doc doc doc-pdf
