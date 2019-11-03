SAGE=sage

SAGEFILES =									\
	cutgeneratingfunctionology/igp/bug_examples.sage			\
	cutgeneratingfunctionology/igp/compendium_procedures.sage		\
	cutgeneratingfunctionology/igp/continuous_case.sage			\
	cutgeneratingfunctionology/igp/discontinuous_case.sage			\
	cutgeneratingfunctionology/igp/discrete_case.sage			\
	cutgeneratingfunctionology/igp/extreme_functions_in_literature.sage	\
	cutgeneratingfunctionology/igp/extreme_functions_sporadic.sage		\
	cutgeneratingfunctionology/igp/intervals.sage				\
	cutgeneratingfunctionology/igp/real_number_field.sage			\
	cutgeneratingfunctionology/igp/fast_linear.sage				\
	cutgeneratingfunctionology/igp/functions.sage				\
	cutgeneratingfunctionology/igp/parametric.sage				\
	cutgeneratingfunctionology/igp/semialgebraic_mathematica.sage		\
	cutgeneratingfunctionology/igp/simple_extremality_test.sage		\
	cutgeneratingfunctionology/igp/survey_examples.sage			\
	cutgeneratingfunctionology/igp/extreme_functions_mlr_cpl3.sage		\
	cutgeneratingfunctionology/igp/quasi_periodic.sage			\
	cutgeneratingfunctionology/igp/crazy_perturbation_examples.sage		\
	cutgeneratingfunctionology/igp/crazy_perturbation.sage			\
	cutgeneratingfunctionology/igp/kslope_ppl_mip.py			\
	cutgeneratingfunctionology/igp/vertex_enumeration.py			\
	cutgeneratingfunctionology/igp/kslope_pattern.sage			\
	cutgeneratingfunctionology/igp/2q_mip.sage				\
	cutgeneratingfunctionology/igp/kslope_mip.sage				\
	cutgeneratingfunctionology/igp/animation_2d_diagram.sage		\
	cutgeneratingfunctionology/igp/lifting_project.sage			\
	cutgeneratingfunctionology/igp/procedures/injective_2_slope_fill_in_proof.py

## Don't test; currently broken
# 	cutgeneratingfunctionology/igp/parametric_cpl.sage			\


# Dual feasible functions
SAGEFILES +=								    \
	cutgeneratingfunctionology/dff/dff_functions.sage		    \
	cutgeneratingfunctionology/dff/dff_test_plot.sage		    \
	cutgeneratingfunctionology/dff/discontinuous_dff.sage		    \
	cutgeneratingfunctionology/dff/computer_based_search_naive_dff.sage \
	cutgeneratingfunctionology/dff/gdff_linear_test.sage		    \
	cutgeneratingfunctionology/dff/Gomory_conversion.sage

# Multirow
SAGEFILES +=								\
	cutgeneratingfunctionology/multirow/piecewise_functions.sage	\
	cutgeneratingfunctionology/multirow/lifting_region.sage

SAGEFILES +=								\
	cutgeneratingfunctionology/spam/basic_semialgebraic.py		\
	cutgeneratingfunctionology/spam/semialgebraic_maple.py		\
	cutgeneratingfunctionology/spam/semialgebraic_mathematica.py	\
	cutgeneratingfunctionology/spam/semialgebraic_qepcad.py		\
	cutgeneratingfunctionology/spam/big_cells.py			\
	cutgeneratingfunctionology/spam/big_cells_impl.py		\
	cutgeneratingfunctionology/spam/real_set.py			\
	cutgeneratingfunctionology/spam/polyhedral_complex.py

all:
	@echo "No need to 'make' anything. Just run it in Sage; see README.rst"

install:
	@echo "No need to install anything. Just run it in Sage; see README.rst"

CHECK_PARALLEL=4

check: check-encoding
	PYTHONPATH=`pwd` $(SAGE) -tp $(CHECK_PARALLEL) --force_lib $(SAGEFILES)

check-long: check-encoding
	cp .check-long-timings.json .tmp_check-long-timings.json
	PYTHONPATH=`pwd` $(SAGE) -tp $(CHECK_PARALLEL) --force_lib --long --stats-path .tmp_check-long-timings.json $(SAGE_CHECK_FLAGS) $(SAGEFILES)
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

demo.ipynb: demo.rst
	$(SAGE) -rst2ipynb $< > $@

.PHONY: all build install test coverage sdist pip-install pip-uninstall pip-develop clean clean-doc doc doc-pdf
