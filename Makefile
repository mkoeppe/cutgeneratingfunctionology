SAGE=sage

SAGEFILES =									\
	cutgeneratingfunctionology/igp/bug_examples.sage			\
	cutgeneratingfunctionology/igp/compendium_procedures.sage		\
	cutgeneratingfunctionology/igp/continuous_case.sage			\
	cutgeneratingfunctionology/igp/discontinuous_case.sage			\
	cutgeneratingfunctionology/igp/discrete_case.sage			\
	cutgeneratingfunctionology/igp/extreme_functions_in_literature.sage	\
	cutgeneratingfunctionology/igp/extreme_functions_sporadic.sage		\
	cutgeneratingfunctionology/igp/real_number_field.sage			\
	cutgeneratingfunctionology/igp/functions.sage				\
	cutgeneratingfunctionology/igp/parametric.sage				\
	cutgeneratingfunctionology/igp/parametric_cpl.sage			\
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
	cutgeneratingfunctionology/igp/faster_subadditivity_test.sage

# Separate modules under igp
SAGEFILES +=										\
	cutgeneratingfunctionology/igp/class_call.py					\
	cutgeneratingfunctionology/igp/fast_linear.py					\
	cutgeneratingfunctionology/igp/fast_piecewise.py				\
	cutgeneratingfunctionology/igp/intervals.py					\
	cutgeneratingfunctionology/igp/move_semigroup.py				\
	cutgeneratingfunctionology/igp/parametric_family.py				\
	cutgeneratingfunctionology/igp/procedures/injective_2_slope_fill_in_proof.py

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

# SPAM
SAGEFILES +=									\
	cutgeneratingfunctionology/spam/basic_semialgebraic.py			\
	cutgeneratingfunctionology/spam/basic_semialgebraic_formal_closure.py	\
	cutgeneratingfunctionology/spam/basic_semialgebraic_linear_system.py	\
	cutgeneratingfunctionology/spam/basic_semialgebraic_intersection.py	\
	cutgeneratingfunctionology/spam/basic_semialgebraic_local.py		\
	cutgeneratingfunctionology/spam/semialgebraic_predicate.py		\
	cutgeneratingfunctionology/spam/semialgebraic_qepcad.py			\
	cutgeneratingfunctionology/spam/big_cells.py				\
	cutgeneratingfunctionology/spam/big_cells_impl.py			\
	cutgeneratingfunctionology/spam/parametric_real_field_element.py	\
	cutgeneratingfunctionology/spam/real_set.py				\
	cutgeneratingfunctionology/spam/polyhedral_complex.py			\
	cutgeneratingfunctionology/spam/examples/relu.py


ifneq ($(shell command -v math > /dev/null 2>&1),)
SAGEFILES += \
	cutgeneratingfunctionology/spam/semialgebraic_mathematica.py
endif

ifneq ($(shell command -v maple > /dev/null 2>&1),)
SAGEFILES += \
	cutgeneratingfunctionology/spam/semialgebraic_maple.py
endif

CHECK_PARALLEL=4

-include Makefile.conf

all:
	@echo "No need to 'make' anything. Just run it in Sage; see README.rst"

install:
	$(SAGE) -pip install .

check: check-encoding
	PYTHONPATH=`pwd` $(SAGE) -tp $(CHECK_PARALLEL) --force_lib --warn-long 10 $(SAGE_CHECK_FLAGS) $(SAGEFILES)

check-long: check-encoding
	cp .check-long-timings.json .tmp_check-long-timings.json
	PYTHONPATH=`pwd` $(SAGE) -tp $(CHECK_PARALLEL) --force_lib --long --warn-long 300 --stats-path .tmp_check-long-timings.json $(SAGE_CHECK_FLAGS) $(SAGEFILES)
	rm .tmp_check-long-timings.json

.PHONY: fixdoctests
fixdoctests: $(SAGEFILES:%=%.fixdoctests)

%.fixdoctests: %
	PYTHONPATH=`pwd` $(SAGE) -fixdoctests --long $<
check-bib:
	$(MAKE) -C open-optimization-bibliography/test check

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
