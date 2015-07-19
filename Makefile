SAGE=sage

SAGEFILES =					\
	bug_examples.sage			\
	compendium_procedures.sage		\
	continuous_case.sage			\
	discontinuous_case.sage			\
	extreme_functions_in_literature.sage	\
	extreme_functions_sporadic.sage		\
	functions.sage				\
	simple_extremality_test.sage		\
	survey_examples.sage 			\
	extreme_functions_mlr_cpl3.sage		\
	quasi_periodic.sage

all:
	@echo "No need to 'make' anything. Just run it in Sage; see README.rst"

install:
	@echo "No need to install anything. Just run it in Sage; see README.rst"

check: check-encoding
	$(SAGE) -tp 4 $(SAGEFILES)

check-encoding:
	@if LC_ALL=C grep -v -n '^[ -}]*$$' $(SAGEFILES) ; then echo "Offending characters found."; exit 1; else echo "All Sage files are ASCII and have no tabs. Good."; exit 0; fi

tags: $(SAGEFILES)
	etags $(SAGEFILES)
