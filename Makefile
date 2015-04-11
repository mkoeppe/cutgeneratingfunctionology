SAGE=sage

SAGEFILES =					\
	bug_examples.sage			\
	compendium_procedures.sage		\
	continuous_case.sage			\
	discontinuous_case.sage			\
	extreme_functions_in_literature.sage	\
	functions.sage				\
	simple_extremality_test.sage		\
	survey_examples.sage


all:
	@echo "No need to 'make' anything. Just run it in Sage; see README.md"

install:
	@echo "No need to install anything. Just run it in Sage; see README.md"

check:
	$(SAGE) -tp 4 $(SAGEFILES)

tags: $(SAGEFILES)
	etags $(SAGEFILES)
