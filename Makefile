SAGE=/Users/mkoeppe/s/sage/sage-5.11/sage

SAGEFILES =					\
	bug_examples.sage			\
	compendium_procedures.sage		\
	continuous_case.sage			\
	debug_examples.sage			\
	discontinuous_case.sage			\
	extreme_functions_in_literature.sage	\
	functions.sage				\
	survey_examples.sage

check:
	$(SAGE) -tp 4 $(SAGEFILES)

