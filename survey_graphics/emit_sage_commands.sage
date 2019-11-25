## Generate the file sage-commands.tex,
## which is useful in papers about the code.

destdir=""

emitted_names = set()

def emit_tex_sage_command(name):
    if name not in emitted_names:
        sage_commands.write(r'\pgfkeyssetvalue{/sagefunc/' + name + r'}{\href{\githubsearchurl?q=\%22def+' + name.replace('\\', r'\\') + '(\\%22}{\\sage{' + name.replace(r'_', r'\underscore{}') + '}}}%)' + '\n')
        emitted_names.add(name)


with open(destdir + "sage-commands.tex", "w") as sage_commands:

    # extreme functions
    for f in dir(extreme_functions):
        if not f.startswith("_"):
            emit_tex_sage_command(f)

    # compendium procedures (from survey):
    for f in [ 'automorphism', 'multiplicative_homomorphism', 'projected_sequential_merge',
               'restrict_to_finite_group', 'interpolate_to_infinite_group', 'two_slope_fill_in' ]:
        emit_tex_sage_command(f)

    # survey examples:
    for f in ['generate_example_e_for_psi_n', 'chen_3_slope_not_extreme', 'psi_n_in_bccz_counterexample_construction', 'gomory_fractional', 'not_minimal_2', 'not_extreme_1', 'kzh_2q_example_1', 'zhou_two_sided_discontinuous_cannot_assume_any_continuity']:
        emit_tex_sage_command(f)


    # code mentioned in survey:
    for f in ['extremality_test', 'plot_2d_diagram']:
        emit_tex_sage_command(f)

    # code mentioned in param-abstract.tex:
    for f in ['nice_field_values', 'ParametricRealFieldElement', 'ParametricRealField']:
        emit_tex_sage_command(f)
