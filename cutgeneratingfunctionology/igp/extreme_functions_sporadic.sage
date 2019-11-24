from six.moves import range
## "Sporadic" extreme functions (found by computer search, not part of a parametric family).
## (query-replace-regexp "^def \\(\\sw*\\)(\\([^)]*\\)):\\(.*\\)\n    r\"\"\"\n" "def \\1(\\2):\\3\n    r\"\"\"\n    .. PLOT::\n\n        from cutgeneratingfunctionology.igp import *\n        h = \\1()\n        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))\n        sphinx_plot(g)\n\n")


def hildebrand_5_slope_22_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_5_slope_22_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    One of Hildebrand's 5-slope functions.

    They held the world record regarding the number of slopes until
    functions with more slopes were discovered in 2014.

    From Hildebrand (2013, unpublished).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_22_1()
        sage: extremality_test(h, False)
        True
    """
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11, \
            1/2, 6/11, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 3/4, 1/4, 1/4, 1/2, 1/2, 3/4, 3/4, 1/4, 1, \
              1/2, 1/2, 3/4, 1/2, 3/4, 1/4, 1/2, 1/4, 1/2, 1/2, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def hildebrand_5_slope_24_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_5_slope_24_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    One of Hildebrand's 5-slope functions.

    They held the world record regarding the number of slopes until
    functions with more slopes were discovered in 2014.

    From Hildebrand (2013, unpublished).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_24_1()
        sage: extremality_test(h, False)
        True
    """
    bkpt = [0, 1/24, 1/12, 1/8, 1/6, 1/3, 3/8, 5/12, 11/24, 1/2, \
            13/24, 7/12, 5/8, 2/3, 17/24, 19/24, 5/6, 7/8, 11/12, 23/24, 1]
    values = [0, 1/2, 1/2, 3/4, 1/2, 1/2, 1/4, 1/2, 1/2, 1, \
              1/4, 1/2, 1/4, 1/4, 3/4, 1/4, 3/4, 3/4, 1/2, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def hildebrand_5_slope_28_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_5_slope_28_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    One of Hildebrand's 5-slope functions.

    They held the world record regarding the number of slopes until
    functions with more slopes were discovered in 2014.

    From Hildebrand (2013, unpublished).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_5_slope_28_1()
        sage: extremality_test(h, False)
        True
    """
    bkpt = [0, 1/28, 1/14, 3/28, 5/28, 3/14, 2/7, 9/28, 11/28, 3/7, 13/28, 1/2, \
            15/28, 4/7, 17/28, 9/14, 19/28, 5/7, 11/14, 23/28, 6/7, 25/28, 13/14, 27/28, 1]
    values = [0, 1/2, 1/2, 3/4, 3/4, 1/2, 1/2, 1/4, 1/4, 1/2, 1/2, 1, \
              1/4, 3/4, 1/2, 1/2, 1/4, 1/2, 1/2, 3/4, 1/2, 1/2, 1/4, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def hildebrand_2_sided_discont_1_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_2_sided_discont_1_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    The first known example of function that is discontinuous on both
    sides of the origin but is also extreme.

    Constructed by Robert Hildebrand (2013, unpublished).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_2_sided_discont_1_slope_1()
        sage: extremality_test(h, False)
        True
    """
    return FastPiecewise([singleton_piece(0, 0),
                          open_piece((0, 4/8), (1/8, 6/8)),
                          closed_piece((1/8, 2/8), (3/8, 6/8)),
                          open_piece((3/8, 2/8), (4/8, 4/8)),
                          singleton_piece(4/8, 1),
                          open_piece((4/8, 2/8), (5/8, 4/8)),
                          closed_piece((5/8, 2/8), (7/8, 6/8)),
                          open_piece((7/8, 4/8), (8/8, 6/8)),
                          singleton_piece(1, 0)])


def hildebrand_2_sided_discont_2_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_2_sided_discont_2_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    The second known example of function that is discontinuous on both
    sides of the origin but is also extreme.  This one has 2 slopes.

    Constructed by Robert Hildebrand (2013, unpublished).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: h = hildebrand_2_sided_discont_2_slope_1()
        sage: extremality_test(h, False)
        True
    """
    return FastPiecewise([singleton_piece(0, 0),
                          open_piece((0, 4/8), (1/8, 2/8)),
                          singleton_piece(1/8, 6/8),
                          open_piece((1/8, 2/8), (3/8, 6/8)),
                          singleton_piece(3/8, 2/8),
                          open_piece((3/8, 6/8), (4/8, 4/8)),
                          singleton_piece(4/8, 1),
                          open_piece((4/8, 4/8), (5/8, 6/8)),
                          singleton_piece(5/8, 2/8),
                          open_piece((5/8, 6/8), (7/8, 2/8)),
                          singleton_piece(7/8, 6/8),
                          open_piece((7/8, 2/8), (1, 4/8)),
                          singleton_piece(1, 0)])

def hildebrand_discont_3_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = hildebrand_discont_3_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    This is a very new discontinuous 3-slope function that is extreme.

    Constructed by Robert Hildebrand (2013, unpublished).

    A detailed extremality proof appears as an example in :cite:`hong-koeppe-zhou:software-paper`.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: psi = hildebrand_discont_3_slope_1()
        sage: extremality_test(psi, False)
        True

    In :cite:`koeppe-zhou:discontinuous-facets`, it is shown that this function (called `\psi`)
    is neither a weak facets nor a facet, by showing that the function
    ` \psi' = ``discontinuous_facets_paper_example_psi_prime()`` ` has a larger additivity
    domain `E` (sans limits)::

        sage: psi = hildebrand_discont_3_slope_1()
        sage: E_psi = set(generate_additive_faces_sans_limits(psi))
        sage: psi_prime = discontinuous_facets_paper_example_psi_prime(merge=False)
        sage: E_psi_prime = set(generate_additive_faces_sans_limits(psi_prime))
        sage: E_psi.issubset(E_psi_prime)
        True
        sage: sorted(E_psi_prime.difference(E_psi), key=lambda F: F.minimal_triple)
        [<Face ([0, 1/8], [0, 1/8], [1/8])>, <Face ([0, 1/8], [0, 1/8], [1/8, 1/4])>, ...]

    In fact, if one uses only the faces of `\psi` that are additive sans limits, then there is
    one covered component only; two intervals remain uncovered::

        sage: show_plots=False
        sage: show_plots=True      # not tested
        sage: sorted(generate_covered_components_strategically(psi, show_plots=show_plots,
        ....:                                                  additive_faces=E_psi))
        [[<Int(0, 1/8)>, <Int(3/8, 1/2)>],
         [<Int(1/8, 1/4)>, <Int(1/4, 3/8)>],
         [<Int(5/8, 3/4)>, <Int(3/4, 7/8)>]]

    """
    return FastPiecewise([right_open_piece((0, 0), (1/8, 6/8)),
                          closed_piece((1/8, 2/8), (3/8, 6/8)),
                          left_open_piece((3/8, 2/8), (4/8, 1)),
                          left_open_piece((4/8, 4/8), (5/8, 6/8)),
                          closed_piece((5/8, 6/8), (7/8, 2/8)), 
                          right_open_piece((7/8, 2/8), (1, 4/8)), 
                          singleton_piece(1, 0)])

def discontinuous_facets_paper_example_psi_prime(merge=True):
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = discontinuous_facets_paper_example_psi_prime()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    Example from the paper cite:`koeppe-zhou:discontinuous-facets`.
    See ``hildebrand_discont_3_slope_1``.

    If ``merge=False``, the function is set up on the same breakpoints
    as ``hildebrand_discont_3_slope_1`` for an easier comparison of the
    additive faces; if ``merge=True``, it is set up on a coarser set of breakpoints.

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: logging.disable(logging.INFO)             # Suppress output in automatic tests.
        sage: psi_prime = discontinuous_facets_paper_example_psi_prime()
        sage: extremality_test(psi_prime, False)
        True
    """
    return discontinuous_interpolation([0,   1/8, 3/8, 1/2, 5/8, 7/8],
                                       [0,   1/4, 3/4, 1,   3/4, 1/4],
                                       [0,   1/4, 3/4, 1/2, 3/4, 1/4],
                                       [1/2, 1/4, 3/4, 1,   3/4, 1/4], merge=merge)

def kzh_5_slope_fulldim_1(): #renamed from extreme_5slope_no_0d_1d_1
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A continuous 5-slope extreme function without any 0-d or 1-d
    maximal additive faces except for the symmetry reflection or x=0
    or y=0.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_1()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    bkpt = [0, 2/37, 3/37, 5/37, 7/37, 10/37, 11/37, 12/37, 13/37, 14/37, 15/37,
            18/37, 20/37, 22/37, 23/37, 25/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 1]
    values = [0, 59/90, 7/9, 64/135, 4/9, 73/90, 43/54, 29/45, 16/45, 11/54, 17/90,
              5/9, 71/135, 2/9, 31/90, 1, 19/45, 73/270, 23/90, 67/90, 197/270, 26/45, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_2(): #renamed from extreme_5slope_no_0d_1d_2
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A continuous 5-slope extreme function without any 0-d or 1-d
    maximal additive faces except for the symmetry reflection or x=0
    or y=0.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_2()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    bkpt = [0, 1/37, 2/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37,
            23/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 8/15, 4/13, 10/13, 47/78, 152/195, 239/390, 151/390, 43/195, 31/78,
              3/13, 9/13, 7/15, 1, 151/195, 539/780, 121/260, 149/390, 241/390, 139/260, 241/780, 44/195, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_3(): #renamed from extreme_5slope_no_0d_1d_3
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_3()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A continuous 5-slope extreme function without any 0-d or 1-d
    maximal additive faces except for the symmetry reflection or x=0
    or y=0.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_3()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37,
            22/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 8/15, 11/30, 49/60, 13/20, 77/100, 181/300, 119/300, 23/100, 7/20,
              11/60, 19/30, 7/15, 1, 119/150, 71/100, 151/300, 21/50, 29/50, 149/300, 29/100, 31/150, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_4(): #renamed from extreme_5slope_no_0d_1d_4
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_4()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37,
            22/37, 24/37, 25/37, 26/37, 27/37, 28/37, 29/37, 33/37, 34/37, 35/37, 36/37, 1]
    values = [0, 41/63, 61/126, 191/252, 149/252, 197/252, 155/252, 97/252, 55/252,
              103/252, 61/252, 65/126, 22/63, 1, 97/126, 173/252, 115/252, 47/126, 79/126, 137/252, 79/252, 29/126, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_5(): #renamed from extreme_5slope_no_0d_1d_5
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_5()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function without any 0-d or 1-d maximal additive faces
    except for the symmetry reflection or x=0 or y=0.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_5()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    bkpt = [0, 1/37, 3/37, 5/37, 7/37, 10/37, 12/37, 13/37, 15/37, 18/37, 20/37,
            22/37, 24/37, 25/37, 28/37, 29/37, 33/37, 34/37, 1]
    values = [0, 145/221, 6/13, 10/13, 127/221, 347/442, 261/442, 181/442, 95/442,
              94/221, 3/13, 7/13, 76/221, 1, 101/221, 159/442, 283/442, 120/221, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_1(): #renamed from extreme_5slope_no_transrefl or from fulldim_covers_5slope_q22_6()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.

    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_1()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 22; f = 10
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11,
            1/2, 6/11, 13/22, 7/11, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 23/32, 3/4, 7/16, 13/16, 3/16, 9/16, 1/4, 9/32, 1,
              11/32, 3/8, 3/4, 7/16, 9/16, 1/4, 5/8, 21/32, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_2(): # renamed from fulldim_covers_5slope_q22_1()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_2()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 22, f = 1
    bkpt = [0, 1/22, 3/22, 2/11, 5/22, 3/11, 4/11, 9/22, 5/11, 1/2, \
            6/11, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 10/11, 1]
    values = [0, 1, 1/4, 1/3, 31/48, 1/2, 2/3, 25/48, 5/6, \
              11/16, 5/16, 1/6, 23/48, 1/3, 1/2, 17/48, 2/3, 3/4, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_3(): # renamed from fulldim_covers_5slope_q22_2()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_3()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_3()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 22, f = 1
    bkpt = [0, 1/22, 3/22, 2/11, 5/22, 3/11, 4/11, 9/22, 5/11, 1/2, \
            6/11, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 10/11, 1]
    values = [0, 1, 31/113, 34/113, 81/113, 62/113, 68/113, 49/113, 96/113, 77/113, \
                36/113, 17/113, 64/113, 45/113, 51/113, 32/113, 79/113, 82/113, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_4(): # renamed from fulldim_covers_5slope_q22_3()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_4()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 22, f = 3
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 5/22, 3/11, 7/22, 9/22, 5/11, 1/2, \
            6/11, 13/22, 7/11, 15/22, 8/11, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 15/23, 8/23, 1, 26/69, 49/69, 50/69, 29/69, 31/69, 10/69, 11/23, \
                4/23, 19/23, 12/23, 59/69, 38/69, 40/69, 19/69, 20/69, 43/69, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_5(): # renamed from fulldim_covers_5slope_q22_4()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_5()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_5()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 22, f = 7
    bkpt = [0, 1/22, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 1/2, \
            13/22, 7/11, 15/22, 8/11, 9/11, 10/11, 21/22, 1]
    values = [0, 7/16, 5/8, 3/8, 9/16, 1, 13/32, 5/32, 11/16, \
               3/16, 9/32, 23/32, 13/16, 5/16, 27/32, 19/32, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_fulldim_covers_6(): # renamed from fulldim_covers_5slope_q22_5()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_fulldim_covers_6()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    5-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_5_slope_fulldim_covers_6()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []
    """
    # q = 22, f = 8
    bkpt = [0, 1/22, 1/11, 3/22, 5/22, 3/11, 7/22, 4/11, 9/22, \
            6/11, 13/22, 7/11, 8/11, 17/22, 9/11, 21/22, 1]
    values = [0, 7/12, 5/9, 2/9, 7/9, 4/9, 5/12, 1, 13/36, \
              5/18, 5/9, 2/9, 7/9, 4/9, 13/18, 23/36, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_fulldim_covers_1(): # renamed from fulldim_covers_6slope_q25_1()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    6-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_fulldim_covers_1()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 25; f = 8
    bkpt = [0, 1/25, 2/25, 6/25, 7/25, 8/25, 9/25, 2/5, 11/25, 12/25, \
            13/25, 3/5, 16/25, 17/25, 18/25, 4/5, 21/25, 22/25, 23/25, 24/25, 1]
    values = [0, 17/36, 1/4, 3/4, 19/36, 1, 37/144, 5/24, 121/288, 107/288, 7/12, \
                35/72, 19/72, 53/72, 37/72, 5/12, 181/288, 167/288, 19/24, 107/144, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_fulldim_covers_2(): # renamed from fulldim_covers_6slope_q26_1()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    6-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_fulldim_covers_2()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 26, f = 13
    bkpt = [0, 1/26, 3/26, 5/26, 3/13, 7/26, 4/13, 5/13, 6/13, 1/2, \
            7/13, 8/13, 9/13, 19/26, 10/13, 21/26, 23/26, 25/26, 1]
    values = [0, 3/7, 5/7, 1/7, 2/7, 5/7, 6/7, 2/7, 4/7, 1, \
              4/7, 2/7, 6/7, 5/7, 2/7, 1/7, 5/7, 3/7, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_fulldim_covers_3(): # renamed from fulldim_covers_6slope_q38_1()
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_3()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    6-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example has a similar 2d-diagram as that of kzh_6_slope_fulldim_covers_2()
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_fulldim_covers_3()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 38, f = 19
    bkpt = [0, 1/38, 5/38, 7/38, 9/38, 5/19, 6/19, 7/19, 9/19, 1/2, \
            10/19, 12/19, 13/19, 14/19, 29/38, 31/38, 33/38, 37/38, 1]
    values = [0, 7/17, 11/17, 3/17, 5/17, 12/17, 14/17, 6/17, 10/17, 1, \
            10/17, 6/17, 14/17, 12/17, 5/17, 3/17, 11/17, 7/17, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_fulldim_covers_4():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    6-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example has a similar 2d-diagram as that of kzh_6_slope_fulldim_covers_2()
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_fulldim_covers_4()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 27, f = 5
    bkpt = [0, 1/27, 2/27, 1/9, 4/27, 5/27, 2/9, 7/27, 1/3, 11/27, 13/27, 14/27, 2/3, 19/27, 7/9, 23/27, 25/27, 26/27, 1]
    values = [0, 47/127, 40/127, 87/127, 80/127, 1, 39/127, 86/127, 45/127, 58/127, 98/127, 209/254, 45/254, 29/127, 69/127, 82/127, 41/127, 88/127, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_fulldim_covers_5():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_5()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_fulldim_covers_5()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    6-slope extreme function whose extremality proof does not depend
    on lower-dimensional additive faces.  All intervals are directly covered.
    This is in contrast to ``hildebrand_5_slope_22_1`` etc., whose extremality proof
    requires to translate and reflect covered intervals.
    This example has a similar 2d-diagram as that of kzh_6_slope_fulldim_covers_2()
    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_fulldim_covers_5()
        sage: extremality_test(h)
        True
        sage: uncovered_intervals_from_covered_intervals(generate_directly_covered_intervals(h))
        []

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 27, f = 7
    bkpt = [0, 1/27, 2/27, 1/9, 4/27, 5/27, 2/9, 7/27, 8/27, 10/27, 11/27, 4/9, 13/27, 16/27, 2/3, 7/9, 22/27, 23/27, 8/9, 26/27, 1]
    values = [0, 43/76, 59/76, 3/4, 1/4, 17/76, 33/76, 1, 1/2, 17/38, 25/38, 3/4, 37/76, 29/38, 9/38, 39/76, 1/4, 13/38, 21/38, 1/2, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

#def kzh_5_slope_q22_f10_0():
#    return hildebrand_5_slope_22_1()

def kzh_5_slope_q22_f10_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_q22_f10_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    """
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11, 1/2, 6/11, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 29/48, 7/24, 5/24, 7/12, 5/12, 19/24, 17/24, 19/48, 1, 11/16, 29/48, 3/4, 2/3, 17/48, 31/48, 1/3, 1/4, 19/48, 5/16, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_q22_f10_2():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_q22_f10_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    """
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11, 1/2, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 21/22, 1]
    values = [0, 21/34, 10/17, 4/17, 9/17, 8/17, 13/17, 7/17, 13/34, 1, 11/34, 9/34, 19/34, 7/34, 27/34, 15/34, 25/34, 23/34, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_q22_f10_3():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_q22_f10_3()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    """
    bkpt = [0, 1/22, 1/11, 3/22, 7/22, 4/11, 9/22, 5/11, 1/2, 6/11, 13/22, 7/11, 15/22, 17/22, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 23/32, 3/4, 7/16, 9/16, 1/4, 9/32, 1, 11/32, 3/8, 3/4, 7/16, 13/16, 3/16, 9/16, 1/4, 5/8, 21/32, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_5_slope_q22_f10_4():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_q22_f10_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    """
    bkpt = [0, 1/22, 1/11, 3/22, 2/11, 3/11, 7/22, 4/11, 9/22, 5/11, 1/2, 6/11, 15/22, 17/22, 10/11, 21/22, 1]
    values = [0, 5/6, 3/4, 7/16, 17/48, 31/48, 9/16, 1/4, 1/6, 1, 11/24, 3/8, 13/16, 3/16, 5/8, 13/24, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

#def kzh_5_slope_q22_f10_5():
#    return kzh_5_slope_fulldim_covers_1()

def kzh_5_slope_q22_f2_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_5_slope_q22_f2_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    """
    bkpt = [0, 1/11, 3/22, 2/11, 5/22, 3/11, 7/22, 4/11, 5/11, 1/2, 13/22, 7/11, 8/11, 17/22, 9/11, 19/22, 10/11, 21/22, 1]
    values = [0, 1, 2/5, 8/15, 3/10, 13/30, 1/5, 7/10, 3/5, 11/15, 4/15, 2/5, 3/10, 4/5, 17/30, 7/10, 7/15, 3/5, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_7_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_7_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 7-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_7_slope_1()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 33; f = 11
    bkpt = [0, 1/33, 2/33, 1/11, 4/33, 7/33, 8/33, 3/11, 10/33, 1/3, 4/11, 13/33, 14/33, \
            5/11, 16/33, 19/33, 25/33, 28/33, 29/33, 10/11, 31/33, 32/33, 1]
    values = [0, 14/33, 2/11, 20/33, 4/11, 7/11, 13/33, 9/11, 19/33, 1, 14/33, 23/66, \
              20/33, 13/66, 5/11, 5/22, 17/22, 6/11, 53/66, 13/33, 43/66, 19/33, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_7_slope_2():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_7_slope_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 7-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_7_slope_2()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 34; f = 12
    bkpt = [0, 1/34, 1/17, 3/34, 2/17, 5/34, 7/34, 4/17, 9/34, 5/17, 11/34, 6/17, 13/34, 7/17, 15/34, 8/17, 1/2, \
            19/34, 10/17, 21/34, 25/34, 13/17, 27/34, 29/34, 15/17, 31/34, 16/17, 33/34, 1]
    values = [0, 59/232, 21/58, 143/232, 21/29, 159/232, 73/232, 8/29, 89/232, 37/58, 173/232, 1, 3/8, 73/116, 137/232, \
              111/232, 153/232, 101/232, 1/4, 4/29, 25/29, 3/4, 131/232, 79/232, 121/232, 95/232, 43/116, 5/8, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_7_slope_3():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_7_slope_3()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 7-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_7_slope_3()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 31; f = 10
    bkpt = [0, 2/31, 3/31, 4/31, 6/31, 7/31, 8/31, 10/31, 11/31, 13/31, 14/31, 15/31, 16/31, 17/31, \
        18/31, 19/31, 20/31, 21/31, 22/31, 23/31, 24/31, 25/31, 26/31, 27/31, 28/31, 30/31, 1]
    values = [0, 77/106, 19/53, 15/53, 38/53, 34/53, 29/106, 1, 18/53, 51/106, 55/212, 101/212, 29/53, 19/106, \
              1/4, 99/212, 57/106, 49/106, 113/212, 3/4, 87/106, 24/53, 111/212, 157/212, 55/106, 35/53, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_7_slope_4():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_7_slope_4()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 7-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_7_slope_4()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 29; f = 10
    bkpt = [0, 1/29, 2/29, 3/29, 4/29, 6/29, 7/29, 8/29, 9/29, 10/29, 11/29, 12/29, 15/29, \
            16/29, 17/29, 18/29, 19/29, 20/29, 21/29, 22/29, 23/29, 24/29, 27/29, 28/29, 1]
    values = [0, 65/99, 19/66, 7/33, 19/33, 14/33, 26/33, 47/66, 34/99, 1, 16/33, 26/99, 47/99, \
              25/99, 61/99, 107/198, 17/99, 82/99, 91/198, 38/99, 74/99, 52/99, 73/99, 17/33, 0]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def pattern0_sym_fn(l, sv):
    f = 18 * l + 11
    q = 2 * f
    s=[0, 1]
    for k in range(1, l + 1):
        s += [k, k + 1] * 3
    s += [l + 1] * 3
    for k in range(l, 0, -1):
        s += [k, k + 1, k]
    sk = s + [s[0]] + s[-1::-1]
    v_n = [0] * (f + 1)
    for k in range(f):
        v_n[k + 1] = v_n[k] + sv[sk[k]]
    v_d = v_n[-1]
    bkpt = [x / q for x in range(q + 1)]
    values = [y / v_d for y in v_n + v_n[-2::-1]]
    return piecewise_function_from_breakpoints_and_values(bkpt, values)

def kzh_6_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_6_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 6-slope extreme function.

    Its two-dimensional polyhedral complex includes special patterns.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_6_slope_1()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 58; f = 29
    l = 1
    sv = (11, 7, -5)
    return pattern0_sym_fn(l, sv)

def kzh_10_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_10_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 10-slope extreme function.

    Its two-dimensional polyhedral complex includes special patterns.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_10_slope_1()
        sage: extremality_test(h)
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 166; f = 83
    l = 4
    sv= (13, 12, 11, 9, -10, -12)
    return pattern0_sym_fn(l, sv)

def kzh_28_slope_1():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_28_slope_1()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 28-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    This was briefly the world record, until Basu, Conforti, Di Summa, and Paat gave a
    construction of a family of functions with an arbitrary number of slopes
    (see ``bcdsp_arbitrary_slope``).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_28_slope_1()
        sage: number_of_slopes(h)
        28
        sage: extremality_test(h) # long time
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 778; f = 389
    l = 21
    sv= (71755, 71715, 71715, 71655, 71655, 71595, 71595, 71595, 70995, 70995, 70611, 70611, \
         69267, 69267, -47799, -69939, -70835, -70835, -71555, -71555, -71635, -71685, -71725)
    return pattern0_sym_fn(l, sv)

def kzh_28_slope_2():
    r"""
    .. PLOT::

        from cutgeneratingfunctionology.igp import *
        h = kzh_28_slope_2()
        g = plot_with_colored_slopes(h, show_legend=False, aspect_ratio=0.125, figsize=(8, 1.5), thickness=2, **only_f_ticks_keywords(h))
        sphinx_plot(g)

    A 28-slope extreme function.

    This example was found by computer-based search
    described in Koeppe--Zhou [KZh2015a].

    This was briefly the world record, until Basu, Conforti, Di Summa, and Paat gave a
    construction of a family of functions with an arbitrary number of slopes
    (see ``bcdsp_arbitrary_slope``).

    EXAMPLES::

        sage: from cutgeneratingfunctionology.igp import *
        sage: h = kzh_28_slope_2()
        sage: number_of_slopes(h)
        28
        sage: extremality_test(h) # long time
        True

    Reference:

    [KZh2015a] M. Koeppe and Y. Zhou, New computer-based search strategies for
    extreme functions of the Gomory--Johnson infinite group problem, 2015,
    e-print http://arxiv.org/abs/1506.00017 [math.OC].
    """
    # q = 778; f = 389
    l = 21
    sv= (27455, 27447, 27447, 27435, 27435, 27423, 27423, 27423, 27303, 27303, 27303, 26919, \
         26343, 26343, -18171, -26631, -27271, -27271, -27415, -27415, -27431, -27441, -27449)
    return pattern0_sym_fn(l, sv)

###
### Add more functions here
###

