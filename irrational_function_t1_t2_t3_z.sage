### Modification of ...general_number_field.

### Interesting: For this function, the orbit of the "generic" seed 11/36 (length 12)
### is much smaller than the orbit of a "contradiction of signs" seed 19/60 (length 90260)

h = bhk_irrational(delta = (1/60,
                            1981*(sqrt(2))/100000,    # 3*(sqrt(2))/200
                            (sqrt(2))/1000000))

# def try_bfs(seed):
#     print "####### Seed ", seed, "######"
#     try:
#         stab_int, walk_dict = find_stability_interval_with_deterministic_walk_list(seed, generate_uncovered_intervals(h), generate_moves(h), h, max_num_it=500000)
#         print "Length: ", len(walk_dict)
#     except MaximumNumberOfIterationsReached:
#         print "Exceeded iteration limit"

# if minimality_test(h):
#     try_bfs(11/36)                             # "generic"
#     try_bfs(19/60)                             # contradition of signs



## z_stab_int, z_walk_list = find_stability_interval_with_deterministic_walk_list(19/60, generate_uncovered_intervals(h), generate_moves(h), h)

## crazy_x = sorted(z_walk_list.keys())[1]
   
