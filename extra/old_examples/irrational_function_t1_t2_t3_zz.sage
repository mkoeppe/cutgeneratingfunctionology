### Modification of ...general_number_field.
# length: 70028

del1 = 1/60
del2 = 2828/100000 #2000*(sqrt(2))/100000
del3 = (sqrt(3))/1000000

h = bhk_irrational(delta=(del1, del2, del3))

## def try_bfs(seed):
##     print "####### Seed ", seed, "######"
##     try:
##         stab_int, walk_dict = find_stability_interval_with_deterministic_walk_list(seed, generate_uncovered_intervals(h), generate_moves(h), h, max_num_it=500000)
##         print "Length: ", len(walk_dict)
##     except MaximumNumberOfIterationsReached:
##         print "Exceeded iteration limit"

## if minimality_test(h):
##     try_bfs(11/36)                             # "generic"
##     try_bfs(19/60)                             # contradition of signs



## z_stab_int, z_walk_list = find_stability_interval_with_deterministic_walk_list(19/60, generate_uncovered_intervals(h), generate_moves(h), h)

## crazy_x = sorted(z_walk_list.keys())[1]
   
