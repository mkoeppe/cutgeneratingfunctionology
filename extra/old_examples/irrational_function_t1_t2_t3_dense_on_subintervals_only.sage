## This one has three dense subintervals in each of the flat intervals.
## but is not extreme.

t1, t2, t3 = 9/400, 1/4000*sqrt(2) + 9/200, 19/400
h = bhk_irrational(delta=[t1, t2-t1, t3-t2])
