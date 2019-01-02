### This example is a rational approximation of `irrational_function_t1_t2_t3_general_numberfield.sage'
### The breadth first search ends with a reachable orbit of 600000 elements. 
### Then finding epsilon takes forever.

del1 = 1/60
del2 = 21213 / 1000000 # approximates 3*(sqrt(2))/200
del3 = 1732 /  1000000 # approximates (sqrt(3))/1000

load("extreme_functions_in_literature.sage")
h = bhk_irrational(delta=(del1, del2, del3))

