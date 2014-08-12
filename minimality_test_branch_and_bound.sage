# Tree of closed intervals [a, b] to speed up Delta pi calculations

# Type1: We check x breakpoint, y breakpoint against x + y arbitrary point.
# So use discrete branching on x, y.

# Persistent tree: In one variable.
# Saves maximum and minimum information of the function over nodes of tree.

tree_threshold = 0  # Means allow branch down to singleton intervals

class IntervalNode:
    # An interval can be broken into three subintervals:
    # [a, v] (v, w) [w, b] 
    # where a, v, w, b are breakpoints 
    # and (v, w) is a subinterval that contains no breakpoint.
    # The indices a_index and b_index are the indices of the endpoints of the function
    # that are a and b
    
    def __init__(self, fn, a, b, a_index, b_index):
        self.fn = fn
        self.a = a
        self.b = b
        self.a_index = a_index
        self.b_index = b_index
        self.end_node = True
        self.min_max()

    def __repr__(self):
        return self._repr_with_level_(0)

    def _repr_with_level_(self, level, recurse = False):
        rep = "IntervalNode a: %s, b: %s, a_index: %s, b_index: %s, min: %s, max: %s" \
              % (self.a, self.b, self.a_index, self.b_index, self._min, self._max)
        if self.end_node:
            return "  " * level + "<" + rep + ">"
        else:
            rep += ", v: %s, w: %s " % (self.v, self.w)
            rep = "  " * level + "<" + rep 
            if recurse:
                rep += "\n" \
                       + self.left._repr_with_level_(level + 1, recurse=recurse) + "\n" \
                       + self.right._repr_with_level_(level + 1, recurse=recurse)
            rep += ">"
            return rep

    def print_tree(self):
        print self._repr_with_level_(0, True)

    def children(self):
        return self.left, self.right

    def min_max(self):
        mini = getattr(self, '_min', None)
        maxi = getattr(self, '_max', None)
        if mini is None:
            fn = self.fn
            end_points = fn.end_points()
            a_index = self.a_index
            b_index = self.b_index
            a = self.a
            b = self.b
            if b_index - a_index > tree_threshold:
                # Branch
                v_index = (b_index + a_index) // 2
                w_index = v_index + 1
                assert a_index <= v_index < w_index <= b_index
                v = self.v = end_points[v_index]
                w = self.w = end_points[w_index]
                self.left = IntervalNode(fn, a, v, a_index, v_index)
                self.right = IntervalNode(fn, w, b, w_index, b_index)
                left_min, left_max = self.left.min_max()
                right_min, right_max = self.right.min_max()
                mini = min(left_min, right_min)
                maxi = max(left_max, right_max)
                self.end_node = False
            else:
                mini = maxi = fn(end_points[self.b_index])
                for i in range(self.a_index, self.b_index):
                    value = fn(end_points[i])
                    if value < mini:
                        mini = value
                    elif value > maxi:
                        maxi = value
            self._min = mini
            self._max = maxi
        return mini, maxi
    
    def min_max_on_interval(self, c, d):
        ## later we want to bail out earlier based on a given range as well.
        #print "Enter node ", self
        c = max(self.a, c)
        d = min(self.b, d)

        if c > d:
            return None, None
        if c == self.a and d == self.b:
            return self._min, self._max

        minmax = [None, None]

        def update_minmax(small, big):
            #print "update: from ", minmax, "with ", small, big,  
            if small is not None and (minmax[0] is None or small < minmax[0]):
                minmax[0] = small
            if big is not None and (minmax[1] is None or big > minmax[1]):
                minmax[1] = big
            #print "to ", minmax 

        fn = self.fn
        if self.end_node:
            # Test at c and d, which may or may not be endpoints of the function.
            fnc = fn(c)
            update_minmax(fnc, fnc)
            fnd = fn(d)
            update_minmax(fnd, fnd)
            end_points = fn.end_points()
            for i in range(self.a_index, self.b_index):
                if c <= end_points[i] <= d:
                    fnx = fn(end_points[i])
                    update_minmax(fnx, fnx)
        else:    
            left, right = self.children()
            if c <= self.v:
                left_min, left_max = left.min_max_on_interval(c, d)
                update_minmax(left_min, left_max)
            if self.w <= d:
                right_min, right_max = right.min_max_on_interval(c, d)
                update_minmax(right_min, right_max)
            c = max(self.v, c)
            d = min(self.w, d)
            if c <= d:
                if self.v < c:
                    fnc = fn(c)
                    update_minmax(fnc, fnc)
                if c < d < self.w:
                    fnd = fn(d)
                    update_minmax(fnd, fnd)
        #print "min_max: ", minmax
        return minmax[0], minmax[1]

def generate_min_max_tree(fn):
    return IntervalNode(fn, 0, 1, 0, len(fn.end_points()) - 1)

# Branching tree: In two variables.


## class BranchNode:
    
##     def variable_intervals(self):
##         return self._variable_intervals

    

## def branch


