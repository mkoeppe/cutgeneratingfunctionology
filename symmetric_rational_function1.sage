slope = [10/3,0,10/3,0,10/3,-10/3,0,-10/3,0,-10/3]
interval_length = [1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10,1/10]

bkpt = []
bkpt.append(0)
for i in range(len(interval_length)):
    bkpt.append(bkpt[i]+interval_length[i])
    
bkpt2 = []
for i in range(len(bkpt)-1):
    bkpt2.append(bkpt[i])
for i in range(len(bkpt)):
    bkpt2.append(bkpt[i]+1)

function_values = [0]
  
for i in range(1,len(bkpt)-1):
    function_values.append(function_values[i-1] + slope[i - 1] * (bkpt[i] - bkpt[i-1]))

pieces = [[(bkpt[i],bkpt[i+1]),lambda x,i=i: function_values[i]+slope[i]*(x - bkpt[i])] for i in range(len(bkpt)-1)] 
h = Piecewise(pieces)