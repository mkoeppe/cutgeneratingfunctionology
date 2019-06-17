from pyomo.environ import *
import numpy as np
import tensorflow as tf
from tensorflow import keras

def nn_model(trained_nn,user_model,linked_variable_index,bounds_type="lp"):
    user_model.nn_block=generate_block(trained_nn,bounds_type=bounds_type)
    user_model.linked_variable_index=linked_variable_index
    def link_variable(m,i):
        return m.nn_block.x[i]==m.x[i]
    user_model.link=Constraints(user_model.linked_variable_index,rule=link_variable)
    return user_model
    
def generate_block(trained_nn,bounds_type="lp"):
    dims,weights,bias=get_dims_weights_bias(trained_nn)
    initial_lb,initial_ub=[0]*dims[0],[1]*dims[0]
    if bounds_type=="1d":
        bounds=compute_1dim_bounds(dims,initial_lb,initial_ub,weights,bias)
    elif bounds_type=="lp":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=True)
    elif bounds_type=="ip":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=False)
    else:
        raise NameError('Wrong bound type')
    lbs=bounds[0]
    ubs=bounds[1]
    num_layer = len(dims)
    b = Block(concrete=True)

    varindex={}
    for i in range(num_layer):
        varindex[i]=list(range(dims[i]))
    b.vi=Set(initialize=list((key,j) for key in varindex.keys() for j in varindex[key]))

    # variables
    b.x=Var(b.vi,domain=Reals)
    b.a=Var(b.vi,domain=Reals)
    b.z=Var(b.vi,domain=Binary)

    # constraints
    b.cons=ConstraintList()
    for i in range(1,num_layer):
        if i == num_layer-1:
            for j in range(dims[i]):
                b.cons.add(b.x[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*b.x[(i-1,k)] for k in range(dims[i-1])))
        else:
            for j in range(dims[i]):
                b.cons.add(0<=b.x[(i,j)])
                b.cons.add(b.a[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*b.x[(i-1,k)] for k in range(dims[i-1])))
                b.cons.add(b.a[(i,j)]<=b.x[(i,j)])
                lb,ub=compute_one_node_bounds(lbs[i-1],ubs[i-1],weights[i-1][j],bias[i-1][j])
                b.cons.add(b.x[(i,j)]<=b.z[(i,j)]*ub)
                b.cons.add(b.x[(i,j)]<=b.a[(i,j)]-lb*(1-b.z[(i,j)]))
    return b

def nn_model(model,obj_coefs,input_bounds,sense=maximize,bounds_type="lp",relaxed=True):
    """
    Return a model which computes the maximum(minimum) of c*y subject to y=f(x) and L<=x<=U.
    Function f is given by the trained model, c=obj_coefs, and L,U are given by input_bounds.
    The bounds parameter defines how to compute the bigM used in the formulation.
    The relaxed parameter defines the type of the formulation (MIP or LP).
    """
    dims,weights,bias=get_dims_weights_bias(model)
    initial_lb,initial_ub=input_bounds[0],input_bounds[1]
    if bounds_type=="1d":
        bounds=compute_1dim_bounds(dims,initial_lb,initial_ub,weights,bias)
    elif bounds_type=="lp":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=True)
    elif bounds_type=="ip":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=False)
    else:
        raise NameError('Wrong bound type')
    m=nn_bigm_formulation(dims,weights,bias,bounds,obj_coefs,sense=maximize,relaxed=relaxed)
    return m

def generate_toy_trained_model():
    """
    Generate a trained neural network for testing.
    """
    fashion_mnist = keras.datasets.fashion_mnist
    (train_images, train_labels), (test_images, test_labels) = fashion_mnist.load_data()
    train_images = train_images / 255.0
    test_images = test_images / 255.0
    model = keras.Sequential([keras.layers.Flatten(input_shape=(28, 28)),keras.layers.Dense(32, activation=tf.nn.relu),keras.layers.Dense(16, activation=tf.nn.relu),keras.layers.Dense(10, activation=tf.nn.softmax)])
    model.compile(optimizer=tf.train.AdamOptimizer(),loss='sparse_categorical_crossentropy',metrics=['accuracy'])
    model.fit(train_images, train_labels, epochs=5)
    return model

def get_dims_weights_bias(model):
    """
    Return the dimensions, weights and bias of a trained model obtained in tensorflow.
    """
    weights=[]
    bias=[]
    n=int(len(model.weights)/2)
    for i in range(n):
        weights.append(model.get_weights()[2*i].transpose())
        bias.append(model.get_weights()[2*i+1])
    dims=[len(weights[0][0])]
    for i in range(len(bias)):
        dims.append(len(bias[i]))
    return dims,weights,bias

def nn_image(dims,weights,bias,input,input_layer,output_layer):
    """
    Given a trained nn (using only relu), input in the input_layer, compute the output in the output_layer. If input_layer=0 and output_layer equals the total number of layers, then this function is the prediction function (for logits).
    """
    cur=np.asarray(input)
    for i in range(output_layer-input_layer):
        W=np.asarray(weights[input_layer+i])
        b=np.asarray(bias[input_layer+i])
        if input_layer+i==len(dims)-2:
            cur=np.matmul(W,cur)+b
        else:
            cur=np.maximum(0,np.matmul(W,cur)+b)
    return cur

def nn_predict(dims,weights,bias,datapoint):
    """
    Return the prediction/class of the given datapoint.
    """
    logits=nn_image(dims,weights,bias,datapoint,0,len(dims)-1)
    return np.argmax(logits)

def compute_1dim_bounds(dims,initial_lb,initial_ub,weights,bias):
    """
    Compute the (loose) bounds of every node by 1-dimensional update layer by layer.
    """
    num_layer=len(dims)
    lb=[initial_lb]
    ub=[initial_ub]
    for i in range(1,num_layer):
        prev_lb=lb[-1]
        prev_ub=ub[-1]
        cur_lb=[]
        cur_ub=[]
        for k in range(dims[i]):
            l,u=compute_one_node_bounds(prev_lb,prev_ub,weights[i-1][k],bias[i-1][k])
            if i !=num_layer-1:
                l=max(0,l)
                u=max(0,u)
            cur_lb.append(l)
            cur_ub.append(u)
        lb.append(cur_lb)
        ub.append(cur_ub)
    return lb,ub

def compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=True):
    """
    Compute the bounds of every node by solving LPs of MIPs.
    """
    num_layer = len(dims)
    lb=[initial_lb]
    ub=[initial_ub]
    for i in range(1,num_layer):
        cur_lb=[]
        cur_ub=[]
        for k in range(dims[i]):
            l=nn_one_node_bigm(dims,weights,bias,initial_lb,initial_ub,i,k,solver='cplex',sense=minimize,relaxed=relaxed)
            u=nn_one_node_bigm(dims,weights,bias,initial_lb,initial_ub,i,k,solver='cplex',sense=maximize,relaxed=relaxed)
            cur_lb.append(l)
            cur_ub.append(u)
        lb.append(cur_lb)
        ub.append(cur_ub)
    return lb,ub

def compute_one_node_bounds(prev_lb,prev_ub,w,b):
    """
    Compute the lower and upper bound of a node from previous bounds of the previous layer and the affine transformation.
    """
    lb=0
    ub=0
    for i in range(len(w)):
        wi=w[i]
        if wi>0:
            lb+=wi*prev_lb[i]
            ub+=wi*prev_ub[i]
        else:
            lb+=wi*prev_ub[i]
            ub+=wi*prev_lb[i]
    lb+=b
    ub+=b
    return lb,ub

def nn_one_node_bigm(dims,weights,bias,initial_lb,initial_ub,layer_index,node_index,solver='cplex',sense=minimize,relaxed=False,solve=True,max_nodes=1000000):
    """
    Compute the lower or upper bound of one node indexed by (layer_index,node_index) by solving an LP or MIP.
    """
    num_layer = len(dims)
    model = ConcreteModel()
    varindex={}
    for i in range(layer_index+1):
        varindex[i]=list(range(dims[i]))
    model.vi=Set(initialize=list((key,j) for key in varindex.keys() for j in varindex[key]))

    # variables
    model.x=Var(model.vi,domain=Reals)
    model.a=Var(model.vi,domain=Reals)
    if relaxed:
        model.z=Var(model.vi,bounds=(0,1),domain=NonNegativeReals)
    else:
        model.z=Var(model.vi,domain=Binary)

    # objective function
    model.OBJ=Objective(expr=model.x[(layer_index,node_index)],sense=sense)

    # constraints
    model.cons=ConstraintList()
    for i in range(dims[0]):
        model.cons.add(model.x[(0,i)]>=initial_lb[i])
        model.cons.add(model.x[(0,i)]<=initial_ub[i])
    prev_lb=initial_lb
    prev_ub=initial_ub
    # only need first layer_index number of layers
    for i in range(1,layer_index+1):
        cur_lb=[]
        cur_ub=[]
        if i==num_layer-1:
            for j in range(dims[i]):
                model.cons.add(model.x[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*model.x[(i-1,k)] for k in range(dims[i-1])))
        else:
            for j in range(dims[i]):
                model.cons.add(0<=model.x[(i,j)])
                model.cons.add(model.a[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*model.x[(i-1,k)] for k in range(dims[i-1])))
                model.cons.add(model.a[(i,j)]<=model.x[(i,j)])
                lb,ub=compute_one_node_bounds(prev_lb,prev_ub,weights[i-1][j],bias[i-1][j])
                model.cons.add(model.x[(i,j)]<=model.z[(i,j)]*ub)
                model.cons.add(model.x[(i,j)]<=model.a[(i,j)]-lb*(1-model.z[(i,j)]))
                cur_lb.append(max(0,lb))
                cur_ub.append(max(0,ub))
            prev_lb=cur_lb
            prev_ub=cur_ub

    if not solve:
        return model
    # use cplex to solve for now
    if solver=='cplex':
        solver=SolverFactory('cplex',executable='/Applications/CPLEX_Studio128/cplex/bin/x86-64_osx/cplex')
        solver.options['mip limits nodes']=max_nodes
        solver.options['mip strategy bbinterval']=1
        solver.solve(model)
    return model.OBJ()

def nn_bigm_formulation(dims,weights,bias,bounds,obj_coefs,sense=minimize,relaxed=False):
    """
    Return a (relaxed) model of nn given weights, bias and precomputed bounds of every node. 
    """
    lbs=bounds[0]
    ubs=bounds[1]
    num_layer = len(dims)
    model = ConcreteModel()
    varindex={}
    for i in range(num_layer):
        varindex[i]=list(range(dims[i]))
    model.vi=Set(initialize=list((key,j) for key in varindex.keys() for j in varindex[key]))

    # variables
    model.x=Var(model.vi,domain=Reals)
    model.a=Var(model.vi,domain=Reals)
    if relaxed:
        model.z=Var(model.vi,bounds=(0,1),domain=NonNegativeReals)
    else:
        model.z=Var(model.vi,domain=Binary)
    
    # objective function
    model.OBJ=Objective(expr=sum(obj_coefs[i]*model.x[(num_layer-1,i)] for i in range(dims[num_layer-1])),sense=sense)
    
    # constraints
    model.cons=ConstraintList()
    for j in range(dims[0]):
        model.cons.add(model.x[(0,j)]>=lbs[0][j])
        model.cons.add(model.x[(0,j)]<=ubs[0][j])
    for i in range(1,num_layer):
        if i == num_layer-1:
            for j in range(dims[i]):
                model.cons.add(model.x[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*model.x[(i-1,k)] for k in range(dims[i-1])))
        else:
            for j in range(dims[i]):
                model.cons.add(0<=model.x[(i,j)])
                model.cons.add(model.a[(i,j)]==bias[i-1][j]+sum(weights[i-1][j][k]*model.x[(i-1,k)] for k in range(dims[i-1])))
                model.cons.add(model.a[(i,j)]<=model.x[(i,j)])
                lb,ub=compute_one_node_bounds(lbs[i-1],ubs[i-1],weights[i-1][j],bias[i-1][j])
                model.cons.add(model.x[(i,j)]<=model.z[(i,j)]*ub)
                model.cons.add(model.x[(i,j)]<=model.a[(i,j)]-lb*(1-model.z[(i,j)]))
    return model

def generate_epsilon_neighborhood(datapoint,epsilon):
    """
    Generate bounds of the epsilon neighborhood cented at datapoint. Default: input space is [0,1]^n0.
    """
    d=len(datapoint)
    lb=[]
    ub=[]
    for v in datapoint:
        lb.append(max(0,v-epsilon))
        ub.append(min(1,v+epsilon))
    return lb,ub

def generalization_error_model(model,datapoint,epsilon,comparing_class,bounds_type="lp",relaxed=False):
    """
    Generate a model which computes the maximum of (logit of comparing_class - logit of true_class) in the neighborhood of datapoint.
    True_class is the label of the datapoint. Ideally for well-traied network, the maximum should be negative if epsilon is small.
    """
    dims,weights,bias=get_dims_weights_bias(model)
    true_class=nn_predict(dims,weights,bias,datapoint)
    if true_class==comparing_class:
        return
    obj_coefs=[0]*dims[-1]
    obj_coefs[true_class]=-1
    obj_coefs[comparing_class]=1
    initial_lb,initial_ub=generate_epsilon_neighborhood(datapoint,epsilon)
    if bounds_type=="1d":
        bounds=compute_1dim_bounds(dims,initial_lb,initial_ub,weights,bias)
    elif bounds_type=="lp":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=True)
    elif bounds_type=="ip":
        bounds=compute_lp_bounds(dims,initial_lb,initial_ub,weights,bias,relaxed=False)
    else:
        raise NameError('Wrong bound type')
    m=nn_bigm_formulation(dims,weights,bias,bounds,obj_coefs,sense=maximize,relaxed=relaxed)
    return m




