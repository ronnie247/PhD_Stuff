import pandas as pd
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import time
from matplotlib.ticker import PercentFormatter
import warnings
warnings.filterwarnings("ignore")
from matplotlib import inline

#SIMPLE RANDOM WALK
def simpl_rnd_walk(n):
    """
    Generates a simple random walk of length n
    
    Args:
        n (int): the length of the walk
    Returns:
        (x, y) (list, list): Random walk of length n
    """
    x, y = [0], [0]
    deltas = [(1,0), (0,1), (-1,0), (0,-1)] #possible directions for the 2D lattice
    for _ in range(n):
        dx, dy = deltas[np.random.randint(0,4)]   #direction chosen at random
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
    return x, y

def plot_simpl_rnd_walk(n):
    """
    Plots a simple random walk of length n
    
    Args:
        n (int): the length of the walk
    Returns:
        Plot of a simple random walk of length n
    """
    x, y = simpl_rnd_walk(n)
    plt.figure(figsize = (8, 8))
    plt.plot(x, y, 'bo-', linewidth = 1)
    plt.plot(0, 0, 'go', ms = 12, label = 'Start')
    plt.plot(x[-1], y[-1], 'ro', ms = 12, label = 'End')
    plt.axis('equal')
    plt.legend(fontsize=15)
    plt.title('Figure 1: Simple random walk of length ' + str(n), fontsize=14, fontweight='bold', y = 1.05)
    plt.show()

plot_simpl_rnd_walk(40)

#SELF AVOIDING RANDOM WALK
#FIRST SIMULATION ALGORITHM
#PRINCIPLE
def myopic_saw(n):
    """
    Tries to generate a SAW of length n using the myopic algorithm
    
    Args:
        n (int): the length of the walk
    Returns:
        (x, y, stuck, steps) (list, list, bool, int): 
            (x,y) is a SAW of length <= n
            stuck is 1 if the walk could not terminate
            steps is the number of sites of the final walk
    """    
    x, y = [0], [0]
    positions = set([(0,0)])  #positions is a set that stores all sites visited by the walk
    stuck = 0
    for i in range(n):
        deltas = [(1,0), (0,1), (-1,0), (0,-1)]
        deltas_feasible = []  #deltas_feasible stores the available directions 
        for dx, dy in deltas:
            if (x[-1] + dx, y[-1] + dy) not in positions:  #checks if direction leads to a site not visited before
                deltas_feasible.append((dx,dy))
        if deltas_feasible:  #checks if there is a direction available
            dx, dy = deltas_feasible[np.random.randint(0,len(deltas_feasible))]  #choose a direction at random among available ones
            positions.add((x[-1] + dx, y[-1] + dy))
            x.append(x[-1] + dx)
            y.append(y[-1] + dy)
        else:  #in that case the walk is stuck
            stuck = 1
            steps = i+1
            break  #terminate the walk prematurely
        steps = n+1
    return x, y, stuck, steps

def plot_saw(n):
    """
    Plots the output of the myopic algorithm
    
    Args:
        n (int): the length of the walk
    Returns:
        Plot of the output of the myopic algorithm
    """    
    x, y, stuck, steps = myopic_saw(n)
    plt.figure(figsize = (8, 8))
    plt.plot(x, y, 'bo-', linewidth = 1)
    plt.plot(0, 0, 'go', ms = 12, label = 'Start')
    plt.plot(x[-1], y[-1], 'ro', ms = 12, label = 'End')
    plt.axis('equal')
    plt.legend()
    if stuck:
        plt.title('Figure 2: Walk stuck at step ' + str(steps), fontsize=14, fontweight='bold', y = 1.05)
    else:
        plt.title('Figure 2: SAW of length ' + str(n), fontsize=14, fontweight='bold', y = 1.05)
    plt.show()
    
plot_saw(80)

#ANALYSIS
def myopic_benchmark(n, samples):
    """
    Samples walks from the myopic algorithm and stores their length
    
    Args:
        n (int): the length of the walk
        samples (int): number of walks to sample from the myopic algorithm
    Returns:
        (fails, lengths) (list, list): 
            fails is a list of 0's and 1's. 1 if the walk was stuck
            lengths stores the length of each generated walk
    """    
    fails = []
    lengths = []
    for _ in range(samples):
        _, _, stuck, steps = myopic_saw(n)
        fails.append(stuck)
        lengths.append(steps)
    return (fails, lengths)

n_list = [10,30,50,70,100,200,500]
results_success, results_lengths = [], []
for n in n_list:
    fails, lengths = myopic_benchmark(n, 10000)
    results_success.append(100*(1-np.mean(fails)))
    results_lengths.append(lengths)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

bp = ax1.boxplot(results_lengths, meanline = True, showmeans = True)
ax1.set_xticklabels([str(n) for n in n_list])
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
ax1.set_ylim(0, 200)
ax1.set_title('Statistics on lengths of generated SAWs, for different n' 
              + '\n orange line = median, dotted green line = mean ')
#ax1.set_title('Statistiques sur les longueurs des MAAE générées, pour différents n' 
#             + '\n en orange la médiane, en pointillé vert la moyenne ')

ax2.set_xticks(n_list)
lp = ax2.plot(n_list, results_success, 'o-')
ax2.set_title('Proportion of generated SAWs that were not stuck, for different n')
#ax2.set_title('Proportion de MAAE non coincées, pour différents n')
ax2.yaxis.set_major_formatter(PercentFormatter())

fig.suptitle('Figure 3: Statistics for the myopic algorithm \n\n', fontsize=14, fontweight='bold', y = 1.05)

#SIMULATION VIA REJECTION SAMPLING
# PRINCIPLE
def is_saw(x, y, n):
    """
    Checks if walk of length n is self-avoiding
    
    Args:
        (x,y) (list, list): walk of length n
        n (int): length of the walk
    Returns:
        True if the walk is self-avoiding
    """    
    return n+1 == len(set(zip(x,y)))  #creating a set removes duplicates, so it suffices to check the size of the set

def rejection_sampling(n, samples):
    """
    Compute acceptance rate of rejection sampling
    
    Args:
        n (int): length of generated walks
        samples (int): number of simple walks generated 
    Returns:
        (float): acceptance rate
    """    
    accepted = 0
    for _ in range(samples):
        x, y = simpl_rnd_walk(n) #generates a simple random walk
        if is_saw(x, y, n): #checks if it is a SAW
            accepted += 1
    return accepted/samples

#ANALYSIS
def plot_acceptance(n, samples):
    """
    Plots the acceptance rates for lengths 1 to n
    
    Args:
        n (int): greatest length to consider
        samples (int): number of simple walks to generate for each length
    Returns:
        Plot of the acceptance rate as a function of n
    """    
    acceptance = []
    for i in range(1, n):
        acceptance.append(rejection_sampling(i, samples)*100)
    plt.figure(figsize = (8, 8))
    plt.plot(range(1, n), acceptance, 'o-')
    ax = plt.axes()
    ax.set_xticks(list(range(1,n)))
    ax.yaxis.set_major_formatter(PercentFormatter())
    plt.title('Figure 4: Acceptance rate of the rejection sampling method, for different n', 
              fontsize=14, fontweight='bold', y = 1.05)
    plt.show()

plot_acceptance(20,10000)

#ESTIMATION OF END TO END DISTANCE VIA IMPORTANCE SAMPLING
#PRINCIPLE
def saw_weights(n):
    """
    Tries to generate a SAW of length n using the myopic algorithm and keeps track of its weight
    
    Args:
        n (int): the length of the walk
    Returns:
        (x, y, stuck, weights) (list, list, bool, list): 
            (x,y) is a SAW of length <= n
            stuck is 1 if the walk could not terminate
            weights is a list of inverse probabilities 
    """
    #code is almost identical to myopic_saw(n)
    x, y = [0], [0] 
    positions = set([(0,0)])
    weights = []
    stuck = 0
    for i in range(n):
        deltas = [(1,0), (0,1), (-1,0), (0,-1)]
        deltas_feasible = []
        for dx, dy in deltas:
            if (x[-1] + dx, y[-1] + dy) not in positions:
                deltas_feasible.append((dx,dy))
        if deltas_feasible:
            weights.append(len(deltas_feasible))  #add inverse probability to weights
            dx, dy = deltas_feasible[np.random.randint(0,len(deltas_feasible))]
            positions.add((x[-1] + dx, y[-1] + dy))
            x.append(x[-1] + dx)
            y.append(y[-1] + dy)
        else:
            stuck = 1
            break
    return x, y, stuck, weights


def importance_sampling(n, samples):
    """
    Estimates squared end-to-end distance of SAWs of length n by importance sampling
    
    Args:
        n (int): length of walks
        samples (int): number of SAWs to sample from myopic algorithm
    Returns:
        (estim, std, ci_inf, ci_sup) (float, float, float ,float):
            estim is the mean squared end-to-end distance in the weighted sample
            std is the approximate standard error (see expression above)
            ci_inf and ci_sup are the bounds of the 95% approximate confidence interval
    """
    saw = 0
    weight_list = []  #stores weight of each SAW
    norm_list = []  #stores squared end-to-end distance of each SAW
    while saw < samples:  
        x, y, stuck, weights = saw_weights(n)
        if not stuck:
            weight_list.append(np.prod(weights))  #computes final weight of the walk as product of its weights
            norm_list.append(x[-1]**2 + y[-1]**2)  #computes squared end-to-end distance
            saw += 1
    weight_list = np.array(weight_list)
    norm_list = np.array(norm_list)
    estim = np.sum(weight_list * norm_list) / np.sum(weight_list)  #computes estimate (vectorized)
    std = np.sqrt(np.sum( (weight_list / np.sum(weight_list))**2 * (norm_list - estim)**2 )) #computes std (vectorized)
    ci_inf = estim - 1.96 * np.sqrt(1/samples * std)
    ci_sup = estim + 1.96 * np.sqrt(1/samples * std)
    return (estim, std, ci_inf, ci_sup)

#ANALYSIS
sq_dis = [2.66667,4.55556,7.04,9.56338,12.5744,15.5562,19.0128,22.4114,26.2425,30.0177,34.187,38.3043,42.7864,
         47.2177,51.9925,56.7164,61.7665,66.7658,72.0765] # harcoded list of true values of the mean squared distance

relative_error= []
for i in range(19):
    esti, _, _, _ = importance_sampling(i+2, 10000)
    relative_error.append(abs(esti-sq_dis[i])/sq_dis[i]*100)
plt.figure(figsize = (8, 8))
plt.plot(range(2, 21), relative_error, 'o-')
ax = plt.axes()
ax.set_xticks(list(range(2,21)))
ax.yaxis.set_major_formatter(PercentFormatter())
plt.title('Figure 5: Relative errors of the estimate of the squared end-to-end distance, for different n', 
              fontsize=14, fontweight='bold', y = 1.05)
plt.show()

#SIMULATION VIA PIVOT ALGORITHM
#DIMERIZATION

def small_saw(n): 
    """
    Generates a SAW of length n by rejection sampling, with early stopping (checks at each step if walk is 
    non-intersecting. If it intersects itself, stops prematurely)
    Will be used for n<=10, hence the name small_saw
    
    Args:
        n (int): the length of the walk
    Returns:
        (x, y) (list, list): SAW of length n
    """
    #early stopping, hence more lines than previous rejection sampling
    deltas = [(1,0), (0,1), (-1,0), (0,-1)]
    not_saw = 1
    while not_saw:
        x, y = [0], [0]
        positions = set([(0,0)])
        abort = 0
        i = 0
        while i < n and not(abort):
            dx, dy = deltas[np.random.randint(0,4)]
            if (x[-1] + dx, y[-1] + dy) in positions:
                abort = 1
                break
            else:
                x.append(x[-1] + dx)
                y.append(y[-1] + dy)
                positions.add((x[-1] + dx, y[-1] + dy))
                i = i+1
        if not(abort):
            not_saw = 0
    return x, y


def dimer(n):
    """
    Generates a SAW of length n by dimerization
    
    Args:
        n (int): the length of the walk
    Returns:
        (x, y) (list, list): SAW of length n
    """
    if n <= 3:
        x, y, _, _ = myopic_saw(n) #base case uses the myopic algorithm
        return x, y
    else:
        not_saw = 1
        while not_saw:
            (x_1, y_1) = dimer(n//2)  #recursive call
            (x_2, y_2) = dimer(n-n//2)  #recursive call
            x_2 = [x + x_1[-1] for x in x_2]  #translates the second walk to the end of the first one
            y_2 = [y + y_1[-1] for y in y_2]  #translates the second walk to the end of the first one
            x_concat, y_concat = x_1 + x_2[1:], y_1 + y_2[1:]  #performs concatenation
            if is_saw(x_concat, y_concat, n):   #if walk obtained is SAW, stop
                not_saw = 0
        return x_concat, y_concat
    
def plot_dimer(n):
    """
    Plots the output of the dimerization method
    
    Args:
        n (int): the length of the walk
    Returns:
        Plot of the output of the dimerization algorithm
    """    
    x, y = dimer(n)
    plt.figure(figsize = (8, 8))
    plt.plot(x, y, 'b.-', linewidth = 1)
    plt.plot(0, 0, 'go', ms = 12, label = 'Start')
    plt.plot(x[-1], y[-1], 'ro', ms = 12, label = 'End')
    plt.axis('equal')
    plt.legend()
    plt.title('Figure 5: SAW of length ' + str(n)+ ' generated by dimerization', 
              fontsize=14, fontweight='bold', y = 1.05)
    plt.show()
    

start_time = time.time()
plot_dimer(300)
print("--- %s seconds ---" % round((time.time() - start_time), 4))

#THE PIVOT ALGORITHM

#definition of elements of the dihedral group
m1, m2, m3, m4 = np.array([[1,0],[0,1]]), np.array([[0,-1],[1,0]]), np.array([[-1,0],[0,-1]]), np.array([[0,1],[-1,0]])
m5, m6, m7, m8 = np.array([[1,0],[0,-1]]), np.array([[0,1],[1,0]]), np.array([[-1,0],[0,1]]), np.array([[0,-1],[-1,0]])
m_list = [m1, m2, m3, m4, m5, m6, m7, m8]

def pivot(n, samples, method):
    """
    Generates SAWs of length n using the pivot algorithm
    For this algorithm, SAWs are represented by a 2*(n+1) matrix. 
    This makes it very simple to pivot walks (via matrix multiplication)
    
    Args:
        n (int): length of the SAWs
        samples (int): length of the Markov chain to be generated
        method (string): either dimerization or a straight rod
    Returns:
        results (list): list of the resulting SAWs
    """
    results = []
    #initialization
    if method == 'dimer':
        x, y = dimer(n)
    else:
        x, y = range(0,n+1), [0 for _ in range(n+1)]  #straight horizontal rod
    walk = np.c_[x,y].T
    results.append(walk)
    #main step
    for _ in range(samples):
        i = np.random.randint(0,n)  #choose the location of the pivot uniformly at random
        m = m_list[np.random.randint(0,8)]  #choose transformation in the dihedral group uniformly at random
        walk_1, walk_2 = walk[:, 0:i+1], walk[:, i+1:]  #split walk in 2 parts
        pivot = np.vstack(walk[:, i])  #site that will be used as pivot
        walk_2 = np.dot(m, (walk_2 - pivot)) + pivot   #transform second part of the walk around the pivot
        walk_piv = np.c_[walk_1,walk_2]  #concatenate
        if np.unique(walk_piv, axis=1).shape[1] == (n+1):   #check if resulting walk is SAW
            results.append(walk_piv)  
            walk = walk_piv
        else:
            results.append(walk)
    return results

def plot_pivot(n, samples, method, i, j):
    """
    Plots SAWs generated by the pivot algorithm
    
    Args:
        n (int): length of SAWs
        samples (int): length of the Markov chain
        method (string): either dimerization or a straight rod
        i, j (int, int): plot all walks between the i-th and the j-th 
    Returns:
        Plot of SAWs generated by the pivot algorithm
    """    
    plt.figure(figsize = (10, 8))
    results = pivot(n, samples, method)
    for k in range(i,j):
        walk = results[k]
        x, y = walk[0], walk[1] 
        plt.plot(x, y, '.-', linewidth = 1)
        plt.plot(x[-1], y[-1], 'ro', ms = 12)
    plt.plot(0, 0, 'ro', ms = 12, label = 'End')
    plt.plot(0, 0, 'go', ms = 12, label = 'Start')
    plt.axis('equal')
    plt.legend()
    if method == 'dimer':
        plt.title(str(j-i) + ' SAWs of length ' + str(n) + ' generated by pivot,\n'+
                  'initialization by dimerization', 
              fontsize=14, fontweight='bold', y = 1.05)
    else:
        plt.title(str(j-i) + ' SAWs of length ' + str(n) + ' generated by pivot,\n'+
                  'initialization by a horizontal rod', 
              fontsize=14, fontweight='bold', y = 1.05)
    plt.show()

plot_pivot(100, 50, 'dimer', 0, 20)

plot_pivot(100, 10000, 'rod', 9950, 9970)

#Use of Pivot for Estimation

def pivot_estimate(n, samples, method):
    """
    Estimates squared end-to-end distance of SAWs of length n by pivot
    
    Args:
        n (int): length of walks
        samples (int): length of the Markov chain
    Returns:
        results (list): list that stores squared end-to-end distances
    """
    results = []
    if method == 'dimer':
        x, y = dimer(n)
        dist = x[-1]**2 + y[-1]**2
    else:
        x, y = range(0,n+1), [0 for _ in range(n+1)]
        dist = x[-1]**2 + y[-1]**2
    walk = np.c_[x,y].T
    results.append(dist)
    for _ in range(samples):
        i = np.random.randint(0,n)
        m = m_list[np.random.randint(0,8)]
        walk_1, walk_2 = walk[:, 0:i+1], walk[:, i+1:]
        pivot = np.vstack(walk[:, i])
        walk_2 = np.dot(m, (walk_2 - pivot)) + pivot
        walk_piv = np.c_[walk_1,walk_2]
        if np.unique(walk_piv, axis=1).shape[1] == (n+1):
            end = walk_piv[:,-1]
            results.append(np.inner(end, end))
            walk = walk_piv
        else:
            results.append(results[-1])
    return results

sq_dis = [2.66667,4.55556,7.04,9.56338,12.5744,15.5562,19.0128,22.4114,26.2425,30.0177,34.187,38.3043,42.7864,
         47.2177,51.9925,56.7164,61.7665,66.7658,72.0765] # harcoded list of true values of the mean squared distance

relative_error= []
for i in range(19):
    esti = np.mean(pivot_estimate(i+2, 10000, 'dimer'))
    relative_error.append(abs(esti-sq_dis[i])/sq_dis[i]*100)
plt.figure(figsize = (8, 8))
plt.plot(range(2, 21), relative_error, 'o-')
ax = plt.axes()
ax.set_xticks(list(range(2,21)))
ax.yaxis.set_major_formatter(PercentFormatter())
plt.title('Figure 6: Relative errors of the estimate of the squared end-to-end distance using pivot, for different n', 
              fontsize=14, fontweight='bold', y = 1.05)
plt.show()

