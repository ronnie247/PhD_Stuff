import numpy as np
import plotly.graph_objects as go

def grid(m, n):
    y, x = np.indices((m, n))
    return x, y

def set_edges(x, y, p_given=0.5):
    if x.shape != y.shape:
        raise ValueError('x and y should have the same shape')
    m, n = x.shape
    I, J = np.where(np.random.rand(m-1, n-1)>=p_given)
    K, L = np.where(np.random.rand(m-1, n-1)>=p_given)

    x_edges = []
    y_edges = []
    for i, j  in zip(I, J):
        x_edges.extend([x[i,j], x[i, j+1], None])
        y_edges.extend([y[i,j], y[i,j], None])

    for k, l in zip(K,L):
        x_edges.extend([x[k,l], x[k, l], None])
        y_edges.extend([y[k,l], y[k+1, l], None])    
             
    return x_edges, y_edges       
    

x, y = grid(12, 12)
nodes = go.Scatter(name='node', 
                   x=x.flatten(), 
                   y=y.flatten(), 
                   mode='markers', 
                   marker_size=7, 
                   marker_color='black')

x_edges, y_edges = set_edges(x, y, p_given=0.8)
edge_trace = go.Scatter(name='edge', 
                        x=x_edges, 
                        y=y_edges, 
                        mode='lines', 
                        line_width=2, 
                        line_color='black')

fig = go.Figure(data=[edge_trace, nodes])
fig.update_layout(width=500, height=500, template='none', xaxis_visible=False, yaxis_visible=False)
fig.show()