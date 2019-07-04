import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_table

import matplotlib.pyplot as plt

import networkx as nx
import plotly.graph_objs as go

import pandas as pd
import numpy as np

import pickle

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

data_dir = './data/'
################################
# dataframe
data_fname = data_dir + 'data.csv'
df = pd.read_csv(data_fname)

def obtain_mean_std(ID):
  m = np.nan; std = np.nan
  m = df[df['strainID']==ID]['Assay Raw'].mean()
  std = df[df['strainID']==ID]['Assay Raw'].std()
  return m, std

##################################
# ancestry graph

ancestry_fname = data_dir + 'graph.txt'
G = nx.read_edgelist(ancestry_fname, create_using=nx.DiGraph())

def obtain_ancestry_path(G, strainID):
  """Given a strain ID, obtain its ancestors in reverse chronical order along with 
  the plasmids to a list.
  If the strainID is not in the graph G, then an empty list [] will be returned."""
  ancestry_path = []
  
  try:  
    ancestor_dic = nx.dfs_successors(G, strainID)

    key = strainID
    while key:
      v = ancestor_dic.get(key)
      if v:  # This key does have an ancestor
        plamsids_dic = G.edges[key, v[0]]
        plasmids = []
        if bool(plamsids_dic): # check the dic is not empty
          plasmids = plamsids_dic['plasmids']
        ancestry_path.append((plasmids, v[0]))
        key = v[0]
      else:
        key = v
  except KeyError:
      pass  
  
  return ancestry_path


#######################################
# model

f = open(data_dir + 'model.pkl', 'rb')
model = pickle.load(f)
feature_columns = pickle.load(f)
f.close()

pos_plasmids_filter = model.coef_[11:] > 0
pos_plasmids = np.array(feature_columns[11:])[pos_plasmids_filter]
pos_plasmid_coefs = model.coef_[11:][pos_plasmids_filter]
ppc, pp = zip(*sorted(zip(pos_plasmid_coefs, pos_plasmids), reverse=True))

def get_yield(strain, plasmid_lst, feature_columns):
  features = [0] * len(feature_columns)
  features[feature_columns.index(strain)] = 1
  for x in plasmid_lst:
    ind = feature_columns.index(x)
    features[ind] = 1
  features = np.array(features).reshape(1, -1)
  y = model.predict(features)[0]
  return y

#######################################
# engineering plan
def get_origin_plasmids(G, ID):
  ancestry_path = obtain_ancestry_path(G, ID)
  plasmids = []
  for x in ancestry_path:
    plasmids += x[0]
  origin = None
  if len(ancestry_path) > 0:
    origin = ancestry_path[-1][1]
  return origin, plasmids

def obtain_engineering_plan(G, ID, pp, ppc):
  current_predicted_yield = np.nan
  columns = ['genetic component', 'increased yield (use this component)', 
            'increased yield (use all previous + this component)', 
            'total yield']
  df = pd.DataFrame(columns=columns)

  origin, plasmids = get_origin_plasmids(G, ID)

  if origin in feature_columns[:11]:
    plasmid_lst = []
    for p in plasmids:
      if p in feature_columns[11:]:
        plasmid_lst.append(p)  
    current_predicted_yield = get_yield(origin, plasmid_lst, feature_columns)

    this_pp = []; this_ppc = []
    for i, p in enumerate(pp):
      if p not in plasmids:
        this_pp.append(p)
        this_ppc.append(ppc[i])
          
    df['genetic component'] = this_pp
    df['increased yield (use this component)'] = this_ppc
    
    c1 = []; c2 = []
    prev = 0; prev_total = current_predicted_yield
    for x in this_ppc:
      this_yield = x + prev
      this_total = x + prev_total
      c1.append(this_yield)
      c2.append(this_total)
      prev = this_yield
      prev_total = this_total
        
    df['increased yield (use all previous + this component)'] = c1
    df['total yield'] = c2

    df = df.round(3)
  
  return current_predicted_yield, df


#######################################
# candidate plan
f = open(data_dir + 'dic.pkl', 'rb')
dic = pickle.load(f)
f.close()

def obtain_combos(n_plasmids, dic):
  yields = sorted(dic.get(n_plasmids, {}).keys(), reverse=True)
  df = pd.DataFrame(columns=[
    'background microbe', 'genetic components', 'expected yield'])
  c1 = []; c2 = []; c3 = []
  for y in yields:
    combos = dic[n_plasmids][y]
    for combo in combos:
      c1.append(combo[0])
      c2.append(', '.join(combo[1]))
      c3.append(y)
  df['background microbe'] = c1
  df['genetic components'] = c2
  df['expected yield'] = c3
  return df.round(3)

######################################
# app
app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([

  html.H1('Your Microbe Dashboard'),

  html.Div([
    html.Div([
      html.H2('Microbe Comparison'), 
      dcc.Markdown('Input microbe ID here:'), 
      dcc.Input(id='input1', value='M756', type='text'),
      dcc.Input(id='input2', value='M806', type='text'), 
      html.Div(id='output1'), 
      html.Div(id='output2'), 
      dcc.Graph(id='comparison-bar-chart')
    ], className='four columns'), 

    html.Div([
      html.H2('Microbe History'), 
      dcc.Markdown('Input microbe ID here:'), 
      dcc.Input(id='input3', value='M12', type='text'),
      dcc.Graph(
        id='ancestry-path', 
      ), 
    ], className='four columns'), 

    html.Div([
      html.H2('Which microbe to build next?'), 
      dcc.Markdown('Input the background microbe here:'), 
      dcc.Input(id='bg-strain', value='M995', type='text'), 
      html.Div('', style={'padding': 10}), 
      dcc.Markdown('Input the genetic components you would like to introduce to the background microbe, separated by "," '), 
      dcc.Input(id='plasmids', value='G72,G96', type='text'), 
      html.H6(id='expected-yield', style={'padding': 10}),
    ], className='four columns')
  ], className='row'), 

  html.Div([
    html.Div([
      html.H2('Microbe Engineering Plan'), 
      html.H4('-- Limited amount of genetic components'),
      dcc.Markdown('Input the number of genetic components that you would like to introduce into an existing microbe:'), 
      dcc.Input(id='n-plasmids', value=1, type='number'), 
      html.Div('', style={'padding': 10}), 
      dcc.Markdown('Input the number of top results to be returned:'),   
      dcc.Input(id='n-results', value=10, type='number'), 
      html.Div('', style={'padding': 10}),    
      dash_table.DataTable(
        id='table-plasmid', 
        style_cell_conditional=[
          {'if': {'column_id': 'background microbe'},
          'width': '200px'},
          {'if': {'column_id': 'expected yield'},
          'width': '150px'},
        ], 
        style_table={'width': '1000px'}, 
      )
    ], className='row'),

    html.Div([
      html.H2('Microbe Engineering Plan'), 
      html.H4('-- Start from a current microbe'),
      dcc.Markdown('Input microbe ID here:'), 
      dcc.Input(id='plan-input', value='M995', type='text'), 
      html.Div('', style={'padding': 10}), 
      dcc.Markdown('Input the number of top results to be returned:'), 
      dcc.Input(id='n-results2', value=10, type='number'), 
      html.Div('', style={'padding': 10}), 
      html.Div(id='predicted_yield'), 
      html.Div('', style={'padding': 10}), 
      dash_table.DataTable(
        id='engineer-table', 
        style_table={'width': '1000px'}, 
      )
    ], className='row'),
  ], className='row')
   

], )

@app.callback(
  Output(component_id='output1', component_property='children'),
 
  inputs=[Input('input1', 'value')] 
)
def update_output1(input1):
  m, std = obtain_mean_std(input1)
  return u'{}: the overall yield is {:.3f} \u00B1 {:.3f}.'.format(input1, m, std)

@app.callback(
  Output(component_id='output2', component_property='children'),
  
  inputs=[Input('input2', 'value')] 
)
def update_output2(input2):
  m, std = obtain_mean_std(input2)
  return u'{}: the overall yield is {:.3f} \u00B1 {:.3f}.'.format(input2, m, std)

@app.callback(
  Output('comparison-bar-chart', 'figure'), 
  inputs=[Input('input1', 'value'), Input('input2', 'value')]
)
def create_bar_chart(input1, input2):
  m1, std1 = obtain_mean_std(input1)

  m2, std2 = obtain_mean_std(input2)

  data = [go.Bar(
    x = [input1, input2], 
    y = [m1, m2], 
    error_y=dict(
      type='data', 
      array=[std1, std2], 
      visible=True
    )
  )]

  fig = go.Figure(
    data=data, 
    layout=go.Layout(
      title='Yield comparison'
    )
  )

  return fig

@app.callback(
  Output('ancestry-path', 'figure'), 
  inputs=[Input('input3', 'value')]
)
def update_figure(strainID):
  ancestry_path = obtain_ancestry_path(G, strainID)
  n = len(ancestry_path) + 1

  trace = go.Scatter(
    x=[1] * n, 
    y=list(range(1, n+1)), 
    mode='markers+text', 
    marker=dict(size=[30]*n),
    text=[strainID] + [elt[1] for elt in ancestry_path], 
    textposition='middle right'
  )

  # edges
  x0 = [1] * (n-1)
  y0 = list(range(2, n+1))
  x1 = [1] * (n-1)
  y1 = list(range(1, n))

  edge_creation = [
        dict(ax=x0[i], ay=y0[i], axref='x', ayref='y',
          x=x1[i], y=y1[i], xref='x', yref='y') for i in range(0, len(x0))
      ]

  shift = 0
  edge_annotations = go.Scatter(
    x=[1 + shift] * len(x0), 
    y=[(y0[i] + y1[i])/2 for i in range(len(x0))], 
    mode='text', 
    text=[', '.join(elt[0]) for elt in ancestry_path], 
    textposition = 'middle right'
  )

  fig = go.Figure(
    data=[trace, edge_annotations],
    layout=go.Layout(
      title='The microbe history', 
      hovermode='closest', 
      # annotations = edge_annotations + text_annotations,
      annotations = edge_creation, 
      xaxis=dict(showgrid=False, zeroline=False, showticklabels=False), 
      yaxis=dict(showgrid=False, zeroline=False, showticklabels=False), 
      showlegend=False
    )
  ) 

  return fig

@app.callback(
  Output('expected-yield', 'children'), 
  inputs=[Input('bg-strain', 'value'), Input('plasmids', 'value')], 
)
def update_yield(strain, plasmids):
  origin = strain
  plasmid_lst = []
  for p in plasmids.split(','):
    if p in feature_columns[11:]:
      plasmid_lst.append(p)    

  if strain not in feature_columns[:11]:
    origin, original_plasmids = get_origin_plasmids(G, strain)
    for p in original_plasmids:
      if p in feature_columns[11:]:
        plasmid_lst.append(p)
  
  if origin:
    y = get_yield(origin, plasmid_lst, feature_columns)
  else:
    y = np.nan
  return 'The expected yield for this combination is {:.4f}.'.format(y)

@app.callback(
  [Output(component_id='table-plasmid', component_property='columns'), 
  Output(component_id='table-plasmid', component_property='data')], 
  [Input('n-plasmids', 'value'), Input('n-results', 'value')], 
)
def get_candidate_table(n_plasmids, n_results):
  d = obtain_combos(n_plasmids, dic)
  columns=[{"name": i, "id": i} for i in d.columns]
  data = d.head(n_results).to_dict('records')
  return columns, data

@app.callback(
  [Output(component_id='engineer-table', component_property='columns'), 
  Output(component_id='engineer-table', component_property='data'), 
  Output('predicted_yield', 'children')], 
  [Input('plan-input', 'value'), Input('n-results2', 'value')], 
)
def get_engineering_table(ID, n_results):
  current_y, d = obtain_engineering_plan(G, ID, pp, ppc)
  columns=[{"name": i, "id": i} for i in d.columns]
  data = d.head(n_results).to_dict('records')
  return columns, data, 'The current predicted yield of the input microbe is {:.3f}.'.format(current_y)

if __name__ == '__main__':
  app.run_server()