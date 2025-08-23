import argparse
import glob
import itertools
import pickle

from dash import Dash, html, dash_table, dcc, Output, Input, State
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import pandas as pd


# Only column names listed here will be displayed.
# Consult an example of the clusterProfiler output file for the
# full set of available names
GO_DISPLAY_COLUMNS = [
    'ID',
    'Description',
    'GeneRatio',
    'BgRatio',
    'p.adjust'
]

INIT_PADJ_THRESHOLD = 0.05
global_padj_threshold = INIT_PADJ_THRESHOLD

# default styling for the nodes
BASE_STYLESHEET = [
    {
        'selector': '[sig_down = 1]',
        'style': {'background-color': "#788bb8", 'label': 'data(label)'}
    },
    {
        'selector': '[sig_up = 1]',
        'style': {'background-color': "#cc761f", 'label': 'data(label)'}
    },
    {
        'selector': '[not_sig = 1]',
        'style': {'background-color': "#a1a1a1", 'label': 'data(label)'}
    },
    {
        'selector': '[not_found = 1]',
        'style': {'background-color': "#ffffff", 'label': 'data(label)'}
    },
    {
        'selector': 'edge',
        'style': {'line-color': '#A3A3A3', 'width': 2}
    }
]


def parse_cl_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', 
                        '--results_dir', 
                        required=True, 
                        help='Path to the results directory (the timestamped directory)')
    parser.add_argument('-t', 
                        '--network_type', 
                        choices = ['dge_only', 'full'],
                        required=True, 
                        help='Path to file of differential expression results')
    parser.add_argument('--debug', action='store_true')
    return parser.parse_args()


def nx_to_cytoscape_elements(G):
    '''
    Helper function that takes a nx.Graph object and
    converts it to cytoscape-compatible format
    '''
    elements = []

    # Add nodes
    for node, data in G.nodes(data=True):
        elements.append({'data': {'id': str(node), 'label': str(node), **data}})
    # Add edges
    for source, target, data in G.edges(data=True):
        elements.append({
            'data': {
                'source': str(source),
                'target': str(target),
                **data
            }
        })
    return elements


def parse_gmt(gmt_file):
    '''
    Parses a GMT-format file (e.g. MsigDB format) 
    and returns a dict giving the name of the gene set
    pointing at a list of gene symbols
    '''
    d = {}
    for line in open(gmt_file):
        line = line.strip()
        contents = line.split('\t')
        gene_set_name = contents[0]
        gene_set = contents[2:]
        d[gene_set_name] = gene_set
    return d


def update_elements(elements):
    '''
    Sets attributes on the nodes. Those attributes ultimately
    determine how the nodes are displayed.
    '''
    for item in elements:
        # item is a dict with 'data' and 'position' keys
        if global_padj_threshold is None:
            padj_threshold = INIT_PADJ_THRESHOLD
        else:
            padj_threshold = global_padj_threshold

        info = item['data']
        if 'label' in info:
            try:
                info['sig_up'] = int((info['padj'] <= padj_threshold) & (info['log2FoldChange'] > 0))
                info['sig_down'] = int((info['padj'] <= padj_threshold) & (info['log2FoldChange'] < 0))
                info['not_sig'] = int(info['padj'] > padj_threshold)
            except KeyError:
                info['sig_up'] = 0
                info['sig_down'] = 0
                info['not_found'] = 1


def prep_app(dge_file, gmt_file, graph_pickle, h5_go_file):
    '''
    Sets up the Dash app. This includes reading and reformatting
    input files as well as setting up the application structure and
    callbacks. Returns a Dash.app instance
    '''

    # First we need to read the input files and massage them a bit to be
    # useful for the app

    dge_df = pd.read_table(dge_file, index_col=0)
    gmt_dict = parse_gmt(gmt_file)

    with open(graph_pickle, 'rb') as fin:
        G = pickle.load(fin) 

    # a list of all the genes in all the communities 
    gene_list = list(itertools.chain.from_iterable(gmt_dict.values()))

    # Define a reverse lookup table so we can find the community that is
    # associated with a gene.
    # Since we are working with communities, each gene is in only a single
    # community and hence the gene-to-set is 1:1
    lookup_dict = {v:gs for gs, values in gmt_dict.items() for v in values}

    gene_set_names = list(gmt_dict.keys())
    gene_set_options = [
        {'label':f'{gene_set_name} ({len(gmt_dict[gene_set_name])} genes)', 'value': gene_set_name} 
        for gene_set_name in gene_set_names]

    # 
    for node in G.nodes:
        try:
            gene_info = dge_df.loc[node].to_dict()
            gene_info['sig_up'] = int((gene_info['padj'] <= INIT_PADJ_THRESHOLD) & (gene_info['log2FoldChange'] > 0))
            gene_info['sig_down'] = int((gene_info['padj'] <= INIT_PADJ_THRESHOLD) & (gene_info['log2FoldChange'] < 0))
            gene_info['not_sig'] = int(gene_info['padj'] > INIT_PADJ_THRESHOLD)
        except:
            gene_info = {}
        G.nodes[node].update(gene_info)
   
    # Now that the input data is loaded and ready, we initialize the app
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    # Define the layout/elements of the app here:
    app.layout = [
        html.H1('Network Browser'),

        # A section containing dropdowns to search by the community or by a gene symbol
        html.Div([
            html.Label('Select a community:', htmlFor='set_select'),
            dcc.Dropdown(options=gene_set_options, value=gene_set_names[0], id='set_select'),
            html.Div('-or-'),
            html.Label('Search for a gene and display the associated network', htmlFor='gene_dropdown'),
            dcc.Dropdown(
                id='gene_dropdown',
                options=[{'label': i, 'value': i} for i in gene_list],
                placeholder='Select or type a gene',
                searchable=True,
                clearable=True
            ),
        ], style={'width': '25%', 'padding': '20px'}),

        # The cytoscape viz, which includes some controls/tooltips
        html.Div([
            # the actual network
            dcc.Loading([
            cyto.Cytoscape(
                id='cytoscape',
                responsive=True,
                elements=[],
                boxSelectionEnabled=True, 
                style={'width': '900px', 
                    'height': '600px', 
                    'margin-right': '10px',
                    'padding': '10px',
                    'border': '3px solid gray',
                        'borderRadius': '10px'},
                layout={'name': 'cose'},
                stylesheet=BASE_STYLESHEET
            )], overlay_style={"visibility":"visible", "opacity": .5, "backgroundColor": "white"}),
            # a div containing controls
            html.Div([
                html.Label('Significance threshold (adjusted p-value):', htmlFor='padj-filter'),
                html.Div(id='padj-filter', children=[
                    dcc.Input(id='padj-filter-input', type='number', value=INIT_PADJ_THRESHOLD, min=0, max=1.0),
                    dbc.Button('Apply filter', id='submit-padj')
                ], style={'display': 'flex', 'flexDirection': 'row'}),
                html.Label('Hover info:', htmlFor='hover-tooltip'),
                html.Div(id='hover-tooltip', style={
                                                'padding': '12px', 
                                                'minHeight': '40px',
                                                'background': '#f8f8f8', 
                                                'border': '1px solid #ccc',
                                                'borderRadius': '5px', 
                                                'marginTop': '16px', 
                                                'margin': 0,
                                                'display': 'inline-block'}),
                html.Label('Selection info:', htmlFor='selected-tooltip'),
                html.Div(id='selected-tooltip', style={
                                    'padding': '12px', 
                                    'minHeight': '40px',
                                    'background': '#f8f8f8', 
                                    'border': '1px solid #ccc',
                                    'borderRadius': '5px', 
                                    'marginTop': '16px', 
                                    'margin': 0,
                                    'display': 'inline-block'}),
            ], style={'display': 'flex', 'flexDirection': 'column'}),
            ], style={
                'display': 'flex',
                'flexDirection': 'row',
                'alignItems': 'flex-start',
                'justifyContent': 'flex-start',
                'padding': '0',
                'margin': '0',
                'gap': '0'
        }),
        dbc.RadioItems(id='ontology_select', options=['BP', 'MF', 'CC'], value='BP', inline=True),
        dash_table.DataTable(id='go-table', data=[], page_size=10, style_table={'width': '50%', 'margin-top': '20px'})
    ]

    # When the user hovers over a node, give some additional info
    @app.callback(
        Output('hover-tooltip', 'children'),
        Input('cytoscape', 'mouseoverNodeData'),
    )
    def display_node_hover_tooltip(data):
        if data is None:
            return "Hover over a node to see more information."

        return html.Div([
            html.B(data.get('label','')),
            html.Br(),
            f'Log fold-change: {data["log2FoldChange"]:.2f}',
            html.Br(),
            f'Adjusted p-value: {data["padj"]:.5f}'
        ])

    # For updating the GO enrichment table
    @app.callback(
        Output(component_id='go-table', component_property='data', allow_duplicate=True),
        Input(component_id='set_select', component_property='value'),
        prevent_initial_call=True
    )
    def update_table(selected):
        if selected is not None:
            with pd.HDFStore(h5_go_file) as h5:
                key = f'/BP/{selected}'
                if key in h5:
                    _df = h5[key]
                    return _df[GO_DISPLAY_COLUMNS].to_dict('records')
                else:
                    return pd.DataFrame().to_dict('records')


    # For updating the GO enrichment table by ontology
    @app.callback(
        Output(component_id='go-table', component_property='data'),
        Input(component_id='ontology_select', component_property='value'),
        State('set_select', 'value'),
    )
    def update_ontology(selected_ontology, seletected_set):
        if selected_ontology is not None:
            with pd.HDFStore(h5_go_file) as h5:
                key = f'/{selected_ontology}/{seletected_set}'
                if key in h5:
                    _df = h5[key]
                    return _df[GO_DISPLAY_COLUMNS].to_dict('records')
                else:
                    return pd.DataFrame().to_dict('records')


    # For updating the network graph
    @app.callback(
        Output(component_id='cytoscape', component_property='elements', allow_duplicate=True),
        Input(component_id='set_select', component_property='value'),
        prevent_initial_call=True
    )
    def update_fig(selected):
        if selected is not None:
            gene_set = gmt_dict[selected]
            G2 = G.subgraph(gene_set)
            elements = nx_to_cytoscape_elements(G2)
            update_elements(elements)
            return elements
        else:
            return []

    # Once a user selects a gene, this resets the community selection
    # dropdown. That, in turn, will change the network graph
    @app.callback(
        Output('set_select', component_property='value'),
        Input('gene_dropdown', 'value')
    )
    def search_gene(gene):
        if gene is not None:
            selected_gene_set = lookup_dict[gene]
            return selected_gene_set

    @app.callback(
        Output('selected-tooltip', 'children'),
        Input('cytoscape', 'tapNodeData'),
    )
    def display_node_click_tooltip(data):
        if data is None:
            return "Click a node to see more information."

        return html.Div([
            html.B(data.get('label','')),
            html.Br(),
            f'Log fold-change: {data["log2FoldChange"]:.2f}',
            html.Br(),
            f'Adjusted p-value: {data["padj"]:.5f}'
        ])

    # When a user clicks on a node, highlight the node and edges
    @app.callback(
        Output('cytoscape', 'stylesheet'),
        Input('cytoscape', 'tapNodeData')
    )
    def highlight_edges_on_node_click(node):
        if node is None:
            return BASE_STYLESHEET
        
        node_id = node['id']

        # Highlight the edges
        highlight_styles = [
            {
                'selector': f'edge[source = "{node_id}"], edge[target = "{node_id}"]',
                'style': {
                    'line-color': 'black',
                    'width': 3,
                    'z-index': 9999
                }
            }
        ]
        # Highlight the node itself
        highlight_styles += [
            {
                'selector': f'node[id = "{node_id}"]',
                'style': {'border-width': 3, 'border-color': 'black'}
            }
        ]
        return BASE_STYLESHEET + highlight_styles

    # callback for when a user sets a adjusted p-value threshold.
    # Will update the node colors to indicate significance at that threshold
    @app.callback(
        Output(component_id='cytoscape', component_property='elements'),
        Input('submit-padj', 'n_clicks'),
        State('padj-filter-input', 'value'),
        State('cytoscape', 'elements')
    )
    def update_padj(n_clicks, value, cytoscape_elements):
        global global_padj_threshold
        global_padj_threshold = value
        update_elements(cytoscape_elements)
        return cytoscape_elements
    
    return app


def find_and_assert_singleton(path_glob):
    files = glob.glob(path_glob)
    if len(files) == 1:
        return files[0]
    else:
        raise Exception(f'Found more than one file matching the pattern: {path_glob}.'
                        ' We only expect a single file to match.')


if __name__ == '__main__':
    args = parse_cl_args()
    results_dir = args.results_dir
    dge_file = find_and_assert_singleton(f'{args.results_dir}/*/dge_results/*.symbol_remapped.tsv')
    gmt_file = find_and_assert_singleton(f'{args.results_dir}/*/ppi_networks/*.{args.network_type}.gmt')
    graph_pickle = find_and_assert_singleton(f'{args.results_dir}/*/ppi_networks/{args.network_type}.pkl')
    h5_go_file = find_and_assert_singleton(f'{args.results_dir}/*/ppi_networks/*.{args.network_type}.h5')
    app = prep_app(dge_file, gmt_file, graph_pickle, h5_go_file)
    app.run(host="0.0.0.0", port=8050, debug=args.debug)

