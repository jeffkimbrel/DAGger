import networkx as nx


nodes = {'M': ['A-T-G'],
        'Q': ['C-A-A', 'C-A-G'],
        'T': ['A-C-A', 'A-C-C', 'A-C-G', 'A-C-T'],
        'L': ['C-T1-A', 'C-T1-G', 'C-T1-T', 'C-T1-C', 'T-T2-A', 'T-T2-G'],
        'K': ['A-A-G', 'A-A-A'],
        'V': ['G-T-A', 'G-T-C', 'G-T-G', 'G-T-T'],
        'D': ['G-A-C', 'G-A-T']
}

ys = {'A': 3,
    'C': -1,
    'T': -3,
    'G': 1,
    'T1': -3.5,
    'T2': 2

}

peptide = 'MQTLKVD'

G = nx.DiGraph()

last = []
pos = {}

for aa_pos, aa in enumerate(peptide):
    nt_pos = aa_pos * 3 + 1
    

    new_last = []

    for codon in nodes[aa]:
        nt = codon.split("-")

        pos[f'{nt_pos+0}{nt[0]}'] = (nt_pos+0, ys[nt[0]])
        pos[f'{nt_pos+1}{nt[1]}'] = (nt_pos+1, ys[nt[1]])
        pos[f'{nt_pos+2}{nt[2]}'] = (nt_pos+2, ys[nt[2]])

        if len(last) > 0:
            for last_node in last:
                G.add_edge(f'{last_node}', f'{nt_pos+0}{nt[0]}')

        G.add_edge(f'{nt_pos+0}{nt[0]}', f'{nt_pos+1}{nt[1]}')
        G.add_edge(f'{nt_pos+1}{nt[1]}', f'{nt_pos+2}{nt[2]}')

        new_last.append(f'{nt_pos+2}{nt[2]}')
    
    last = new_last




import matplotlib.pyplot as plt

# print(pos)





#pos = nx.spring_layout(G)
nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'))
nx.draw_networkx_labels(G, pos)
nx.draw_networkx_edges(G, pos, arrows=True)
# nx.draw_networkx_edges(G, pos, edgelist=black_edges, arrows=False)
plt.show()