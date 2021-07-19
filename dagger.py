import matplotlib.pyplot as plt
import networkx as nx
import re


def pep_to_graph(peptide, shift=0):

    G = nx.DiGraph()

    last = []
    pos = {}

    for aa_pos, aa in enumerate(peptide):
        nt_pos = aa_pos * 3 + 1 + shift

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

    for i in list(G.nodes()):
        G.nodes[i]['NT'] = re.findall(r"\D", i)

    return (G, pos)


nodes = {'A': ['G-C-A', 'G-C-C', 'G-C-G', 'G-C-T'],
         'C': ['T-G-T', 'T-G-C'],
         'D': ['G-A-C', 'G-A-T'],
         'E': ['G-A-A', 'G-A-G'],
         'F': ['T-T-C', 'T-T-T'],
         'G': ['G-G-A', 'G-G-C', 'G-G-G', 'G-G-T'],
         'H': ['C-A-C', 'C-A-T'],
         'I': ['A-T-A', 'A-T-C', 'A-T-T'],
         'K': ['A-A-G', 'A-A-A'],
         # 'L': ['C-T1-A', 'C-T1-G', 'C-T1-T', 'C-T1-C', 'T-T2-A', 'T-T2-G'],
         'L': ['C-T-A', 'C-T-G', 'C-T-T', 'C-T-C', 'T-T-A', 'T-T-G'],
         'M': ['A-T-G'],
         'N': ['A-A-C', 'A-A-T'],
         'P': ['C-C-A', 'C-C-C', 'C-C-G', 'C-C-T'],
         'Q': ['C-A-A', 'C-A-G'],
         # 'R': ['C-G1-A', 'C-G1-C', 'C-G1-G', 'C-G1-T', 'A-G2-G', 'A-G2-A'],
         'R': ['C-G-A', 'C-G-C', 'C-G-G', 'C-G-T', 'A-G-G', 'A-G-A'],
         'S': ['T-C-A', 'T-C-C', 'T-C-G', 'T-C-T', 'A-G-C', 'A-G-T'],
         'T': ['A-C-A', 'A-C-C', 'A-C-G', 'A-C-T'],
         'V': ['G-T-A', 'G-T-C', 'G-T-G', 'G-T-T'],
         'W': ['T-G-G'],
         'Y': ['T-A-C', 'T-A-T'],
         }

ys = {'A': 3,
      'C': -1,
      'T': -3,
      'G': 1,
      'T1': -3.5,
      'T2': 2,
      'G1': 0,
      'G2': 2

      }

def getMCS(g1, g2):
    # https://stackoverflow.com/questions/43108481/maximum-common-subgraph-in-a-directed-graph

    matching_graph = nx.Graph()

    for n1, n2 in g2.edges():
        if g1.has_edge(n1, n2):
            matching_graph.add_edge(n1, n2)
            # print(n1, n2)

    components = nx.connected_components(matching_graph)

    largest_component = max(components, key=len)
    # print(len(largest_component))
    return nx.induced_subgraph(matching_graph, largest_component), largest_component
    # return nx.induced_subgraph(g1, largest_component)
    # return matching_graph


first_pep = 'MSGTGFSKLKANWKKLSRKLFSTDGSDSLNNATSKNNGKKSKGKQSDSSVSEIDQEILALIAENDFLVNGKVSLVNLGRIKEKLGDKWPKYSDFVHEFAEKVIERRITSQDLFYRVGEDVYVFVFARLTEEEAIIKCSLIAKEIGEQIFGDAWSSDEFGASIAVTKTDGSIVFEEKSIKDSIANSLLNASTVNPKSALQSVTPETATKTLEDISAKIDALGLPADEVDSEGDPEKILHNFTEMMNGVDEIVSQFNDTTELMKLNKSTPKWESFVHESQNPGTVSVSSLSQKLDALISDTEKIYEDIQETVLPSLNLQESAQDWDDSEEGAGGPEMEFCYWPVLQPAVSSIHSYRLSAEFMLEGSIWSIEELPEELEPGTIAVLDRLLLRRAIVDLLDCRDKGLLNIIIIPVHFTTLNISSLRQAYVRICAGIPKDLRKLIMWEIIDSGAGLWHSQLQTAVSAVKQFGRLVALAMESKNPRFNDLKAIGMDVVGFDCQDLGISDQEARSRLANFKKRAGQAGLRSYVFGLSSKPLLFAALKAGFDFVSGPTVAENVKQPEGVKEYHDLTGEFTLNE'
compare_pep = 'MTDTSGPNKEPPLEHWHSREPFDPYSVEKFTPEQERFYMASQWQMMWWKFRRHKVAVISGIVLLLFYISILISEFLAPYDLQTRHTKYIFAPPQSIHFFHEGEFVGPFVYGYKQRLNRETLKREYEVNTNKPYPLRFFCAGDEYEFWGLFKANLHLVCPDDKKLASFFWLGTDRLGRDILSRIIYGARISLSIGLIGIAISFTLGIILGGISGYYGGWVDNVIQRTIELLKSLPQLPLWLALSAALPVTWSPILVFFGLTVILGLLDWPGLARAVRSKFLSLREEDFCTAAQLMGAKPRRIIGRHLLPSFASHLIASASLAIPSMILGETALSFLGLGLRPPITSWGVLLNETTNINAVATTPWLIYPVLPVIIVVLAFNFLGDGLRDAADPYK'

G1 = pep_to_graph(first_pep,
                  shift=0)

best = {'pos': 0,
        'score': 0,
        'seq' : ''}

for i in range(1-len(first_pep), len(first_pep)):
    print(i)
    G2 = pep_to_graph(compare_pep,
                      shift=i)
    g, s = getMCS(G1[0], G2[0])

    if len(s) > best['score']:
        best = {'pos': i,
                'score': len(s),
                'seq': s}

G2 = pep_to_graph(compare_pep,
                  shift=best['pos'])

g, s = getMCS(G1[0], G2[0])
print(s)

plt.figure(3)
nx.draw(g, with_labels=True)
plt.show()






# draw both
# plt.figure(1)
# nx.draw_networkx_nodes(G1[0], G1[1])
# nx.draw_networkx_labels(G1[0], G1[1])
# nx.draw_networkx_edges(G1[0], G1[1], arrows=True)

# plt.figure(2)
# nx.draw_networkx_nodes(G2[0], G2[1])
# nx.draw_networkx_labels(G2[0], G2[1])
# nx.draw_networkx_edges(G2[0], G2[1], arrows=True)