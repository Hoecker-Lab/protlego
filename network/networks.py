from protlego.database.data import Result, Hit
from protlego.definitions import logger
from graph_tool.all import *
from protlego.database.data import fetch_id
from protlego.network.colors import *
import os
from protlego.builder.chimera import Chimera
from protlego.builder.chimera import get_SCOP_domain
from protlego.builder.builder import Builder
import numpy as np
from typing import Tuple

import logging
logger = logging.getLogger('protlego')

class Network:
    """
    Graph constructor and visualizer
    Examples
    --------
    hits = fetch_subspace(args)
    g = Network(hits)
    g.create_graph()
    g.plot_graph(labels=['ids','domains','folds'])
    g.get_fragments()
    g.view_vertex(v1)
    g.view_edge(v1,v2)
    g.view_component(n)
    """

    def __init__(self, hits: Result) -> None:
        self.hits = hits
        self.graph = None
        self.comp = None
        self.degrees = None

    def _add_new_vertex_universe(self, graph: Graph, o: Hit, boolean: int) -> Vertex:

        v = graph.add_vertex()
        if boolean == 0:
            graph.vp.cluster[v] = o.q_cluster
            graph.vp.domain[v] = o.query
            graph.vp.scopclass[v] = o.q_sufam_id.split('.')[0]
            graph.vp.fold[v] = o.q_fold_id
            graph.vp.sufam[v] = o.q_sufam_id
            graph.vp.fam[v] = o.q_scop_id
            graph.vp.start[v].append(o.q_start)
            graph.vp.end[v].append(o.q_end)
        else:
            graph.vp.cluster[v] = o.s_cluster
            graph.vp.domain[v] = o.sbjct
            graph.vp.scopclass[v] = o.s_sufam_id.split('.')[0]
            graph.vp.fold[v] = o.s_fold_id
            graph.vp.sufam[v] = o.s_sufam_id
            graph.vp.fam[v] = o.s_scop_id
            graph.vp.start[v].append(o.s_start)
            graph.vp.end[v].append(o.s_end)
        return v

    def create_network(self) -> Graph:
        """
        Creates a network based on the hits. It draws a node for every unique cluster
        and links two cluster when they have a fragment (hit) in common.
        :return: A Graph object
        """

        graph = Graph(directed=False)
        # Nodes properties
        graph.vp.cluster = graph.new_vertex_property("string")
        graph.vp.domain = graph.new_vertex_property("string")
        graph.vp.scopclass = graph.new_vertex_property("string")
        graph.vp.fold = graph.new_vertex_property("string")
        graph.vp.sufam = graph.new_vertex_property("string")
        graph.vp.fam = graph.new_vertex_property("string")
        graph.vp.start = graph.new_vertex_property("vector<int>")
        graph.vp.end = graph.new_vertex_property("vector<int>")

        # Edge properties
        graph.ep.prob = graph.new_edge_property("int")
        graph.ep.id = graph.new_edge_property("int")
        graph.ep.ident = graph.new_edge_property("int")
        graph.ep.no = graph.new_edge_property("int")
        graph.ep.cols = graph.new_edge_property("int")
        graph.ep.rmsd = graph.new_edge_property("float")

        for o in self.hits.hits:
            # vertex 1
            if o.query == o.sbjct: continue  # do not include links to themselves
            vertices_in_v1 = find_vertex(graph, graph.vp.cluster, o.q_cluster)
            if vertices_in_v1:
                v1 = vertices_in_v1[0]
                graph.vp.start[v1].append(o.q_start)
                graph.vp.end[v1].append(o.q_end)
            else:
                v1 = self._add_new_vertex_universe(graph, o, 0)
            # Vertex2
            vertices_in_v2 = find_vertex(graph, graph.vp.cluster, o.s_cluster)
            if vertices_in_v2:
                v2 = vertices_in_v2[0]
                graph.vp.start[v2].append(o.s_start)
                graph.vp.end[v2].append(o.s_end)
            else:
                v2 = self._add_new_vertex_universe(graph, o, 1)

            # Add edge
            ae = graph.edge(v1, v2)
            if not ae:
                ae = graph.add_edge(v1, v2)
                graph.ep.prob[ae] = o.prob
                graph.ep.id[ae] = o.id
                graph.ep.ident[ae] = o.ident
                graph.ep.cols[ae] = o.cols
                graph.ep.no[ae] = o.no
                graph.ep.rmsd[ae] = o.rmsd_tm_pair

        # save the position
        pos = sfdp_layout(graph)
        graph.vp.pos = graph.new_vertex_property("vector<double>")
        graph.vp.pos = pos
        self.graph = graph
        return graph

    def plot_graph(self, graph, fill, **keyword_parameters):
        """
        This function plots the graph
        with customized labels
        :param graph: A graph-tool object to plot
        :param fill: To choose between ['fam','sufam','fold','class']
        :param keyword_parameters: labels, output (filename). labels can be: "domain","fam","sufam","fold","scopclass"
        :return: A plot of the computed with the customized colors
        """
        if fill == 'class':
            graph.vertex_properties['class_color'] = graph.new_vertex_property("string")
            for v in graph.vertices():
                graph.vertex_properties['class_color'][v] = fillcolors['class_color'][graph.vp.scopclass[v]]
        if fill == 'fold' or fill == 'sufam' or fill == 'fam':
            graph.vertex_properties[f'{fill}_color'] = graph.new_vertex_property("string")
            for v in graph.vertices():
                graph.vertex_properties[f'{fill}_color'][v] = fillcolors[f'{fill}_color'][
                    graph.vertex_properties[fill][v]]
        if not keyword_parameters:
            graph_draw(graph, graph.vp.pos, output_size=(1500, 1500),
                       vertex_fill_color=graph.vertex_properties[f"{fill}_color"])
        if 'labels' in keyword_parameters:
            graph_draw(graph, graph.vp.pos, output_size=(1500, 1500),
                       vertex_fill_color=graph.vertex_properties[f"{fill}_color"],
                       vertex_text=graph.vertex_properties[keyword_parameters['labels']])
        if 'output' in keyword_parameters:
            graph_draw(graph, graph.vp.pos, output_size=(1500, 1500),
                       vertex_fill_color=graph.vertex_properties[f"{fill}_color"],
                       vertex_font_size=8, output=keyword_parameters['output'])

    @property
    def fragments(self):
        """
        This function creates and prints out the number of fragments in the graph
        :rtype: integer, number of fragments
        """
        if not self.graph:
            logger.info(" You need to create a network first before computing its sizes."
                        " Call create_network(). Exiting...")
            return
        self.comp, hist = label_components(self.graph)
        self.numFrags = max(self.comp.a) + 1
        logger.info("There are ", self.numFrags, " fragments")
        return self.comp

    def vertex_of_fragment(self, frag: int) -> list:
        if not self.comp:
            _ = self.fragments()
        frag = [i for i, x in enumerate(self.comp.a) if x == frag]
        return frag

    def get_representative_domain(self, frag: int) -> int:
        """
        Gets the most representative domain of a fragment (or graph component)
        :param frag: the index of the fragment
        :return: A vertex object
        """
        if not self.comp:
            _ = self.fragments()
        degrees=[]
        for vertex in self.graph.vertices():
            v1 = self.graph.vertex(vertex)
            degrees.append(v1.out_degree())
        degree = [(x, i) for i, x in enumerate(degrees) if self.comp.a[i] == frag]
        most_connected_vertices = [x[1] for x in degree if x[0] == max(degree)[0]]
        # if there are several vertices equally connected, select that one
        # whose length is the closest to the average.
        if len(most_connected_vertices) > 1:
            target = self.sizes[frag][0]
            verts_len = []
            for vertex in most_connected_vertices:
                start = int(round(np.mean(self.graph.vp.start[vertex])))
                end = int(round(np.mean(self.graph.vp.end[vertex])))
                length = int(end) - int(start)
                verts_len.append(length)
            vertice = most_connected_vertices[min(range(len(verts_len)), key=lambda i: abs(verts_len[i] - target))]
        else:
            vertice = most_connected_vertices[0]
        return vertice

    @property
    def sizes(self):
        """
        Computes the average size for each fragment.
        :return: list of tuples with the average and the stds of sizes (float)
        """
        if not self.comp:
            _ = self.fragments
        sizes = []
        for fragment_index in range(max(self.comp.a) + 1):
            frag_sizes = []
            frag = [i for i, x in enumerate(self.comp.a) if x == fragment_index]
            vfilt = self.graph.new_vertex_property('bool')
            for i in frag:
                vfilt[i] = True
            sub = GraphView(self.graph, vfilt)
            sizes.append((np.mean([self.graph.ep.cols[edge] for edge in sub.edges()]),
                          np.std([self.graph.ep.cols[edge] for edge in sub.edges()])))
        return sizes

    def show_vertex(self, vertex: Graph.vertex) -> Chimera:
        """
        Shows the protein that corresponds to that specific vertex with the
        fragment colored in red
        :param vertex: A Graph.vertex object. The domain to be shown,
        :return: A Chimera object with an internal representation of the fragment
        """
        graph = self.graph
        domain = graph.vp.domain[vertex]
        start = int(round(np.mean(graph.vp.start[vertex])))
        end = int(round(np.mean(graph.vp.end[vertex])))
        domain_path = get_SCOP_domain(domain)
        mol = Chimera(filename=domain_path, validateElements=False)
        mol.renumberResidues()
        mol.reps.add(sel='protein', style='NewCartoon', color=8)
        mol.reps.add(sel=f"protein and resid '{start}' to '{end}'", style='NewCartoon', color=1)
        mol.view(name=domain)
        return mol

    def show_edge(self, edge: Graph.edge) -> Tuple[Chimera, Chimera]:
        """
        Shows an alignment of the two domains that the edge links with their respective
        fragment they have in common colored in red
        :param edge: A graph tool edge object
        :return: The query and subject chimera objects. It also opens a VMD window with the superimposition.
        """

        graph = self.graph
        id = graph.ep.id[edge]
        hit = fetch_id(id)
        a = Builder(hit)
        aln = a.get_alignment(hit.query, hit.no)
        a._get_pairs(aln)
        a.superimpose_structures(aln) # made with global alignment
        query_mol = a.qaPDB[0].copy()
        sbjct_mol = a.saPDB[0].copy()

        qstart = query_mol.get("resid", sel=f"protein and name CA and same residue as index {a.global_qpairs[0][0]}")[0]
        qend = query_mol.get("resid", sel=f"protein and name CA and same residue as index {a.global_qpairs[0][-1]}")[0]
        sstart = sbjct_mol.get("resid", sel=f"protein and name CA and same residue as index {a.global_spairs[0][0]}")[0]
        send = sbjct_mol.get("resid", sel=f"protein and name CA and same residue as index {a.global_spairs[0][-1]}")[0]

        query_mol.reps.add(sel='protein', style='NewCartoon', color=8)
        query_mol.reps.add(sel=f"protein and resid '{qstart}' to '{qend}'", style='NewCartoon', color=1)
        sbjct_mol.reps.add(sel='protein', style='NewCartoon', color=8)
        sbjct_mol.reps.add(sel=f"protein and resid '{sstart}' to '{send}'", style='NewCartoon', color=1)
        query_mol.view()
        sbjct_mol.view()

        return query_mol, sbjct_mol

    def show_component(self, fragment: int):
        """
        Aligns the fragment(s) of all domains belonging to a component, provided the component
        has less than 50 vertices.
        Note that all domains in the same component present a fragment
        in common in different protein environments
        :param fragment:
        :return:
        """
        graph = self.graph
        if not self.comp:
            comp, _ = label_components(graph)
        else:
            comp = self.comp
        vfilt = graph.new_vertex_property('bool')
        frag = [i for i, x in enumerate(comp.a) if x == fragment]
        for i in frag:
            vfilt[i] = True
        sub = GraphView(graph, vfilt)
        if sub.num_vertices() > 50:
            raise RuntimeError("The component to visualize is too large")
        # TODO
