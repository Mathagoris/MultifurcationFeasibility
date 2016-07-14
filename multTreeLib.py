# Multifurcation Tree Problem

# Created by Matt Dohlen, Chen Pekker, Gabriel Quiroz
# June 2016

# This class is a wrapper for biopython phylo trees. We chose to implement
# it this way in order to hide the implementation details and better illustrate
# that multifurcating trees can be binarized while maintaining (in)feasibility
from Bio import Phylo
import networkx as nx
import re


class Tree:
    def __init__(self, tree_file):
        self.tree = Phylo.read(tree_file, 'newick')
        self.LEG = self.generate_LEG()

    def loci(self, terminal):
        # gene, locus, (optional) individual
        gli = terminal.name.split('_')
        return gli[0] + "_" + gli[1]

    # Assumes there is not a node with 1 child
    def is_multifurcating(self):
        return not self.tree.is_bifurcating()

    def draw(self):
        Phylo.draw_ascii(self.tree)

    def is_feasible(self):
        return True

    def generate_LEG(self):
        return nx.Graph()

    def get_terminals_w_path(self, from_list):
        to_list = list(set(self.tree.get_terminals) - set(from_list))
        has_path = set()
        for from_terminal in from_list:
            for to_terminal in to_list:
                if self.loci(from_terminal) == self.loci(to_terminal):
                    has_path.add(self.loci(from_terminal))
        return has_path

    def expand(self, partition, clade):
        # Biopython's tree class does not let you replace a clade with another so it is
        # necessary to delete everything below the clade then split it out to create a
        # new tree
        while not clade.is_terminal():
            for terminal in clade.get_terminals():
                self.tree.prune(terminal)
        self.connect(partition, clade)

    def connect(self, partition, clade):
        if len(partition) == 1:
            self.sub_expand(partition[0])
        else:
            clade.split()
            self.sub_expand(partition[0], clade[0])
            self.connect(partition[1:], clade[1])

    def sub_expand(self, group, clade):
        if len(group) == 1:
            self.reconstruct_tree(group[0], clade)
        else:
            clade.split()
            self.reconstruct_tree(group[0], clade[0])
            self.sub_expand(group[1:], clade[1])

    def reconstruct_tree(self, subtree, clade):
        if subtree.is_terminal():
            clade.name = subtree.name
        else:
            clade.name = subtree.name
            clade.split(len(subtree))
            for i in xrange(0, len(subtree)):
                self.reconstruct_tree(subtree[i], clade[i])

    def binarize(self):
        self.binarize_rec(self.tree.root)

    def binarize_rec(self, clade):
        if clade.is_terminal():
            return
        if len(clade) > 2:
            partition = {}
            for child in clade:
                paths_on_parent_edge = []
                if child.is_terminal():
                    paths_on_parent_edge = self.get_terminals_w_path([child])
                else:
                    paths_on_parent_edge = self.get_terminals_w_path(child.get_terminals())
                if len(paths_on_parent_edge) == 0:
                    if 'no_path' in partition:
                        partition['no_path'].append(child)
                    else:
                        partition['no_path'] = [child]
                else:
                    # Arbitrarily choose the first loci on the parent edge because all the paths
                    # are in the same connected component regardless
                    cc = tuple(nx.node_connected_component(self.LEG, self.loci(paths_on_parent_edge[0])))
                    if cc in partition:
                        partition[cc].append(child)
                    else:
                        partition[cc] = [child]
            no_path = []
            if 'no_path' in partition:
                no_path = partition['no_path']
                del partition['no_path']
            if len(partition) == 0:
                self.expand([no_path], clade)
            else:
                partition[partition.keys()[0]] + no_path
                partition_list = []
                for key in partition.keys():
                    partition_list.append(partition[key])
                self.expand(partition_list, clade)

        for child in clade:
            self.binarize_rec(child)

