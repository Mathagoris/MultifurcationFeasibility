# Multifurcation Tree Problem

# Created by Matt Dohlen, Chen Pekker, Gabriel Quiroz
# June 2016

# This class is a wrapper for biopython phylo trees. We chose to implement
# it this way in order to hide the implementation details and better illustrate
# that multifurcating trees can be binarized while maintaining (in)feasibility
from Bio import Phylo
import networkx as nx


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

    def draw_tree(self):
        Phylo.draw_ascii(self.tree)

    def draw_LEG(self):
        # nx.draw(self.LEG)
        print self.LEG.nodes()

    def is_feasible(self):
        for cc in nx.connected_components(self.LEG):
            for i in xrange(0, len(cc)):
                for j in xrange(i+1, len(cc)):
                    if list(cc)[i].split('_')[0] == list(cc)[j].split('_')[0]:
                        return False
        return True

    def generate_LEG(self):
        LEG = nx.Graph()
        LEG.add_nodes_from([terminal.name for terminal in self.tree.get_terminals()])
        group_paths = {}
        for group in self.group_terminals():
            group_locus = self.loci(group[0])
            for terminal in group[1:]:
                path = self.tree.trace(group[0], terminal)
                if group_locus in group_paths:
                    group_paths[group_locus].append(set(path))
                else:
                    group_paths[group_locus] = [set(path)]
        group_loci = group_paths.keys()
        for i in xrange(0, len(group_loci)):
            loci_1 = group_loci[i]
            for j in xrange(i+1, len(group_loci)):
                loci_2 = group_loci[j]
                if not self.paths_are_edge_disjoint(group_paths[loci_1], group_paths[loci_2]):
                    LEG.add_edge(loci_1, loci_2)
        return LEG

    def group_terminals(self):
        terminals = sorted(self.tree.get_terminals(), key=lambda terminal: terminal.name)
        groups = []
        while terminals != []:
            first = terminals.pop(0)
            group = [first]
            to_remove = []
            for terminal in terminals:
                if first.name in terminal.name:
                    group.append(terminal)
                    to_remove.append(terminal)
            for terminal in to_remove:
                terminals.remove(terminal)
            groups.append(group)
        return groups


    def paths_are_edge_disjoint(self, paths_1, paths_2):
        for path_1 in paths_1:
            for path_2 in paths_2:
                if len(path_1.intersection(path_2)) > 2:
                    return False
        return True

    def get_terminals_w_path(self, from_list):
        to_list = list(set(self.tree.get_terminals()) - set(from_list))
        has_path = set()
        for from_terminal in from_list:
            for to_terminal in to_list:
                if self.loci(from_terminal) == self.loci(to_terminal):
                    has_path.add(from_terminal)
        return has_path

    def prune(self, clade):
        if clade.is_terminal():
            self.tree.prune(clade)
        else:
            for child in clade:
                self.prune(child)
            self.tree.prune(clade)

    def expand(self, partition, clade):
        # Biopython's tree class does not let you replace a clade with another so it is
        # necessary to delete everything below the clade then split it out to create a
        # new tree
        while not clade.is_terminal():
            for terminal in clade.get_terminals():
                self.tree.prune(terminal)
            clade = self.tree.get_terminals()
        self.connect(partition, clade)

    def connect(self, partition, clade):
        if len(partition) == 1:
            self.sub_expand(partition[0], clade)
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
        # Paths may have been generated within connected components of LEG
        # when binarizing so it is necessary to regenerate the LEG
        self.LEG = self.generate_LEG()

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
                    # Arbitrarily choose the first loci on the parent edge because all the loci with
                    # paths on parent edge are in the same connected component regardless
                    cc = nx.node_connected_component(self.LEG, self.loci(paths_on_parent_edge.pop()))
                    if len(cc) == 1:
                        cc = cc.pop()
                    else:
                        cc = tuple(cc)
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

