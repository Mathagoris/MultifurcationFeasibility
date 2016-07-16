# Multifurcation Tree Problem

# Created by Matt Dohlen
# June 2016

# This class is a wrapper for Prof Wu's tree implementation. I chose to implement
# it this way in order to hide the implementation details and better illustrate
# that multifurcating trees can be binarized while maintaining (in)feasibility
from rasmus import treelib
import networkx as nx
import collections


class Tree(object):
    def __init__(self, tree_file, mapping='sli_'):
        self.tree = treelib.read_newick(tree_file)
        self.labeled = False
        self.mapping = mapping
        self.leg = self.create_leg()

    # Assumes there is not a node with 1 child
    def is_multifurcating(self):
        return not treelib.is_binary(self.tree)

    def draw_tree(self):
        treelib.draw_tree(self.tree)

    def draw_LEG(self):
        # nx.draw(self.LEG)
        print self.leg.nodes()

    def is_feasible(self):
        for cc in nx.connected_components(self.leg):
            loci_dct = collections.defaultdict(set)
            for label in cc:
                species, locus = label
                loci_dct[species].add(locus)
            for species in loci_dct.keys():
                if len(loci_dct[species]) >= 2:
                    return False
        return True

    def group_leaves(self):
        # collect leaves based on species and locus
        groupings = collections.defaultdict(list)
        for leaf in self.tree.leaves():
            if self.mapping == 'sli':
                species, locus, ind = leaf.name.split('-')  # leaf format = "species-locus-ind"
            elif self.mapping == 'sil':
                species, ind, locus = leaf.name.split('-')  # leaf format = "species-ind-locus"
            elif self.mapping == 'sli_':
                species, locus, ind = leaf.name.split('_')  # leaf format = "species_locus_ind"
            elif self.mapping == 'sil_':
                species, ind, locus = leaf.name.split('_')  # leaf format = "species_ind_locus"
            else:
                raise Exception("mapping not supported: %s" % self.mapping)

            label = (species, locus)
            groupings[label].append(leaf)

        return groupings

    def create_plct(self, groupings, new_copy=False):
        tree = self.tree
        self.labeled = not new_copy
        if new_copy:
            tree = tree.copy()

        for node in tree:
            node.data["labels"] = set()

        for label, genes in groupings.iteritems():
            lca = treelib.lca(genes)
            for leaf in genes:
                # follows each leaf up and labels branch by storing color in node data
                travel = leaf
                while travel != lca:
                    # add label to branch
                    travel.data["labels"].add(label)
                    travel = travel.parent
        return tree

    def create_leg(self):
        """Creates leg from plct and groupings."""
        leg = nx.Graph()
        groupings = self.group_leaves()
        plct = self.create_plct(groupings)
        leg.add_nodes_from(groupings.keys())  # nodes = (species, locus)
        for node in plct:
            labels = list(node.data["labels"])  # convert label set to label list
            for i in xrange(len(labels)):
                for j in xrange(i + 1, len(labels)):
                    assert labels[i] in leg and labels[j] in leg
                    leg.add_edge(labels[i], labels[j])
        return leg

    def get_conflicts(self):
        """Find irreconcilable connected components of leg."""
        conflicts = set()  # connected components with conflict
        for cc in nx.connected_components(self.leg):
            # key = species, val = set of loci in species for this cc
            loci_dct = collections.defaultdict(set)

            for label in cc:
                species, locus = label
                loci_dct[species].add(locus)

            for sp, loci in loci_dct.iteritems():
                # conflict if a species has more than one loci in this cc
                if len(loci) >= 2:
                    conflicts.add(tuple(cc))
                    break
        return conflicts

    def annotate(self):
        """Annotate tree."""
        # reconcilable => no labels on this branch are pairwise irreconcilable
        # reconcilable_cc => no labels on this branch are part of irreconcilable cc of leg
        if not self.labeled:
            raise Exception("Cannot annotate because tre is unlabeled.")
        conflicts = self.get_conflicts()
        conflicting_labels = set()
        for cc in conflicts:
            conflicting_labels.update(cc)

        # self.tree must be plct
        for node in self.tree:
            labels = node.data["labels"]
            if not labels:
                continue

            node.data["reconcilable_cc"] = True
            loci_dct = collections.defaultdict(set)
            for label in labels:
                if label in conflicting_labels:
                    node.data["reconcilable_cc"] = False

                species, locus = label
                loci_dct[species].add(locus)

            if node.is_leaf():  # a leaf always has a single label
                continue  # so there are no pairs to consider
            node.data["reconcilable"] = True
            for sp, loci in loci_dct.iteritems():
                # conflict if a species has more than one loci
                if len(loci) >= 2:
                    node.data["reconcilable"] = False

    def binarize(self):
