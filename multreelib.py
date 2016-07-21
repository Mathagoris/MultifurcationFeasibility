# Multifurcation Tree Problem

# Created by Matt Dohlen
# June 2016

# This class is a wrapper for the Rasmus CompBio tree implementation located at
# https://github.com/mdrasmus/compbio. I chose to implement it this way in order
# to hide the implementation details and better illustrate that multifurcating
# trees can be binarized while maintaining (in)feasibility

# A lot of the code for creating the LEG was repurposed for this class from Prof
# Wu's plctlib. get_conflicts and annotate have not yet been tested
from rasmus import treelib
import networkx as nx
import collections

def parse_gene(gene, mapping='sli_'):
    if mapping == 'sli':
        species, locus, ind = gene.split('-')  # leaf format = "species-locus-ind"
    elif mapping == 'sil':
        species, ind, locus = gene.split('-')  # leaf format = "species-ind-locus"
    elif mapping == 'sli_':
        species, locus, ind = gene.split('_')  # leaf format = "species_locus_ind"
    elif mapping == 'sil_':
        species, ind, locus = gene.split('_')  # leaf format = "species_ind_locus"
    else:
        raise Exception("mapping not supported: %s" % mapping)

    return (species, locus, ind)

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
        # print self.tree.nodes

    def draw_leg(self):
        # nx.draw(self.LEG)
        print "Connected Components of LEG:\n" + str(list(nx.connected_components(self.leg)))

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
            gene = parse_gene(leaf.name, self.mapping)

            label = gene[:2]
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

    def get_paths_out(self, from_leaves):
        to_leaves = list(set(self.tree.leaves()) - set(from_leaves))
        has_path = set()
        for from_leaf in from_leaves:
            for to_leaf in to_leaves:
                if parse_gene(from_leaf.name, self.mapping)[:2] == \
                                parse_gene(to_leaf.name, self.mapping)[:2]:
                    has_path.add(from_leaf)
        return has_path

    def expand(self, partition, node):
        connecting_tree = treelib.Tree()
        connecting_tree.make_root(name=node.name)
        self.connect(partition, connecting_tree, connecting_tree.root)
        if node.parent is None:
            self.tree = connecting_tree
        else:
            self.tree.replace_tree(node, connecting_tree)

    def connect(self, partition, connecting_tree, node):
        if len(partition) == 1:
            self.sub_expand(partition[0], connecting_tree, node)
        else:
            left = treelib.TreeNode(name=self.tree.new_name())
            right = treelib.TreeNode(name=self.tree.new_name())
            connecting_tree.add_child(node, left)
            connecting_tree.add_child(node, right)
            self.sub_expand(partition[0], connecting_tree, left)
            self.connect(partition[1:], connecting_tree, right)

    def sub_expand(self, group, connecting_tree, node):
        if len(group) == 1:
            connecting_tree.replace_tree(node, group[0])
        else:
            connecting_tree.add_tree(node, group[0])
            right = treelib.TreeNode(name=self.tree.new_name())
            right = connecting_tree.add_child(node, right)
            self.sub_expand(group[1:], connecting_tree, right)

    def binarize(self):
        self.binarize_rec(self.tree.root)
        # Paths may have been generated within connected components of LEG
        # when binarizing so it is necessary to regenerate the LEG
        # NOTE: For some reason if the tree is not copied here the create_plct
        # method fails when trying to add labels to the data dictionary
        self.tree = self.tree.copy()
        self.leg = self.create_leg()

    def binarize_rec(self, node):
        if node.is_leaf():
            return
        if len(node.children) > 2:
            partition = collections.defaultdict(list)
            for child in node:
                paths_on_parent_edge = []
                # May be possible to use node "labels" here from creating PLCT
                if child.is_leaf():
                    paths_on_parent_edge = self.get_paths_out([child])
                else:
                    paths_on_parent_edge = self.get_paths_out(child.leaves())
                if len(paths_on_parent_edge) == 0:
                    partition['no_path'].append(treelib.subtree(self.tree, child))
                else:
                    # Arbitrarily choose the first loci on the parent edge because all the loci with
                    # paths on parent edge are in the same connected component regardless
                    cc = nx.node_connected_component(self.leg,
                                        parse_gene(paths_on_parent_edge.pop().name, self.mapping)[:2])
                    if len(cc) == 1:
                        cc = cc.pop()
                    else:
                        cc = tuple(cc)
                    partition[cc].append(treelib.subtree(self.tree, child))
            no_path = []
            if 'no_path' in partition:
                no_path = partition['no_path']
                del partition['no_path']
            if len(partition) == 0:
                self.expand([no_path], node)
            else:
                # arbitrarily place children with no path in any partition
                partition[partition.keys()[0]] + no_path
                partition_list = []
                for key in partition.keys():
                    partition_list.append(partition[key])
                self.expand(partition_list, node)

        for child in node:
            self.binarize_rec(child)
