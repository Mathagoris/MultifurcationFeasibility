from sys import argv
from multreelib import Tree


def main(tree_file):
    print "Reading tree in file: " + tree_file
    tree = Tree(tree_file)
    is_mult = tree.is_multifurcating()
    print "Tree type: " + "Multifurcating" if is_mult else "Bifurcating"
    initial_feas = tree.is_feasible()
    print "Feasibility: " + "Feasible" if initial_feas else "Infeasible"
    tree.draw_tree()
    tree.draw_leg()
    if is_mult:
        print "Binarizing multifurcating tree..."
        tree.binarize()
        if tree.is_multifurcating():
            print "Failed to binarize"
            tree.draw_tree()
        else:
            print "The tree has been binarized"
            if initial_feas == tree.is_feasible():
                print "The tree remained " + "feasible" if initial_feas else "infeasible"
            else:
                print "Something went wrong!"
                print "The tree is now " + "infeasible" if initial_feas else "feasible"
            tree.draw_tree()
            tree.draw_leg()
    else:
        print "Tree does not need to be binarized"


if __name__ == '__main__':
    if len(argv) != 2:
        print "Must include path to tree file"
    else:
        main(argv[1])
