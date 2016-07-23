# MultifurcationFeasibility

This project is an attempted implementation of a method for binarizing a multifurcating tree outlined in a paper written by Professor Ran's summer 2016 research group. In its current state the implementation does not work in every case. Originally the algorithm was breaking when the tree's root was multifurcating. In an attempt to fix this problem a handle was added at the top of the tree to ensure that the root was not multifurcating. Although this stopped the algorithm from crashing it still fails to binarize trees with a multifurcating root in the adds-edge-to-leg.nwk file.

Due to a lack of time left in the summer it was not possible to properly debug this to find out what was going wrong. One theory about what might be happening is: after a node is expanded, the algorithm might not be properly traversing and binarizing that nodes children. This could be due to binarize_rec keeping track of the old node that has been expanded and replaced by a new tree when it should be using the root of the new tree install.

Small changes were also made to the Rasmus tree library to make replacing a subtree with a different tree possible.