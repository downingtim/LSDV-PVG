install_missing <- function(pkgs)
{
    for (pkg in pkgs) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, dependencies = TRUE)
        }
    }
}

#
install_missing(c("ggExtra","phangorn","Biostrings","adegenet","ape","stringr","seqinr","ggtree","ggnewscale","restez","rentrez","tibble","irlba","phytools","ggplot2","treeio","Rcpp","RcppArmadillo","microseq","hierfstat"))
#install_missing(c("phangorn","Biostrings","ape","stringr","seqinr","ggtree","ggnewscale","restez","rentrez","tibble","irlba","phytools","ggplot2","treeio","Rcpp","RcppArmadillo","microseq","hierfstat"))
