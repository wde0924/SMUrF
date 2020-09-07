# subroutine to prepare vegetated fraction from VCF, MOD44B
# by DW, 10/29/2019 

# last modify DW, 03/04/2020
prep.mod44b <- function(vcf.path, vcf.pattern, yr, reg.ext) {

    vcf.files <- list.files(vcf.path, vcf.pattern, full.names = T)
    vcf.files <- vcf.files[grepl(yr, vcf.files)]
    if (length(vcf.files) == 0) stop('prep.mod44b(): NO VCF files found, please check\n')

    # 0-100 for veg fractions, 200 for water bodies
    tree.file <- vcf.files[grepl('Tree', vcf.files)]
    nonveg.file <- vcf.files[grepl('NonVeg', vcf.files)]

    # read 250m tree and non-tree veg fractions
    if (length(tree.file) > 1) stop('prep.mod44b(): more than 1 TREE fraction file, please check "vcf.pattern"\n')
    if (length(nonveg.file) > 1) stop('prep.mod44b(): more than 1 NON-veg fraction file, please check "vcf.pattern"\n')
    tree.rt <- crop(raster(tree.file), reg.ext)
    nonveg.rt <- crop(raster(nonveg.file), reg.ext)

    # remove missing value of 200 
    tree.rt[tree.rt == 200] <- NA 
    nonveg.rt[nonveg.rt == 200] <- NA   # total non-vegetated fractions

    veg.rt <- 100 - nonveg.rt           # total vegetated fractions
    nontree.veg.rt <- veg.rt - tree.rt  # non-tree vegetated fractions

    frac.stk <- stack(veg.rt, nonveg.rt, tree.rt, nontree.veg.rt)
    names(frac.stk) <- c('veg_frac', 'nonveg_frac', 'tree_frac', 'nontree_veg_frac')

    return(frac.stk)
}  

# end of script