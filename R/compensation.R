#' Compensate data based on spillover matrix save in FCS files
#'
#' Applies compensation matrices for each channel and saves compensated values into "comp" assay in Seurat Object. Uses [flowCore::compensate()] function.
#'
#' @param fcs_fs FlowSet containing all fcs files and compensation matrices. Must be the same FlowSet that was used for [create_seurat()].
#' @param seu Seurat object used for [create_seurat()].
#' @param comp_matrix If multiple compensation matrices are present in fcs file, indicate which one to use (check with names(spillover(fcs_fs[[x]]))).
#' @return Seurat object with new assay "comp", where compensation values were saved in slot "counts".
#' @export
#' @examples
#' # Check if there are multiple spillover matrices saved in FCS files
#' names(spillover(fcs_fs[[1]]))
#' # Pick the correct spillover matrix for compensation
#' seu = compensate_data(fcs_fs, seu, comp_matrix = 3)
compensate_data <- function(fcs_fs,
                            seu,
                            comp_matrix = 1) {
    
  # check some parameters
  if(is.null(spillover(fcs_fs[[1]])[[comp_matrix]]))
    stop("Spillover matrix empty!", call. = FALSE)
  # get list of compensation matrices (one for each sample) using matrices in 3rd column
  comp <- fsApply(fcs_fs, function(x) spillover(x)[[comp_matrix]], simplify = FALSE)
  # compensate and save in new flowSet object
  fs_comp <- compensate(fcs_fs, comp)
  # generate Seurat object
  seu_comp <- create_seurat(fs_comp, data.frame(seu@misc))
  # add compensated matrix to new assay
  seu[["comp"]] <- CreateAssayObject(counts = GetAssayData(seu_comp),
                                     min.cells = 0,
                                     min.features = 0)
  # make compensated data as default assay
  DefaultAssay(seu) <- "comp"
  return(seu)
}
