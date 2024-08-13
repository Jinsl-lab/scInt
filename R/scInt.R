#' The scInt Class
#'
#' The scInt object is created from a list of single-cell single-cell RNA-seq data matrix with features as rows and cells as columns.
#' The class provides functions for data preprocessing, data integration, and visualization.
#'
#' The key slots used in the scInt object are described below.
#'
#' @slot raw.data List of raw data matrix, one per dataset (features as rows and cells as columns).
#' @slot norm.data List of normalized data matrix (features by cells).
#' @slot meta Data frame storing the information associated with each cell.
#' @slot cells All cells names.
#' @slot vargenes Identified high variable features.
#' @slot var_list Identified high variable features in each dataset.
#' @slot integrated.data Integrated single-cell RNA-seq matrix with rows as selected high variable features and columns as all cells.
#' @slot reduction List of the reduced dimensional space, one per method, e.g., pca/umap/tsne.
#' @slot int.umap UMAP reduced dimensional data frame.
#' @slot variation Learned contrastive biological variation matrix.
#' @slot projection Learned projection matrix.
#' @slot S Learned cell-cell similarity matrix.
#' @slot pre.clusters Preassigned cluster label list for each batch in the integration process.
#' @slot parameters List of parameters used throughout analysis.
#' @slot modal The data modal, such as RNA, ATAC.
#' @slot version Version of package used to create object.
#'
#' @exportClass scInt
#' @name scInt-class
#' @rdname scInt-class
#' @import Matrix
#' @importFrom Rcpp evalCpp
#' @importFrom methods setClass
#' @useDynLib scInt
scInt <- methods::setClass(
  "scInt",
  slots = c(
    raw.data = "list",
    norm.data = "list",
    meta = "data.frame",
    cells = "vector",
    vargenes = "vector",
    var_list = "list",
    integrated.data = "matrix",
    reduction = "list",
    int.umap = "data.frame",
    variation = "matrix",
    projection = "matrix",
    S = 'dgCMatrix',
    pre.clusters = 'list',
    parameters = "list",
    modal = "character",
    version = 'ANY'
  )
)

#' show method for scInt
#'
#' @param object scInt object
#' @name show
#' @docType methods
#' @rdname show-methods
#'
setMethod(
  f = "show",
  signature = "scInt",
  definition = function(object) {
    cat(
      "An object of class",
      class(object),
      "\nwith",
      length(object@raw.data),
      "datasets and\n",
      sum(sapply(data@norm.data, ncol)),
      "total cells."
    )
    invisible(x = NULL)
  }
)

#' Create a new scInt object
#'
#' This function initializes a scInt object of multiple single cell dataset, for further integration and analysis.
#' It takes a list of multiple gene-cell expression matrices as inputs.
#' @param raw.data List of raw data matrices (gene by cell).
#' @param meta The meta data of input datasets. (default NULL)
#' @param filter Whether filter genes and cells. (default TRUE)
#' @param min.features Filter the cells expressed features less than min.features. (default 300)
#' @param min.cells Filter the features expressed in cells less than min.cells. (default 3)
#' @param modal Dataset modal.
#' @param do.sparse Whether use sparse format ('dgCMatirx' class). (default TRUE)
#'
#' @return \code{scInt} object.
#'
#' @import Matrix
#' @import methods
#'
#' @export
create.scInt <- function(raw.data,
                         meta = NULL,
                         filter = TRUE, min.features = 300, min.cells = 3,
                         modal = c("single", "multi"),
                         do.sparse = TRUE) {
  object <- methods::new(
    Class = "scInt"
  )

  if (do.sparse) {
    raw.data <- lapply(raw.data, function(x) {
      as(as.matrix(x), "dgCMatrix")
    })
  }

  if (is.null(names(raw.data))) {
    names(raw.data) <- paste0("Batch_", 1:length(raw.data))
  }

  if (length(Reduce(intersect, lapply(raw.data, colnames))) > 0 && length(raw.data) > 1) {
    stop ('Cells across batches must have unique name.')
  }

  cells = Reduce('c', lapply(raw.data, colnames))
  if (!is.null(meta)) {
    if (is.data.frame(meta)) {
      if (nrow(meta) != length(cells) ||
          length(intersect(cells, rownames(meta))) != length(cells)) {
        stop('Wrong meta informations are provided.')
      } else object@meta <- meta
    } else stop('meta must be a data frame.')
  } else {
    # Initialize meta with nUMI, nGene, batchlb, and cell
    nUMI <- unlist(lapply(raw.data, function(x) Matrix::colSums(x)), use.names = FALSE)
    nGene <- unlist(lapply(raw.data, function(x) Matrix::colSums(x > 0)), use.names = FALSE)
    batchlb <- unlist(lapply(seq_along(raw.data), function(i) {
      rep(names(raw.data)[i], ncol(raw.data[[i]]))
    }), use.names = FALSE)
    cell <- cells
    meta <- data.frame(nUMI, nGene, batchlb, cell)
    rownames(meta) <- cell
    object@meta <- meta
  }

  if (filter) {
    for (i in 1:length(raw.data)) {
      rawmatrix <- raw.data[[i]]
      num.features <- Matrix::colSums(rawmatrix > 0)
      cells.use <- names(x = num.features[which(x = num.features > min.features)])
      raw.data[[i]] <- rawmatrix[, cells.use]
      filter.cell <- setdiff(colnames(rawmatrix), cells.use)
      message(paste0(length(filter.cell), ' cells are filtered in the ', names(raw.data)[i]))
    }

    rawdata <- do.call(cbind, raw.data)
    object@meta <- object@meta[colnames(rawdata), ]
    # filter genes
    allgenes <- rownames(rawdata)
    if (min.cells > 0) {
      num.cells <- Matrix::rowSums(rawdata > 0)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      for (i in 1:length(raw.data)){
        raw.data[[i]] <- raw.data[[i]][genes.use, ]
      }
      filter.gene <- setdiff(allgenes, genes.use)
      message(paste0(length(filter.gene), ' genes are filtered'))
    } else {
      message("Genes are no filtered")
    }
  }

  modal = modal[1]
  object@modal = modal
  if (modal == "single") {
    if (length(unique(sapply(raw.data, nrow))) != 1 ||
        length(Reduce(intersect, lapply(raw.data, rownames))) != nrow(raw.data[[1]])) {
      stop ('Features must be same in the single modal datasets integration.')
    }
  }

  object@cells <- Reduce('c', lapply(raw.data, colnames))

  object@raw.data <- raw.data
  object@norm.data <- raw.data
  object@version <- packageVersion("scInt")

  return(object)
}

#' Normalize scInt object
#'
#' This function normalizes the raw dataset of scInt object.
#'
#' @param object \code{scInt} object.
#' @param scale.factor The scale factor used for each cell. (default 10000)
#' @param do.log Whether do log1p transformation. (default TRUE)
#'
#' @return \code{scInt} object with norm.data slot as normalized data.
#'
#' @import Matrix
#'
#' @export
normalize.Data <- function(object, 
                           scale.factor = 10000, 
                           do.log = TRUE) {
  norm.data <- lapply(X = object@raw.data,
                      FUN = function(x) {
                        library.size <- Matrix::colSums(x)
                        expr_mat <- Matrix::t(Matrix::t(x) / library.size) * scale.factor
                        if (do.log) {
                          expr_mat <-log1p(expr_mat)
                        }
                        expr_mat
                      })
  object@norm.data <- norm.data
  return(object)
}

#' Identify the high variable features and pre-clustered cell identities for each dataset.
#'
#' @param object \code{scInt} object.
#' @param dims The dimension of PCA. (default 30)
#' @param nfeatures The number of selected features in each dataset.
#' @param select.var The model to select features. "union" directly combines all features of each dataset, "select" selects the features with the highest frequency across all datasets.
#' @param vargenes User-specified features. (default NULL)
#' @param resolution The resolution of Leiden algorithm. (default 0.8)
#'
#' @return \code{scInt} object with pre.clusters, var_list and vargenes slots set.
#'
#' @import Seurat
#' @importFrom dplyr `%>%`
#'
#' @export
ident.cellIdentity <- function(object, 
                               dims = 30, 
                               nfeatures = 1000, 
                               select.var = c("union", "select"), 
                               vargenes = NULL, 
                               resolution = 0.8) {
  select.var <- match.arg(select.var)
  batch_list <- object@norm.data
  var_list <- list()
  cluster_iden <- list()
  NN <- sum(sapply(batch_list, ncol))
  # clustering for each batch
  for(i in 1:length(batch_list)) {
    batch_list[[i]] <- Seurat::CreateSeuratObject(batch_list[[i]], min.cells = 0, min.features = 0)
    if (substring(batch_list[[i]]@version, 1, 1) == "5") {
      batch_list[[i]] <- Seurat::NormalizeData(batch_list[[i]])
      batch_list[[i]]@assays$RNA@layers$data <- batch_list[[i]]@assays$RNA@layers$counts
      batch_list[[i]] <- Seurat::FindVariableFeatures(batch_list[[i]], 
                                                      nfeatures = nfeatures, verbose = FALSE, clip.max = sqrt(NN))
      batch_list[[i]] <- batch_list[[i]] %>% Seurat::ScaleData(verbose = FALSE) %>% 
        Seurat::RunPCA(npcs = dims, verbose = FALSE) %>% 
        Seurat::FindNeighbors(dims = 1:dims, verbose = FALSE) %>% 
        Seurat::FindClusters(resolution = resolution, verbose = FALSE)
      var_list[[i]] <- Seurat::VariableFeatures(batch_list[[i]])
    } else if (substring(batch_list[[i]]@version, 1, 1) == "4") {
      batch_list[[i]] <- Seurat::FindVariableFeatures(batch_list[[i]], 
                                                      nfeatures = nfeatures, verbose = FALSE, clip.max = sqrt(NN))
      batch_list[[i]] <- batch_list[[i]] %>% Seurat::ScaleData(verbose = FALSE) %>% 
        Seurat::RunPCA(npcs = dims, verbose = FALSE) %>% 
        Seurat::FindNeighbors(dims = 1:dims, verbose = FALSE) %>% 
        Seurat::FindClusters(resolution = resolution, verbose = FALSE)
      var_list[[i]] <- batch_list[[i]]@assays[["RNA"]]@var.features
    }
    cluster_iden[[i]] <- paste0(as.character(i), "_", as.character(batch_list[[i]]@active.ident))
  }
  names(var_list) <- names(batch_list)
  if (!is.null(vargenes) == TRUE) {
    vargenes <- vargenes
  } else if (select.var == "union") {
    vargenes <- Reduce(union, lapply(batch_list, Seurat::VariableFeatures))
  } else if (select.var == "select") {
    vargenes <- Seurat::SelectIntegrationFeatures(batch_list, nfeatures = nfeatures)
  }

  vargenes <- rownames(object@norm.data[[1]])[which(rownames(batch_list[[1]]) %in% vargenes)]
  var_list <- lapply(var_list,
                     function(x) rownames(object@norm.data[[1]])[which(rownames(batch_list[[1]]) %in% x)])
  names(cluster_iden) <- names(batch_list)

  object@pre.clusters <- cluster_iden
  object@vargenes <- vargenes
  object@var_list <- var_list
  return(object)
}

#' Scale function
#'
#' This function scales a matrix.
#'
#' @param data.x Input data matrix.
#' @param do.center Whether center the row values. (default TRUE)
#' @param do.scale Whether scale the row values. (default TRUE)
#' @param row.means The provided row means to center. (default NULL)
#' @param row.sds The provided row standard deviations to scale. (default NULL)
#'
#' @import Matrix
#' @import matrixStats
scale_data <- function(data.x, 
                       do.center = T, 
                       do.scale = T, 
                       row.means = NULL, 
                       row.sds = NULL) {
  if (do.center) {
    if (is.null(row.means)) {
      data_mean <- Matrix::rowMeans(data.x)
    } else {
      data_mean <- row.means
    }
    data.x <- data.x - sapply(1:ncol(data.x), function(i) data_mean)
  }
  if (do.scale) {
    if (is.null(row.sds)) {
      data_stddev <- matrixStats::rowSds(as.matrix(data.x))
    } else {
      data_stddev <- row.sds
    }
    index <- which(data_stddev > 0)
    data.x[index, ] <- data.x[index, ] / sapply(1:ncol(data.x), function(i) data_stddev[index])
  }
  data.x
}

#' Perform UMAP dimensionality reduction
#'
#' This function runs UMAP on the integrated data matrix of scInt object or a single matrix.
#'
#' @param object \code{scInt} object or a single matrix. The \code{scInt} should be already integrated.
#' @param is.meta Whether add the meta information. (default FALSE)
#'
#' @return \code{scInt} object with int.umap slot as UMAP dimensionality reduction results or a single UMAP dimensionality reduction data frame.
#'
#' @import uwot
#'
#' @export
run.umap <- function(object, 
                     n_neighbors = 30, 
                     learning_rate = 0.5, 
                     init = "laplacian", 
                     n_components = 2,
                     metric = 'cosine', 
                     fast_sgd = FALSE, 
                     n_sgd_threads = 1, 
                     min_dist = .25,
                     n_threads = 4, 
                     ret_model = TRUE,  
                     is.meta = FALSE) {
  if ("scInt" %in% class(object)) {
    data_matrix <- t(object@integrated.data)
    out_umap <- uwot::umap(data_matrix, n_neighbors = n_neighbors, learning_rate = learning_rate, init = init, n_components = n_components,
                           metric = metric, fast_sgd = fast_sgd, n_sgd_threads = n_sgd_threads,
                           min_dist = min_dist, n_threads = n_threads, ret_model = ret_model)
    umap_df<- as.data.frame(out_umap$embedding)
    colnames(umap_df) <- c('UMAP_1', 'UMAP_2')
    rownames(umap_df) <- colnames(object@integrated.data)
    if (is.meta) {
      umap_df <- cbind(umap_df, meta[rownames(umap_df), ])
    }
    object@int.umap <- umap_df
  } else if ("matrix" %in% class(object)) {
    data_matrix <- t(object)
    out_umap <- uwot::umap(data_matrix, n_neighbors = n_neighbors, learning_rate = learning_rate, init = init, n_components = n_components,
                           metric = metric, fast_sgd = fast_sgd, n_sgd_threads = n_sgd_threads,
                           min_dist = min_dist, n_threads = n_threads, ret_model = ret_model)
    umap_df<- as.data.frame(out_umap$embedding)
    colnames(umap_df) <- c('UMAP_1', 'UMAP_2')
    object <- umap_df
  }
  return(object)
}

#' Perform global integration on multiple normalized datasets
#'
#' This function integrates multiple normalized matrices using identified similarity matrix S.
#'
#' @param object \code{scInt} object.
#' @param npcs The dimension of projection matrix. (default 40)
#' @param lambda The regularization parameter of technical variation. (default 5)
#' @param select.lambda Whether to select parameter lambda. (default FALSE)
#' @param cand_lambdas The candidate lambda's. (default c(0.1, 1, 5, 10, 15, 20, 50, 100))
#'
#' @return \code{scInt} object with variation, projection and integrated.data slots set.
#'
#' @import Matrix
#' @import RSpectra
#'
#' @export
run.Integrate <- function(object, 
                          npcs = 40, 
                          lambda = 5,
                          select.lambda = FALSE,
                          cand_lambdas = c(0.1, 1, 5, 10, 15, 20, 50, 100)) {
  Z <- do.call(cbind, lapply(object@norm.data, function(x) x[object@vargenes, ]))

  N <- length(object@norm.data)
  cell.names <- lapply(object@norm.data, colnames)
  nums <- sapply(cell.names, length)

  t1 = Sys.time()
  S <- object@S
  rownames(S) = colnames(S) = Reduce('c', cell.names)
  NNN <- sapply(cell.names, function(i) length(S[i, ]@x))
  sort_N <- order(NNN, decreasing = T)
  S <- Matrix::t(S[Reduce('c', cell.names[sort_N]), Reduce('c', cell.names[sort_N])])
  uS <- upper_tri(S@x, S@p, S@i, ncol(S))
  S <- Matrix::t(sparseMatrix(i = as.vector(uS$i), j = as.vector(uS$j), x = as.vector(uS$x),
                      dims = c(ncol(S), ncol(S)), repr = "C"))
  rownames(S) = colnames(S) = Reduce('c', cell.names[sort_N])
  S <- S[Reduce('c', cell.names), Reduce('c', cell.names)]
  rownames(S) = colnames(S) = Reduce('c', cell.names)
  tcross_S <- lapply(1:N,
                     function(ns) {
                       x <- Matrix(0, nrow = nrow(S), ncol = nrow(S), sparse = T)
                       colnames(x) = rownames(x) <- Reduce('c', cell.names)
                       x[cell.names[[ns]], ] <- S[cell.names[[ns]], ]
                       num <- lapply(cell.names, function(j) Matrix::colSums(x[, j]))
                       x <- Matrix::Diagonal(x = Reduce('c', num)) - x
                       if (sum(sapply(num, sum)) > 0) {
                         x <- x/sqrt(sum(sapply(num, sum)))
                       }
                       Matrix::tcrossprod(x)
                     })
  tcross_S <- Reduce('+', tcross_S)

  Z <- scale_data(Z, do.scale = F)
  if (select.lambda) {
    message("Select optimal lambda")
    K_list <- lapply(cand_lambdas, 
                     function(lam) {
                       Matrix::Diagonal(x = Reduce("c", lapply(nums, function(num) rep(1/num, num)))) - lam * tcross_S
                     })
    d_list <- lapply(K_list, function(K) as(Z %*% Matrix::tcrossprod(K, Z), "dgCMatrix"))
    dd <- sapply(d_list, function(d) length(which(diag(d) > 0)))
    if (all(dd < npcs)) {
      ii <- 1
    } else {
      ii <- min(rev(which(dd > npcs))[1]+1, length(cand_lambdas))
    }
    optimal_lambda <- cand_lambdas[ii]
    message(paste0("Optimal lambda = ", optimal_lambda))
    d <- d_list[[ii]]
    rm(K_list, d_list)
  } else {
    K <- Matrix::Diagonal(x = Reduce("c", lapply(nums, function(num) rep(1/num, num)))) - lambda * tcross_S
    d <- as(Z %*% Matrix::tcrossprod(K, Z), "dgCMatrix")
  }
  res_cpca <- RSpectra::eigs_sym(d, k = npcs, which = "LA")[["vectors"]] %>% as.matrix()
  L = as.matrix(Matrix::crossprod(res_cpca, do.call(cbind, lapply(object@norm.data, function(x) x[object@vargenes, ]))))
  t2 = Sys.time()
  print(t2-t1)
  object@variation <- as.matrix(d)
  object@projection <- t(res_cpca)
  object@integrated.data <- L
  return(object)
}

#' Compute the across-batch cell similarity relationships
#'
#' @param object \code{scInt} object.
#' @param k The number of computed top similarities for each cell. (default 5)
#' @param dims The dimension of PCA. (default 30)
#' @param T_th The similarity threshold for filtering cells in the cPCA space. (default 0.6)
#' @param select.T Whether to select parameter T_th. (default FALSE)
#' @param cand_Ts The candidate T_th's. (default seq(0.5, 0.8, 0.05))
#'
#' @return \code{scInt} object with S slot set.
#'
#' @import Matrix
#' @import Rcpp
#' @importFrom dplyr `%>%`
#'
#' @export
compute.Similarity <- function(object, 
                               k = 5, 
                               dims = 50, 
                               T_th = 0.6,
                               select.T = FALSE,
                               cand_Ts = seq(0.5, 0.8, 0.05)) {
  N <- length(object@norm.data)
  NN <- sum(sapply(object@norm.data, ncol))
  cell.names <- lapply(object@norm.data, colnames)
  cell.nums <- sapply(cell.names, length)
  cluster_iden <- object@pre.clusters

  t1 = Sys.time()
  message("Compute the Adjusted Cross-batch Similarity")

  modis <- list()
  if (select.T) {
    norms_list <- list()
    for (i in 1:N) {
      message(paste0("Find similar cells for batch ", i))
      
      norms_list[[i]] <- lapply(object@norm.data,
                                function(x) cosineNorm(x[object@var_list[[i]], ]))
      pca_list <- run.pca(norms_list[[i]], i, dims)
      
      modis[[i]] <- modify_simi(d_i = pca_list[[1]], d_remained = pca_list[[2]], cluster_labels = cluster_iden[[i]],
                                        dims = dims, k = k)
      rm(pca_list)
    }

    T_th <- select_T(cand_Ts, cluster_iden, modis, norms_list, dims)
    message(paste0("Optimal T_th = ", as.character(T_th)))
    message("Filtering")
    modis <- lapply(1:N, function(i) filter.cpca(modis[[i]], norms_list[[i]], i, dims, T_th))
    message("Done")
    #rm()
  } else {
    for (i in 1:N) {
      message(paste0("Find similar cells for batch ", i))
      
      norms <- lapply(object@norm.data,
                      function(x) cosineNorm(x[object@var_list[[i]], ]))
      pca_list <- run.pca(norms, i, dims)
      
      modi <- modify_simi(d_i = pca_list[[1]], d_remained = pca_list[[2]], cluster_labels = cluster_iden[[i]],
                          dims = dims, k = k) %>% filter.cpca(., norms, i, dims, T_th)
      
      modis[[i]] <- modi
      message("Done!")
      rm(norms, pca_list, modi)
    }
  }

  S <- lapply(1:N,
              function(i){
                ss <- Matrix(0, nrow = NN, ncol = NN, sparse = TRUE)
                rownames(ss) = colnames(ss) = Reduce('c', cell.names)
                remained_names <- Reduce('c', cell.names[-i])[modis[[i]]$indexs[[2]]]
                ss[cell.names[[i]][modis[[i]]$indexs[[1]]], remained_names] <- modis[[i]]$modi_S
                ss[remained_names, cell.names[[i]][modis[[i]]$indexs[[1]]]] <- Matrix::t(modis[[i]]$modi_S)
                ss
              })

  S <- Reduce('+', S) %>% as(., 'dgCMatrix')

  groups <- Reduce('c', lapply(1:N, function(i) rep(i, length(cell.names[[i]])))) - 1
  sparseS <- norm_S_batches(S@x, S@p, S@i, NN, groups, N)
  S <- sparseMatrix(i = as.vector(sparseS$i), j = as.vector(sparseS$j), x = as.vector(sparseS$x), dims = c(NN, NN))

  t2 = Sys.time()
  print(t2-t1)

  object@S <- S
  return(object)
}

select_T <- function(cand_Ts, cluster_iden, modis, norms_list, dims) {
  # select T_th
  message(paste0("Select optimal T_th"))
  # step 1
  cand_Ts <- sort(cand_Ts)
  cand_list <- list()
  cell_identities <- list()
  N <- length(norms_list)
  for (i in 1:N) {
    retain_cells <- list()
    cell_identities[[i]] <- list()
    SimiNum <- c()
    remaining_identities <- Reduce('c', cluster_iden[-i])
    modi_T <- filter.cpca(modis[[i]], norms_list[[i]], i, dims, cand_Ts)
    retain_cells <- lapply(1:length(cand_Ts), function(s) which(colSums(modi_T$modi_S[[s]]) > 0))
    SimiNum <- sapply(1:length(cand_Ts), function(s) length(retain_cells[[s]]))
    cell_identities[[i]] <- lapply(1:length(cand_Ts), function(s) unique(remaining_identities[retain_cells[[s]]]))
    D_T <- as.vector(c(- diff(SimiNum) / (SimiNum[-length(cand_Ts)] + 1e-4), 0.51))
    cand_list[[i]] <- cand_Ts[1:which(D_T >= 0.5)[1]]
  }
  New_cands <- sort(Reduce(intersect, cand_list))
  
  # step 2
  cell_identities <- lapply(1:length(New_cands),
                            function(j) Reduce(union, lapply(cell_identities, function(x) x[j])))
  overlaps <- lapply(1:length(New_cands),
                     function(j) sapply(1:N, function(l) length(intersect(cell_identities[[j]][[l]], cell_identities[[length(New_cands)]][[l]])) /
                                          length(cell_identities[[length(New_cands)]][[l]])))
  T_n <- sapply(1:N, function(n) {
    x = sapply(overlaps, function(o) o[n])
    which(x == 1)[1]
  })
  optimal_T <- median(New_cands[T_n])
  return(optimal_T)
}

#' Compute the cell-cell similarities
#'
#' @param d_i The dataset i in PCA space.
#' @param d_remained Remaining datasets of dataset i in PCA space.
#' @param cluster_labels The pre-clustered labels of dataset i.
#' @param dims The dimension of PCA. (default 30)
#' @param k The number of computed top similarities for each cell. (default 5)
#' @param th The similarity threshold for filtering cells in the PCA space. (default 0.6)
#'
#' @return List of identified cell-cell similarities.
#'
#' @import Matrix
#' @importFrom Biobase rowMax
modify_simi <- function(d_i, 
                        d_remained, 
                        cluster_labels, 
                        dims = 30, 
                        k = 5, 
                        th = 0.6) {
  # cosine similarities
  cos_d_i <- cosineNorm(t(d_i))
  cos_d_remained <- cosineNorm(t(d_remained))
  cos_sim <- eigenMapMatcrossprod(cos_d_i, cos_d_remained)
  rm(cos_d_i, cos_d_remained)

  max_index <- max.col(cos_sim)
  max_sim <- rowMax(cos_sim)
  sim_index <- ifelse(max_sim < th, 0, 1)

  dis.mat <- sqrt(pmax(2*(1 - cos_sim), 0))
  rm(cos_sim)

  sigma <- sapply(1:nrow(d_i), function(j) dis.mat[j, max_index[j]])
  for (lb in unique(cluster_labels)) {
    sigma[which(cluster_labels == lb & sim_index >= th)] <- quantile(sigma[which(cluster_labels == lb & sim_index >= th)], 0.5)
  }

  modi_S <- sim_index * exp(-abs(dis.mat - sigma))
  modi_S <- t(apply(modi_S, 1, function(x) {
    x[-order(x, decreasing = T)[1:k]] <- 0
    x
  }))
  modi_S[which(modi_S < th)] <- 0
  modi_S <- as(modi_S, 'dgCMatrix')
  c_i <- which(Matrix::rowSums(modi_S) > 0)
  c_remained <- which(Matrix::colSums(modi_S) > 0)
  return(list("modi_S" = modi_S[c_i, c_remained],
              "indexs" = list("i" = c_i,
                              "remains" = c_remained)
  ))
}

#' Compute specific PCA dimensionality reduction for each dataset i
#'
#' @param norms The list of normalized matrices.
#' @param i The dataset i.
#' @param dims The dimension of PCA.
#'
#' @return List of PCA results.
#'
#' @importFrom matrixStats rowSds
#' @importFrom irlba irlba
run.pca <- function(norms, 
                    i, 
                    dims) {
  mean_i <- Matrix::rowMeans(norms[[i]])
  sds_i <- matrixStats::rowSds(as.matrix(norms[[i]]))
  scales <- lapply(norms, function(x) as.matrix(scale_data(x, row.means = mean_i, row.sds = sds_i)))
  project_i <- irlba::irlba(t(scales[[i]]), nv = dims) %>% "[["("v")
  pca_list <- lapply(scales, function(x) t(eigenMapMatcrossprod(project_i, x)))
  return(pca_list = list('d_i' = pca_list[[i]],
                         'd_remained' = do.call(rbind, pca_list[-i])))
}

#' Filter cell-cell similarities for each dataset i in cPCA space
#'
#' @param modi The list of cell-cell similarities
#' @param norms The list of normalized matrices.
#' @param i The dataset i.
#' @param dims The dimension of cPCA.
#' @param T_th The (candidate) similarity threshold(s) for filtering cells in the cPCA space.
#'
#' @return List of filtered cell-cell similarities.
#'
#' @importFrom rfunctions crossprodcpp
filter.cpca <- function(modi, 
                        norms, 
                        i, 
                        dims, 
                        T_th) {
  if (length(modi$indexs$i) == 0) return(modi)
  mean_i <- Matrix::rowMeans(norms[[i]][, modi$indexs$i])
  scale_i <- as.matrix(scale_data(norms[[i]][, modi$indexs$i], row.means = mean_i, do.scale = F))
  scale_remained <- do.call(cbind, norms[-i])
  scale_remained <- as.matrix(scale_data(scale_remained[, modi$indexs$remains], row.means = mean_i, do.scale = F))
  cov_i <- rfunctions::crossprodcpp(t(scale_i)) / ncol(scale_i)
  cov_remained <- rfunctions::crossprodcpp(t(scale_remained)) / ncol(scale_remained)
  cpca_project <- (cov_i - cov_remained) %>% RSpectra::eigs_sym(k = dims, which = "LA") %>% "[["("vectors") %>% as.matrix()

  e_i <- cosineNorm(eigenMapMatcrossprod(cpca_project, scale_i))
  e_remained <- cosineNorm(eigenMapMatcrossprod(cpca_project, scale_remained))
  if (length(T_th) > 1) {
    us_list <- lapply(T_th, 
                      function(t_th) filter_T(modi$modi_S@x, modi$modi_S@p, modi$modi_S@i, ncol(modi$modi_S), t_th, e_i, e_remained))
    modi$modi_S <- lapply(us_list,
                          function(us) sparseMatrix(i = as.vector(us$i), j = as.vector(us$j), x = as.vector(us$x), dims = c(nrow(modi$modi_S), ncol(modi$modi_S))))
  } else {
    us <- filter_T(modi$modi_S@x, modi$modi_S@p, modi$modi_S@i, ncol(modi$modi_S), T_th, e_i, e_remained)
    modi$modi_S <- sparseMatrix(i = as.vector(us$i), j = as.vector(us$j), x = as.vector(us$x), dims = c(nrow(modi$modi_S), ncol(modi$modi_S)))
  }
  return(modi)
}

#' Perform reference-based mapping.
#'
#' This function performs reference-based mapping under three mapping models. scInt using "global" model by default to learned contrastive biological variation uniformly.
#'
#' @param ref.data The reference \code{scInt} object.
#' @param query.data The query \code{scInt} object.
#' @param npcs The dimension of projection matrix. (default 40)
#' @param model The mapping model.
#' @param var.model The model to combine features of reference and query.
#' @param resolution The resolution of Leiden algorithm in the "global" model. (default 0.6)
#' @param k The number of computed top similarities for each cell. (default 5)
#' @param T_th The similarity threshold for filtering cells in the cPCA space. (default 0.7)
#' @param lambda The regularization parameter of technical variation. (default 5)
#'
#' @return Reference-based integration result.
#'
#' @import methods
#' @import Matrix
#'
#' @export
run.Mapping <- function(ref.data, 
                        query.data, 
                        npcs = 40, 
                        model = c("global", "projection", "cpca"), 
                        var.model = c("intersect", "union"),
                        resolution = 0.6, 
                        k = 5, 
                        T_th = 0.6, 
                        lambda = 5) {
  model <- match.arg(model)
  var.model <- match.arg(var.model)
  t1 = Sys.time()

  if (var.model == "intersect") {
    if (model == "atlas") stop("var.model must be 'union' under 'atlas' integration model.")
    vargs <- ref.data@vargenes
    # ref.Z
    ref.Z <- do.call(cbind, lapply(ref.data@norm.data, function(x) x[vargs, ]))
    # query.Z
    query.Z <- as(matrix(0, nrow = length(vargs), ncol = ncol(do.call(cbind, query.data@norm.data))), "sparseMatrix")
    colnames(query.Z) <- Reduce("c", lapply(query.data@norm.data, colnames))
    rownames(query.Z) <- vargs
    common_vargs <- intersect(ref.data@vargenes, rownames(query.data@norm.data[[1]]))
    if (length(common_vargs) < 500) warning("It is recommended to use 'atlas' model because of too few common vargenes.")
    query.Z[common_vargs, ] <- Reduce(cbind2, lapply(query.data@norm.data,
                                                     function(x) x[common_vargs, ]))
    if (model == "projection") {
      ref.projection <- ref.data@projection
    }

    if (model == "cpca") {
      ref.var <- ref.data@variation
      ref.mean <- Matrix::rowMeans(ref.Z[common_vargs, ])
      query.matrix <- as.matrix(query.Z)
      query.matrix[common_vargs, ] <- scale_data(query.matrix[common_vargs, ], do.scale = F, do.center = T, row.means = ref.mean)#F, row.means = ref.mean
      query.var <- as.matrix(Matrix::tcrossprod(query.matrix)) / ncol(query.matrix)
    }

    if (model == "global") {
      ref.var <- ref.data@variation
      ref.projection <- ref.data@projection

      if (length(query.data@norm.data) >= 2 && ncol(query.data@variation) != 0) {
        inter_vargs <- intersect(rownames(query.data@variation), vargs)
        if (length(inter_vargs) < 500) warning("It is recommended to use 'atlas' model because of too few common vargenes..")
        query.var <- query.data@variation[inter_vargs, inter_vargs]
        query.projection <- query.var %>%
          RSpectra::eigs_sym(k = npcs, which = "LA") %>%
          "[["("vectors") %>%
          as.matrix() %>%
          t()
      } else {
        if (length(query.data@norm.data) >= 2 && ncol(query.data@variation) == 0) warning("It is recommended to integrate the query data first.")
        query.var <- query.Z %>%
          scale_data(do.scale = F, do.center = T) %>%
          Matrix::tcrossprod() %>%
          as.matrix() %>%
          "/"(ncol(query.Z))
        query.projection <- query.var %>%
          RSpectra::eigs_sym(k = npcs, which = "LA") %>%
          "[["("vectors") %>%
          as.matrix() %>%
          t()
      }
    }
  } else if (var.model == "union") {
    if (length(query.data@vargenes) == 0) stop("query.data must have vargenes.")

    vargs <- union(ref.data@vargenes, query.data@vargenes)
    # ref.Z
    ref.Z <- as(matrix(0, nrow = length(vargs), ncol = ncol(do.call(cbind, ref.data@norm.data))), "sparseMatrix")
    colnames(ref.Z) <- Reduce("c", lapply(ref.data@norm.data, colnames))
    rownames(ref.Z) <- vargs
    ref.vargs <- intersect(rownames(ref.data@norm.data[[1]]), vargs)
    ref.Z[ref.vargs, ] <- Reduce(cbind2, lapply(ref.data@norm.data,
                                                function(x) x[ref.vargs, ]))
    ref.var <- matrix(0, nrow = length(vargs), ncol = length(vargs))
    rownames(ref.var) = colnames(ref.var) = vargs
    ref.var[ref.data@vargenes, ref.data@vargenes] <- ref.data@variation

    # query.Z
    query.Z <- as(matrix(0, nrow = length(vargs), ncol = ncol(do.call(cbind, query.data@norm.data))), "sparseMatrix")
    colnames(query.Z) <- Reduce("c", lapply(query.data@norm.data, colnames))
    rownames(query.Z) <- vargs
    query.vargs <- intersect(rownames(query.data@norm.data[[1]]), vargs)
    query.Z[query.vargs, ] <- Reduce(cbind2, lapply(query.data@norm.data,
                                                    function(x) x[query.vargs, ]))

    if (model == "projection") {
      ref.projection <- ref.var %>%
        RSpectra::eigs_sym(k = npcs, which = "LA") %>%
        "[["("vectors") %>%
        as.matrix() %>%
        t()
    }
    if (model == "cpca") {
      query.var <- matrix(0, nrow = length(vargs), ncol = length(vargs))
      rownames(query.var) = colnames(query.var) = vargs
      query.var[query.data@vargenes, query.data@vargenes] <- query.data@variation
    }
    if (model == "atlas") {
      ref.projection <- ref.var %>%
        RSpectra::eigs_sym(k = npcs, which = "LA") %>%
        "[["("vectors") %>%
        as.matrix() %>%
        t()
      query.var <- matrix(0, nrow = length(vargs), ncol = length(vargs))
      rownames(query.var) = colnames(query.var) = vargs
      query.var[query.data@vargenes, query.data@vargenes] <- query.data@variation
      query.projection <- query.var %>%
        RSpectra::eigs_sym(k = npcs, which = "LA") %>%
        "[["("vectors") %>%
        as.matrix() %>%
        t()
    }
  }

  if (model == "projection") {
    # reference projection
    Z <- cbind(ref.Z, query.Z)
    L <- eigenMapMatMult(ref.projection, as.matrix(Z))
    colnames(L) <- colnames(Z)
    t2 = Sys.time()
    print(t2-t1)
    return(L)
  }

  if (model == "cpca") {
    if (var.model == "union") {
      cpca_project <- (ref.var + query.var) %>% RSpectra::eigs_sym(k = npcs, which = "LA") %>% "[["("vectors") %>% as.matrix() %>% t()
    } else {
      cpca_project <- (ref.var - query.var) %>% RSpectra::eigs_sym(k = npcs, which = "LA") %>% "[["("vectors") %>% as.matrix() %>% t()
    }

    Z <- cbind(ref.Z, query.Z)
    L <- eigenMapMatMult(cpca_project, as.matrix(Z))
    colnames(L) <- colnames(Z)
    t2 = Sys.time()
    print(t2-t1)
    return(L)
  }

  if (model == "global" || model == "atlas") {
    labels_ref <- ref.data@integrated.data %>%
      t() %>%
      FindNeighbors(verbose = FALSE) %>%
      "[["("snn") %>%
      FindClusters(resolution = resolution, verbose = FALSE) %>%
      '[['(1) %>%
      as.character() %>%
      paste0("ref_", .)
    if (ncol(query.data@integrated.data) != 0) {
      labels_query <- query.data@integrated.data %>%
        t() %>%
        FindNeighbors(verbose = FALSE) %>%
        "[["("snn") %>%
        FindClusters(resolution = resolution, verbose = FALSE) %>%
        '[['(1) %>%
        as.character() %>%
        paste0("query_", .)
    } else {
      labels_query <- query.Z %>%
        scale_data() %>%
        as.matrix() %>%
        RunPCA(npcs = npcs, assay = "RNA", verbose = FALSE) %>%
        "@"("cell.embeddings") %>%
        FindNeighbors(verbose = FALSE) %>%
        "[["("snn") %>%
        FindClusters(resolution = resolution, verbose = FALSE) %>%
        '[['(1) %>%
        as.character() %>%
        paste0("query_", .)
    }

    message("Compute similarities under reference projection.")
    d_ref <-as.matrix(ref.projection %*% ref.Z) %>% t()
    d_query <- as.matrix(ref.projection %*% query.Z) %>% t()
    modi_ref <- modify_simi(d_i = d_ref, d_remained = d_query, cluster_labels = labels_ref,
                            dims = npcs, k = k)
    # run cpca
    query.m <- query.Z[, modi_ref$indexs$remains]
    query.v <- query.m %>%
      scale_data(do.scale = F, do.center = T) %>%
      Matrix::tcrossprod() %>%
      as.matrix() %>%
      "/"(length(modi_ref$indexs$remains))
    cpca_project <- (ref.var - query.v) %>% RSpectra::eigs_sym(k = npcs, which = "LA") %>% "[["("vectors") %>% as.matrix() %>% t()
    e_ref <-  cosineNorm(eigenMapMatMult(cpca_project, as.matrix(ref.Z)))
    e_query <- cosineNorm(eigenMapMatMult(cpca_project, as.matrix(query.m)))
    us <- filter_T(modi_ref$modi_S@x, modi_ref$modi_S@p, modi_ref$modi_S@i, ncol(modi_ref$modi_S), T_th, e_ref, e_query)
    modi_ref$modi_S <- sparseMatrix(i = as.vector(us$i), j = as.vector(us$j), x = as.vector(us$x), dims = c(nrow(modi_ref$modi_S), ncol(modi_ref$modi_S)))

    # merge the ref and query similarities
    mS <- Matrix(0, nrow = ncol(ref.Z), ncol = ncol(query.Z), sparse = TRUE)
    mS[modi_ref$indexs$i, modi_ref$indexs$remains] <- modi_ref$modi_S

    groups.ref <- rep(0, ncol(ref.Z))
    sparse_mS <- norm_S_batches(mS@x, mS@p, mS@i, ncol(query.Z), groups.ref, 1)#length(ref.Ns)
    mS.ref <- sparseMatrix(i = as.vector(sparse_mS$i), j = as.vector(sparse_mS$j), x = as.vector(sparse_mS$x), dims = c(ncol(ref.Z), ncol(query.Z)))

    # merge the similarity matrix
    S.ref <- ref.data@S
    if (ncol(query.data@S) == 0) {
      S.query <- Matrix(0, nrow = ncol(query.Z), ncol = ncol(query.Z), sparse = TRUE)
    } else {
      S.query <- query.data@S
    }
    S <- Matrix::bdiag(S.ref, S.query)
    uS <- upper_tri(S@x, S@p, S@i, ncol(S))
    S <- sparseMatrix(i = as.vector(uS$i), j = as.vector(uS$j), x = as.vector(uS$x),
                      dims = c(nrow(S), ncol(S)), repr = "C")
    rownames(S) = colnames(S) = c(colnames(ref.Z), colnames(query.Z))
    S[colnames(ref.Z), colnames(query.Z)] <- mS.ref

    message("Mapping...")
    N <- length(ref.data@norm.data) + length(query.data@norm.data)
    cell.names <- lapply(c(ref.data@norm.data, query.data@norm.data), colnames)
    nums <- sapply(cell.names, length)
    # ref and query tcross_S
    tcross_S <- lapply(1:N,
                       function(ns) {
                         x <- Matrix(0, nrow = nrow(S), ncol = nrow(S), sparse = T)
                         colnames(x) = rownames(x) <- Reduce('c', cell.names)
                         x[cell.names[[ns]], ] <- S[cell.names[[ns]], ]
                         num <- lapply(cell.names, function(j) Matrix::colSums(x[, j]))
                         x <- Matrix::Diagonal(x = Reduce('c', num)) - x
                         if (sum(sapply(num, sum)) > 0) {
                           x <- x/sqrt(sum(sapply(num, sum)))
                         }
                         Matrix::tcrossprod(x)
                       })
    tcross_S <- Reduce('+', tcross_S)

    Z <- scale_data(cbind(ref.Z, query.Z), do.scale = F)
    K <- Matrix::Diagonal(x = Reduce('c', lapply(nums, function(num) rep(1/num, num)))) - lambda * tcross_S
    d <- Z %*% Matrix::tcrossprod(K, Z)

    result_cpca <- RSpectra::eigs_sym(as.matrix(d), k = npcs, which = "LA") %>% '[['('vectors') %>% as.matrix()
    L = as.matrix(Matrix::crossprod(result_cpca, cbind(ref.Z, query.Z)))
    colnames(L) <- c(colnames(ref.Z), colnames(query.Z))
    t2 = Sys.time()
    print(t2-t1)
    return(L)
  }
}


cosineNorm <- function(x) {
  l2 <- sqrt(Matrix::colSums(x^2))
  l2 <- pmax(1e-8, l2)
  mat <- scale(x, center = F, scale = l2)
  mat
}
