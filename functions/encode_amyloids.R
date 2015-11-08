#recover encoding from the string
vector2list <- function(x){
  require(dplyr)
  lapply(x, function(i) {
    pasted_group <- strsplit(i, "_", fixed = TRUE)
    lapply(pasted_group, function(j) strsplit(j, ""))
  }) %>% lapply(first)
}

create_gl <- function(seq)
  lapply(1L:nrow(seq), function(i) {
    res <- do.call(rbind, strsplit(decode_ngrams(seq2ngrams(seq[i, ][!is.na(seq[i, ])], 6, a()[-1])), ""))
    cbind(res, id = paste0("P", rep(i, nrow(res))))
  })


pos_gl <- create_gl(pos_data)
neg_gl <- create_gl(neg_data)

extract_ngrams <- function(gl) {
  require(biogram)
  lapply(aa_groups, function(single_group) 
    lapply(gl, function(single_protein) {
      bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                               ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                               seq = degenerate(single_protein[, -7], single_group),
                                               u = as.character(1L:length(single_group))))
      
      bitrigrams <- bitrigrams > 0
      storage.mode(bitrigrams) <- "integer"
      
      bitrigrams
    })
  )
}

coded_neg <- lapply(aa_groups, function(single_group) {
  lapply(neg_gl, function(single_protein) {
    bitrigrams <- as.matrix(count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)), 
                                             ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)),
                                             seq = degenerate(single_protein[, -7], single_group),
                                             u = as.character(1L:length(single_group))))
    
    bitrigrams <- bitrigrams > 0
    storage.mode(bitrigrams) <- "integer"
    
    bitrigrams
  })
})
