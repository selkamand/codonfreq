
#' Codon Frequency
#'
#' @return dataframe describing pre-calculated rates at which codons occur in the coding genome
#' @export
#'
#' @examples
#' codon_freq()
codon_freq <- function(){
  path_to_codonfreq_file = system.file("extdata/codonusage.genscript.tsv", package = "codonfreq")
  read.csv(header = TRUE, file = path_to_codonfreq_file, sep = "\t", check.names = FALSE) %>%
    dplyr::mutate(`Frequency/Thousand` = Number / sum(Number) * 1000) %>%
    dplyr::tibble()
}

codon_freq_stats <- function(){
  codon_freq() %>%
    ggplot2::ggplot(ggplot2::aes(`Frequency/Thousand`)) +
    ggplot2::geom_density(alpha = 0.5, fill = "steelblue")

  #codon_freq() %>%
    #ggplot2::ggplot(ggplot2::aes(x="bla", y = `Frequency/Thousand`)) +
    #ggplot2::geom_jitter(alpha = 0.5, fill = "steelblue")
}

#' Codon Frequency
#'
#' Take a fasta file and return a dataframe describing how common its codons are in the human genome.
#'
#' @param fasta_path path to fasta file (string)
#' @param write_bedgraph write bedgraph files for each fasta sequence and reading frame? (bool)
#'
#' @return table describing frequency of each codon in the fasta, as calculated using the entire coding genome (dataframe)
#' @export
#'
check_fasta_codon_bias <- function(fasta_path, write_bedgraph = TRUE){
  assertthat::assert_that(file.exists(fasta_path))

  fasta <- seqinr::read.fasta(file = fasta_path, seqtype = "DNA", forceDNAtolower = TRUE)
  fasta_codon_freq_df <- purrr::map(
    fasta, ~
      purrr::map(c(-3, -2, -1, 1, 2, 3), get_codon_frequency, .x) %>%
      magrittr::set_names(c(-3, -2, -1, 1, 2, 3)) %>%
      purrr::map_df(~dplyr::tibble(.x), .id="frame")
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(id = paste0(seqname, "_codon_bias_frame_", frame))


  if(write_bedgraph){
    message("Writing Bedgraph Files in current directory")
    fasta_codon_freq_df %>%
      dplyr::group_by(id) %>%
      dplyr::group_walk(
        ~ .x %>%
          dplyr::select(-frame) %>%
          write.table(x = ., file = paste0(.y, ".bedgraph"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE))
  }
  return(fasta_codon_freq_df)
}
#
# fasta_path="/Users/selkamand/garage/cd74_pdgfrb_fusion/test.fasta"
#
# fasta <- seqinr::read.fasta(file = fasta_path, seqtype = "DNA", forceDNAtolower = TRUE)
#
#
# #fasta <- append(x = fasta, values = purrr::map(fasta, ~ rev(seqinr::comp(.x)) ), after = Inf)
# #fasta
#
# #for (seq in fasta){
# #purrr::map(fasta, ~purrr::map_dfc(0:2, get_codon_frequency, .x))
#  #   print()
# #}
# #
# # fasta$MyTest %>%
# #   seqinr::splitseq(., frame = 2, word = 3)
# #
# # purrr::map(fasta, ~purrr::map_dfc(1:3, get_codon_frequency, .x))
#
#
# fasta_codon_freq_df <- purrr::map(fasta, ~
#              purrr::map(c(-3, -2, -1, 1, 2, 3), get_codon_frequency, .x) %>%
#              magrittr::set_names(c(-3, -2, -1, 1, 2, 3)) %>%
#              purrr::map_df(~dplyr::tibble(.x), .id="frame")
#            ) %>%
#   dplyr::bind_rows() %>%
#   dplyr::mutate(id = paste0(seqname, "_codon_bias_frame_", frame))
#
# fasta_codon_freq_df
#
# fasta_codon_freq_df %>%
#   dplyr::count(seqname)
#
# fasta_codon_freq_df %>%
#   dplyr::group_by(id) %>%
#   dplyr::group_walk(~ .x %>%
#                       dplyr::select(-frame) %>%
#                       write.table(x = ., file = paste0(.y, ".tsv"), sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE))
#
#
# fasta$MyTest2 %>% get_codon_frequency(frame = 1, seq = .)
#
# codon_freq() %>%
#   dplyr::filter(Triplet=="CCA")
#
# codon_freq() %>%
#   dplyr::filter(abs(`Frequency/Thousand`- 22.3) < 0.09)
#


get_codon_frequency <- function(frame = 1, seq){
  seqname = attr(seq, which = "name")
  #message("Sequence Name:", seqname)

  assertthat::assert_that(frame %in% c(-3, -2, -1, 1, 2, 3))
  reverse = FALSE
  if (frame < 0){
    reverse = TRUE
    seq = rev(seqinr::comp(seq))
  }

  seqinr_frame = abs(frame) -1
  codon_freq_df <- codon_freq()

  seqinr::splitseq(seq, frame = seqinr_frame, word = 3) %>%
    match(., tolower(codon_freq_df$Triplet)) %>%
    codon_freq_df$`Frequency/Thousand`[.] %>%
    conditional_reverse(reverse = reverse) %>%
    dplyr::tibble(
      seqname = seqname,
      start_pos = seq(from = abs(frame) - 1, by = 3, length.out=length(.)),
      end_pos = start_pos + 3,
      freq=.
    )
}

conditional_reverse <- function(vector, reverse = FALSE){
  if (reverse){
    return(rev(vector))
  }
  else
    return(vector)
}

