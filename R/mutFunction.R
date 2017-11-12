#' Mutant Calling Function for Sanger Sequencing
#'
#' This function will allow you to automatically call heterozygous mutations from Sanger sequencing
#' @param test.seq Input the sequence you wish to test (as a string pointing to a file in your working directory).
#' @param ref Input the reference sequencing (you did sequence a control, right?)
#' @param length Input the approximate length of your usable sequence as a positive integer.
#' @param trim5 You can manually adjust the amount of trimming done from the 5' end - defaults to 25.
#' @param trim3 You can manually adjust the amount of trimming done from the 3' end - defaults to 25.
#' @export
#' @import sangerseqR
#' @examples
#' mutantCaller()
#' P1               101 ACGACCAGCTCCAGGAGTACGTCCAGCAGACCGTCG--GCCCAGAAACAC    148
#'                      ||||||||||||||||||||||||||||||||||||  ||||||||||||
#' S1               101 ACGACCAGCTCCAGGAGTACGTCCAGCAGACCGTCGCAGCCCAGAAACAC    150
#' This is a two base pair insertion because the secondary sequence is longer than the primary.


mutantCaller <- function(test.seq, ref, length, trim5 = 25, trim3 = 25)
{
  reference <- subseq(primarySeq(readsangerseq(ref), string = T), start = 25, width = length)
  var.seq <- readsangerseq(test.seq)
  var.calls <- makeBaseCalls(var.seq, ratio = 0.2)
  var.alleles <- setAllelePhase(var.calls, reference, trim5 = trim5, trim3 = trim3)
  var.align <- pairwiseAlignment(primarySeq(var.alleles)[1:length], secondarySeq(var.alleles)[1:length], type = 'global-local')
  return(writePairwiseAlignments(var.align))
}
