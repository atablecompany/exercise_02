if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

library(Biostrings)
library(seqinr)

seq_b = readDNAStringSet('fishes.fna.gz')
seq_s = read.fasta('fishes.fna.gz')

length(seq_b)  # pocet sekvenci
width(seq_b[1])  # delka sekvence
names(seq_b)  # nazvy sekvenci z fasta headeru
seq1 = seq_b[1]
seq1_sequence = seq_b[[1]]
seq1_string = toString(seq_b[1])
seq2_string = toString(seq_b[2])

pairwiseAlignment(seq1_string, seq2_string, substitutionMatrix=BLOSUM62, gapOpening=-1, gapExtension=1)


## regular expressions
names_list <- c("anna", "jana", "kamil", "norbert", "pavel", "petr", "stanislav", "zuzana")

grep('jana', names_list, perl=TRUE)
grep('n+', names_list, perl=TRUE)  # contains the letter 'n' at least once
grep('n{2}', names_list, perl=TRUE)  # contains 'nn'
grep('^n', names_list, perl=TRUE)  # starts with 'n'
grep('anna|jana', names_list, perl=TRUE)  # 'anna' or 'jana'
grep('^z.*a$', names_list, perl=TRUE)  # starts with 'z' and ends with 'a'


## task 5
seq = readDNAStringSet('fishes.fna.gz')

y = grep('^ACGAGTGCGT.*ACGCACTCGT$', seq, perl=TRUE)
y
print(seq[25001])

demultiplexer = function(fasta_file, fmids, rmids, labels){
  fasta_file = 'fishes.fna.gz'
  seq = readDNAStringSet(fasta_file)
  res = list()
  results = data.frame(labels,'')
  
  for (i in 1:length(fmids)){
    expr = paste0('^',fmids[i],'.*',reverseComplement(DNAString(rmids[i])),'$')
    y = grep(expr, seq, perl=TRUE)
    
    z = list(y)
    res = append(res, z)
  }
  
  report = c()
  for (i in 1:length(res)){
    report = append(report, paste0(labels[i], ': ', length(res[[i]])))
    report_file = file('report.txt')
    writeLines(report, report_file)
    close(report_file)
  }
  
  for (i in 1:length(labels)){
    seqs_to_write = c()
    for (j in 1:length(res[[i]])){
      seq_to_write = seq[res[[i]][j]]
      seq_to_write = substr(seq_to_write, 1+nchar(fmids[i]), nchar(seq_to_write)-nchar(rmids[i]))
      seqs_to_write = append(seqs_to_write, seq_to_write)  
    }
    seqs_to_write = DNAStringSet(seqs_to_write)
    writeXStringSet(seqs_to_write, paste0(labels[i], '.fasta'))
  }
}

mids = read.csv('fishes_MIDs.csv',sep=';')
fmids = mids$FBarcodeSequence
rmids = mids$RBarcodeSequence
labels = mids$Description

demultiplexer('fishes.fna.gz', fmids, rmids, labels)