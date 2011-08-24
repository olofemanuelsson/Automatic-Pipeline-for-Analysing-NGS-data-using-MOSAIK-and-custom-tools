#Prints nonvarsites for every sample
#USAGE: Rscript get.sites.R varcalls.infile outdir


argv = commandArgs(trailingOnly=TRUE)
varcalls.file = argv[1]
outdir = argv[2]


if(0){ #no batch calling
  datadir = '../res/varcall'
  varcalls.file = 'annovar_master_all_subject_variants.txt'
  varcalls.file = file.path(datadir, varcalls.file)
}

main <- function(varcalls.file, outdir){
  
  varcalls = read.varcalls(varcalls.file)
  sample.nonvarsites = get.nonvarsites(varcalls)
  print.nonvarsites(sample.nonvarsites, outdir)
  
}

print.nonvarsites <- function(sample.nonvarsites, outdir){
  #print nonvarsites for every sample. mpileup -l format: chr startpos-endpos. NB: no "chr" prefix on chromosomes to comply with the reference file chrom-naming.
  samples = names(sample.nonvarsites)
  for(jsample in samples){
    sample.nonvarsites.file = file.path(outdir, paste(jsample, 'nonvarsites.txt', sep='.'))
    jsample.nonvarsites = sample.nonvarsites[[jsample]]
    start.end = jsample.nonvarsites[, c('start', 'stop')]
    start.end = apply(start.end, 1, paste, collapse='-')
    chr = gsub('^chr', '', jsample.nonvarsites[, 'chr'])
    chr.start.end = cbind(chr, start.end)
    write.table(chr.start.end, sep=" ", quote = FALSE, file = sample.nonvarsites.file, row.names = FALSE, col.names = FALSE)
  }
}

read.varcalls <- function(varcalls.file){

  varcalls = read.table(varcalls.file, sep='\t', stringsAsFactors = FALSE)
  #rm last col (empty)
  varcalls = varcalls[, 1:(ncol(varcalls) -1)]

  first.colnames = c('chr', 'start', 'stop', 'refvar', 'observed var')
  last.colnames = c('gene relation', 'gene annot', 'var type', 'aa change', 'ConservedReg', 'SegDup', '1000g', 'Dbsnp_1', 'Dbsnp_2', 'cg46', 'Avsift', 'PP2')

  colnames(varcalls)[1:length(first.colnames)] = first.colnames
  colnames(varcalls)[(ncol(varcalls) - length(last.colnames) + 1):ncol(varcalls)] = last.colnames

  #set colnames for sample columns
  sample.cols = (length(first.colnames) + 1):(ncol(varcalls) - length(last.colnames))
  n.samples = length(sample.cols)

  for(jsample.col in sample.cols){
    var.ind = which(varcalls[, jsample.col] != '-')        
    sample.id = sub('(.*)(Hom|Het)', '\\1', varcalls[var.ind[1], jsample.col])
    sample.id = strsplit(sample.id, ';')[[1]][1]
    colnames(varcalls)[jsample.col] = sample.id
  }

  return(varcalls)
}

get.nonvarsites <- function(varcalls){

  #get sample columns
  sample.cols.start = which(colnames(varcalls) == 'observed var') + 1
  sample.cols.end = which(colnames(varcalls) == 'gene relation') - 1
  samples = colnames(varcalls)[sample.cols.start:sample.cols.end]

  #get sample-specific sites with no variation
  sample.nonvarsites = list()
  length(sample.nonvarsites) = length(samples)
  names(sample.nonvarsites) = samples
  for(jsample in samples){
    sample.nonvarsites[[jsample]] = varcalls[which(varcalls[, jsample] == '-'), c('chr', 'start', 'stop')]   
  }

  return(sample.nonvarsites)  
}


#Execute
main(varcalls.file, outdir)
