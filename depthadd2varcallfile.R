#Reads varcall "master" file and depth-info for non-var sites of every sample and writes new varcall "master" file with added non-var depth info.
#
#USAGE: Rscript depthadd2varcallfile.R varcalls.infile varcalls.outfile
#NB: Assumes that get.siteinfo.pl has been called such that the depth for every sample are found in files named "*.nonvarsites.depth", residing in the same dir as varcalls.infile.

argv = commandArgs(trailingOnly=TRUE)
varcalls.infile = argv[1]
varcalls.outfile = argv[2]

main <- function(varcalls.infile, varcalls.outfile){

  #Put non-var depths from every sample into one single file (skip the base quality column since the read.table fcn of R for some reason then can find tab chars in some lines) 
  datadir = dirname(varcalls.infile)
  nonvars.file = file.path(datadir, 'nonvarsites.depth')
  cmd = sprintf('find %s %s %s', datadir, "-name '*.nonvarsites.info' | xargs -I% awk -F'\t' -v OFS='\t' '{print FILENAME, $1, $2, $4;}' % >", nonvars.file)
  system(cmd)

  #Read varcalls
  varcalls = read.varcalls(varcalls.infile)
  n.vars = nrow(varcalls)
  print(sprintf('n vars: %i', n.vars))

  #add non-var depths to varcalls table
  varcalls = add.nonvar.depth(varcalls, nonvars.file)

  #Write new table with added depths
  write.cols = setdiff(colnames(varcalls), 'chr.start')
  write.table(varcalls[, write.cols], file = varcalls.outfile, sep='\t', row.names=F, quote=F)

}

add.nonvar.depth <- function(varcalls, nonvars.file){
  
  #Read nonvar calls
  nonvars = read.table(nonvars.file, sep='\t', stringsAsFactors=FALSE)
  colnames(nonvars) = c('sample', 'chr', 'start', 'depth')
  nonvars[, 'sample'] = basename(nonvars[, 'sample'])
  nonvars[, 'sample'] = gsub('\\.nonvarsites\\.info', '', nonvars[, 'sample'])

  #add chr.start poskey for merging purposes
  chr.start = apply(nonvars[, c('chr', 'start')], 1, paste, collapse='.')
  chr.start = gsub(' ', '', chr.start)
  nonvars = cbind(nonvars, chr.start, stringsAsFactors=FALSE)

  #Set depth to retreived depth where missing depth.
  #chr.start not unique, ex: 1.100680607, so need a bit cumbersome solution (cant use rownames):
  varcalls.chr.start = varcalls[, 'chr.start']
  samples = unique(nonvars[, 'sample'])
  for(jsample in samples){
    jsample.ind = which(nonvars[, 'sample'] == jsample)
    jsample.keypos = nonvars[jsample.ind, 'chr.start']
    jsample.dp.col = paste(jsample, 'dp', sep='.')
    present.mat = nonvars[jsample.ind, c('chr.start', 'depth')]

    missing.keypos = setdiff(varcalls.chr.start, jsample.keypos)
    n.missing = length(missing.keypos)    
    missing.mat = cbind(missing.keypos, rep(NA, n.missing))
    colnames(missing.mat) = c('chr.start', 'depth')
    all.mat = rbind(present.mat, missing.mat)

    #add retreived depth column to varcalls
    a = merge(varcalls, all.mat, by = 'chr.start')

    #set depth to retreived depth where missing depth
    missing.depth.ind = which(a[, jsample.dp.col] == '-')
    retreived.depth.ind = which(!is.na(a[, ncol(a)]))
    add.depth.ind = intersect(missing.depth.ind, retreived.depth.ind)
    a[add.depth.ind, jsample.dp.col] = a[add.depth.ind, ncol(a)]

    #remove retreived depth column
    varcalls = a[, 1:(ncol(a) -1)]
  }

  return(varcalls)
}

read.varcalls <- function(varcall.file){

  varcalls = read.table(varcall.file, sep='\t', stringsAsFactors = FALSE)
  #rm last col (empty)
  varcalls = varcalls[, 1:(ncol(varcalls) -1)]

  first.colnames = c('chr', 'start', 'stop', 'refvar', 'observed var')
  last.colnames = c('gene relation', 'gene annot', 'var type', 'aa change', 'ConservedReg', 'SegDup', '1000g', 'Dbsnp_1', 'Dbsnp_2', 'cg46', 'Avsift', 'PP2')

  colnames(varcalls)[1:length(first.colnames)] = first.colnames
  colnames(varcalls)[(ncol(varcalls) - length(last.colnames) + 1):ncol(varcalls)] = last.colnames

  ###
  #Split each sample col into three columns (Hom/Het, AF, DP)
  ###
  sample.cols = (length(first.colnames) + 1):(ncol(varcalls) - length(last.colnames))

  n.samples = length(sample.cols)
  n.vars = nrow(varcalls)

  sample.data = list()
  length(sample.data) = n.samples
  sampledata.cols = c('het/hom', 'af', 'dp')
  n.sampledata.cols = length(sampledata.cols)
  for(jsample in 1:n.samples){

    sample.vars.mat = matrix(nrow = n.vars, ncol = n.sampledata.cols)
    colnames(sample.vars.mat) = sampledata.cols
    
    sample.col.ind = sample.cols[jsample]
    var.len = unlist(lapply(strsplit(as.character(varcalls[, sample.col.ind]), ';'), length))
    var.ind = which(var.len == n.sampledata.cols)

    var.data = t(as.data.frame(strsplit(as.character(varcalls[var.ind, sample.col.ind]), ';'), stringsAsFactors=FALSE)) #quite slow
    colnames(var.data) = sampledata.cols
      
    nonvar.ind = setdiff(1:n.vars, var.ind)
    nonvar.data = rep('-', length(nonvar.ind) * n.sampledata.cols)
    dim(nonvar.data) = c(length(nonvar.ind), n.sampledata.cols)
    colnames(nonvar.data) = sampledata.cols
    
    #get sampleID
    sample.id = sub('(.*)(Hom|Het)', '\\1', var.data[1, 1])
    
    #strip prefixes
    var.data[, 'het/hom'] = gsub('.*Hom$', 'Hom', var.data[, 'het/hom'])
    var.data[, 'het/hom'] = gsub('.*Het$', 'Het', var.data[, 'het/hom'])
    var.data[, 'af'] = as.numeric(gsub('(.*)=(.*)', '\\2', var.data[, 'af']))
    var.data[, 'dp'] = as.integer(gsub('(.*)=(.*)', '\\2', var.data[, 'dp']))
        
    sample.vars.mat[var.ind, ] = var.data
    sample.vars.mat[nonvar.ind, ] = nonvar.data

    #add sample data to list
    sample.data[[jsample]] = sample.vars.mat
    names(sample.data)[jsample] = sample.id    
  }
    
  allsamples.vars.mat = as.data.frame(sample.data, stringsAsFactors=FALSE)
  colnames(allsamples.vars.mat) = gsub('^X', '', colnames(allsamples.vars.mat))  
  varcalls = cbind(varcalls[, first.colnames], varcalls[, last.colnames], allsamples.vars.mat, stringsAsFactors = FALSE)

  #Rm first string from PP2 annot
  varcalls[, 'PP2'] = gsub('.*;', '', varcalls[, 'PP2'])
  
  #Rm first string from Sift annot
  varcalls[, 'Avsift'] = gsub('.*;', '', varcalls[, 'Avsift'])
  
  #Rm first string from 1000g
  varcalls[, '1000g'] = gsub('.*;', '', varcalls[, '1000g'])
  
  #Rm first string from cg46
  varcalls[, 'cg46'] = gsub('.*;', '', varcalls[, 'cg46'])

  #Rm 'chr' prefix from chromosome name
  varcalls[, 'chr'] = gsub('chr', '', varcalls[, 'chr'])

  #add chr.start poskey for merging purposes
  chr.start = apply(varcalls[, c('chr', 'start')], 1, paste, collapse='.')
  chr.start = gsub(' ', '', chr.start)
  varcalls = cbind(varcalls, chr.start, stringsAsFactors = FALSE)
  
  return(varcalls)
}

#Execute
main(varcalls.infile, varcalls.outfile)
