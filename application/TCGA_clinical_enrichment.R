# Clinical Enrichment analysis
# From NEMO paper
#clin = read.delim("~/TCGA_data/NAR Data/clinical/sarcoma", sep="\t")

get.fixed.names <- function(patient.names, include.type=F) {
  # fix the TCGA names to only include the patient ids
  return(patient.names)
  if (include.type) {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 15))))
  } else {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 12))))  
  }
}

get.clinical.params <- function(subtype.name) {
  LOC = "~/TCGA_data/NAR Data/clinical/"  
  clinical.data.path = paste(LOC, subtype.name, sep = '')
  clinical.params = read.table(clinical.data.path,
                               sep='\t', header=T, row.names = 1, stringsAsFactors = F)
  rownames.with.duplicates = get.fixed.names(rownames(clinical.params))  
  clinical.params = clinical.params[!duplicated(rownames.with.duplicates),]
  rownames(clinical.params) = rownames.with.duplicates[!duplicated(rownames.with.duplicates)]
  return(clinical.params)
}

check.clinical.enrichment <- function(clustering, subtype.name) {

  clinical.params = get.clinical.params(subtype.name)  
  
  clinical.metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC',
    pathologic_M='DISCRETE', pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE')
  
  pvalues = c()
  
  params.being.tested = c()
  
  for (clinical.param in names(clinical.metadata)) {
    
    if (!(clinical.param %in% colnames(clinical.params))) {
      #print(paste0('WARNING: ', clinical.param, ' does not appear for subtype ', subtype.name))
      next
    }
    
    clinical.values = clinical.params[names(clustering),clinical.param]
    is.discrete.param = clinical.metadata[clinical.param] == 'DISCRETE'
    is.numeric.param = clinical.metadata[clinical.param] == 'NUMERIC'
    stopifnot(is.discrete.param | is.numeric.param)
    
    # skip parameter if many missing values
    
    if (is.numeric.param) {
      numeric.entries = !is.na(as.numeric(clinical.values))
      if (4 * sum(numeric.entries) < length(clinical.values)) {
        #print(paste0('WARNING: skipping on ', clinical.param, ' for subtype ', subtype.name))
        next
      }
    } else {
      not.na.entries = !is.na(clinical.values)
      should.skip = F
      if (4 * sum(not.na.entries) < length(clinical.values)) {
        should.skip = T
      } else if (length(table(clinical.values[not.na.entries])) == 1) {
        should.skip = T
      }
      if (should.skip) {
        #print(paste0('WARNING: skipping on ', clinical.param, ' for subtype ', subtype.name))
        next
      }
    }
    
    params.being.tested = c(params.being.tested, clinical.param)
    
    if (is.discrete.param) {
      #clustering.with.clinical = cbind(clustering, clinical.values)
      #tbl = table(as.data.frame(clustering.with.clinical[!is.na(clinical.values),]))
      #test.res = chisq.test(tbl)
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[!is.na(clinical.values)], clinical.values[!is.na(clinical.values)], T)
      
    } else if (is.numeric.param) {
      #test.res = kruskal.test(as.numeric(clinical.values[numeric.entries]),
      #				clustering[numeric.entries])
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[numeric.entries], as.numeric(clinical.values[numeric.entries]), F)
    }
    
    pvalues = c(pvalues, pvalue)
    
  }
  names(pvalues) = params.being.tested
  return(pvalues)
}

get.empirical.clinical <- function(clustering, clinical.values, is.chisq) {
  #set.seed(42) # what ???
  if (is.chisq) {
      clustering.with.clinical = cbind(clustering, clinical.values)
      tbl = table(as.data.frame(clustering.with.clinical))
      test.res = chisq.test(tbl)
  } else {
    test.res = kruskal.test(as.numeric(clinical.values), clustering)
  }
  orig.pvalue = test.res$p.value
  num.iter = 1000
  total.num.iters = 0
  total.num.extreme = 0
  should.continue = T
  
  while (should.continue) {
    print('another iteration in empirical clinical')
    #perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
    perm.pvalues = as.numeric(lapply(1:num.iter, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
    
      if (is.chisq) {
        clustering.with.clinical = cbind(cur.clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
      } else {
        test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
      }
      cur.pvalue = test.res$p.value
      return(cur.pvalue)
    }))#, mc.cores=60))
    total.num.iters = total.num.iters + num.iter
    total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)
    
    binom.ret = binom.test(total.num.extreme, total.num.iters)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    sig.threshold = 0.05
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if (!is.threshold.in.conf | total.num.iters > 1e5) {
      should.continue = F
    }
  }
  return(cur.pvalue)
}
  