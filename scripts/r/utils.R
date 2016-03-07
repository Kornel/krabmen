
filename.to.tumor <- function(filename) sub('.*gdac.broadinstitute.org_(\\w*)\\..*', '\\1', filename)

asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))

barcode.to.type <- function(patient.codes) return(patient.code.to.type(patient.codes))

barcode.to.participant <- function(barcodes) {
  sapply(barcodes, function(barcode) {
    elems <- unlist(strsplit(barcode, '-'))
    return(elems[3])
  })
}

se <- function(x) {
  sqrt(var(x) / length(x))
}


patient.code.to.type <- function(barcodes) {
  sapply(barcodes, function(barcode) {
    elems <- unlist(strsplit(barcode, '-'))
    sample.vial <- elems[4]
    sample <- as.numeric(substr(sample.vial, 1, 2))
    if (sample == 1) {
      return('tumor')
    } else if (sample == 11) {
      return('healthy')
    } else {
      return('other')
    }
  })
}



safe.f <- function(f) {
  function(x) {
    r <- f(x)
    if (r == 0) {
      return(0.0001)
    } else {
      return(r)
    }
  }
}

safe.mean <- safe.f(mean)
safe.median <- safe.f(median)

