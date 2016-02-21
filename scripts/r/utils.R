
asNumeric <- function(x) as.numeric(as.character(x))
factorsNumeric <- function(d) modifyList(d, lapply(d[, sapply(d, is.factor)], asNumeric))

patient.code.to.type <- function(patient.codes) {
  sapply(patient.codes, function(patient.code) {
    elems <- unlist(strsplit(patient.code, '-'))
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
