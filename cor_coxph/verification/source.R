quiet<-function(x,func,...){
  x2<-try(suppressWarnings(func(x,...)),silent = TRUE)
  if(inherits(x2,'try-error')){
    return(x)
  }else{
    return(x2)
  }
}

naturalize <- function(x, round_digits = NA){
  NAs <- sum(is.na(x)| x=="NaN")
  date_formats<-c("%Y/%m/%d","%Y-%m-%d","%d-%b-%Y")
  if(sum(is.na(quiet(x,as.numeric))) == NAs){
    x <- quiet(x,as.numeric)
    if(!is.na(round_digits)){
      x <- round(x,round_digits)
    }
    return(x)
  }else if(sum(is.na(quiet(x,as.Date,tryFormats = date_formats))) == NAs){
    return(quiet(x,as.Date,tryFormats = date_formats))
  }else{
    return(x)
  }
}

# Equal test to include NA comparisons as T/F
equal <- function(x1,x2, threshold){
  
  if (is.numeric(x1) & threshold != 0) {
    x1 <- abs(x1 - x2) <= (threshold * (abs(x1 + x2)/2))
    x2 <- rep(TRUE,length(x1))
  }
  
  c <- as.character(x1) == as.character(x2)
  c[ is.na(x1) | is.na(x2) ] <- (is.na(x1) == is.na(x2) )[ is.na(x1) | is.na(x2) ]
  c
  
}

compare_datasets <- function(cols,ds1,ds2,index, verbose = TRUE, round_digits = NA, thresholds = 0){
  res <- lapply(cols,function(column_name){
    
    if(is.list(round_digits)){
      if(column_name %in% names(round_digits)){
        round_digit <- round_digits[[column_name]]
      }else{
        round_digit <- round_digits[[".default"]]
      }
    }else{
      round_digit <- round_digits
    }
    
    if(is.list(thresholds)){
      if(column_name %in% names(thresholds)){
        threshold <- thresholds[[column_name]]
      }else{
        threshold <- thresholds[[".default"]]
      }
    }else{
      threshold <- thresholds
    }
    
    
    
    v1 <- naturalize( ds1[ order( ds1[[index]] ), column_name, drop = TRUE], round_digits = round_digit)
    v2 <- naturalize( ds2[ order( ds2[[index]] ), column_name, drop = TRUE], round_digits = round_digit)
    is_equal <- equal( v1, v2 , threshold = threshold)
    if ( !all( is_equal ) & verbose)
      print( paste0("`", column_name, "` is not equal. ", sum(!is_equal),"/", max(length(v1),length(v2)), " mismatches." ) )
    
    
    
    return( list(
      equal = all( is_equal ),
      diffs = data.frame(
        key = sort(ds1[[index]])[!is_equal],
        ds1 = v1[!is_equal],
        ds2 = v2[!is_equal],
        stringsAsFactors = FALSE
      ))
    )
  })
  
  names(res) <- cols
  
  resmatched <- do.call( 'c', lapply(res,`[[`,1))
  if( verbose ){
    message( "There are ", sum(!resmatched), " mismatched fields of ", length(resmatched),'.')
  }
  
  attr(res,"index") <- index
  
  class(res) <- "comparison"
  return( res )
}

print.comparison <- function(x,...){
  n_mismatch <- do.call('c',lapply(x,function(y){
    z <- nrow(y$diffs)
    if(z>0){
      z
    }else{
      NULL
    }
  }))
  
  for(comp in names(n_mismatch)){
    message( paste0("`", comp, "` is not equal. ", n_mismatch[[comp]], " mismatches." ) )
  }
}

rerun_failed_comparisons <- function( x, ds1, ds2){
  mismatched <- names(do.call('c',lapply(x,function(y){
    if(nrow(y$diffs)>0){
      1
    }else{
      NULL
    }
  })))
  compare_datasets(mismatched, ds1, ds2, index = attr(x,"index"))
}

reformat <- function(x, digits = 1){
  format(round(as.numeric(x),digits = digits), nsmall = digits)
}