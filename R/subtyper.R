# Here's a good place to put your top-level package documentation

.onLoad <- function (lib, pkgname="subtyper") {
    ## Put stuff here you want to run when your package is loaded
    invisible()
}


#' Convert nrg format date to R format date
#'
#' @param x nrg format date
#' @return fixed x
#' @author Avants BB
#' @export
nrgDateToRDate <- function( x ) {
  for ( k in 1:length(x)) {
    x[k]=paste0( substr(x[k],1,4), "-", substr(x[k],5,6)  , "-", substr(x[k],7,8))
  }
  as.Date( x, format="%Y-%m-%d" )
}


#' Minimize the difference between two data frames based on t-statistic.
#'
#' This function generates random subsets of a data frame to minimize the difference with another data frame based on a specified set of columns, as measured by the t-statistic.  Authored by Avants and Chat-GPT 3.5.
#'
#' @param df1 Data frame to be subsetted.
#' @param df2 Data frame used as a reference for comparison.
#' @param cols Vector of column names used for matching.
#' @param sample_size the number to sample from df1
#' @param num_iterations Number of random subsets to generate.
#' @param restrict_df1 float lower quantile to restrict df1 based on first col value to match range of df2
#' @param option either random or optimal
#' @param verbose boolean
#'
#' @return rownames of a sub data frame that minimizes the difference with df2 in terms of t-statistic.
#'
#' @examples
#' set.seed(123)
#' df1 <- data.frame(A = rnorm(100), B = factor(sample(1:3, 100, replace = TRUE)), C = rnorm(100))
#' df2 <- data.frame(A = rnorm(50), B = factor(sample(1:3, 50, replace = TRUE)), C = rnorm(50))
#' matching_cols <- c("A", "B")
#' # matched_subset <- match_cohort_pair(df1, df2, matching_cols)
#' # print(matched_subset)
#'
#' @importFrom stats t.test
#' @importFrom optmatch pairmatch match_on
#' @export
match_cohort_pair <- function(df1, df2, cols, sample_size, num_iterations = 1000, restrict_df1=0.05, option='optimal', verbose=TRUE ) {

  stopifnot(  all( cols %in% colnames(df1) ) )
  stopifnot(  all( cols %in% colnames(df2) ) )
  
  if ( restrict_df1 > 0) {
      df1=df1[  
          df1[,cols[1]] > quantile(df2[,cols[1]],restrict_df1,na.rm=T) &
          df1[,cols[1]] < quantile(df2[,cols[1]],1.0-restrict_df1,na.rm=T),  ]
  }

  if ( missing( sample_size ) ) sample_size = min( c(nrow(df2),nrow(df1)))
  if ( sample_size > nrow(df1 ) ) sample_size = nrow( df1 ) - 1
  if ( option != 'optimal' ) {
    best_subset <- NULL
    min_t_statistic <- Inf
    
    # Convert factor columns to character for t-test
    for (col in cols) {
      if (!is.numeric(df1[,col])) {
        df1[,col] <- as.numeric(as.factor(df1[,col]))
      }
      if (!is.numeric(df2[,col])) {
        df2[,col] <- as.numeric(as.factor(df2[,col]))
      }
    }
    t_statistic=rep(NA,length(cols))
    names(t_statistic)=cols
    min_t_statistic=rep(Inf,length(cols))
    names(min_t_statistic)=cols
    for (i in 1:num_iterations) {
      # Randomly subset df1 based on the columns
      subset_indices <- sample(1:nrow(df1), 
        size = sample_size, replace = FALSE )
      # print( head( subset_indices ) )
      subset_df1 <- df1[subset_indices, cols]
      # Calculate t-statistic for the subset
      for (col in cols) {
          t_statistic[col]=abs(t.test(subset_df1[,col], df2[,col])$statistic)
          }
      
      # Check if the current subset has a smaller t-statistic
      if ( mean(t_statistic) < mean(min_t_statistic) ) {
        min_t_statistic <- t_statistic
        best_subset <- subset_indices
        if ( verbose ) print( paste( min_t_statistic ))
      }
    }
    if ( verbose ) print( paste( min_t_statistic ))
    return( rownames(df1)[best_subset] )
  } else {
    df1$my_match_cohort_pair_var = 0
    df2$my_match_cohort_pair_var = 1
    temp = rbind( df1, df2 )
    distances <- list()
    myform = paste0( "my_match_cohort_pair_var ~",  paste(cols,collapse="+") )
    distances$mahal <- match_on( as.formula(myform), data = temp)
    mahal.match <- pairmatch(distances$mahal, data = temp, remove.unmatchables = TRUE )
    matchednames = names( mahal.match )
    return( matchednames[ matchednames %in% rownames(df1) ] )
    return(  list( 
      df1=matchednames[ matchednames %in% rownames(df1) ], 
      df2=matchednames[ matchednames %in% rownames(df2) ] ) )
  }
}


#' Convert left/right variables to a measure of asymmetry
#'
#' @param mydataframe dataframe containing relevant variables
#' @param leftvar left side variable names ie the full names of the variables to asym
#' @param leftname the variable substring indicating left side
#' @param rightname the variable substring indicating right side
#' @param replacer string to replace left with in column names of output
#' @return fixed x
#' @author Avants BB
#' @export
mapAsymVar <-function( mydataframe, leftvar, leftname='left', rightname='right', replacer='Asym' ) {
  rightvar =  gsub( leftname, rightname, leftvar )
  hasright = rightvar %in% colnames(mydataframe)
  temp = mydataframe[,leftvar[hasright]] - mydataframe[,rightvar[hasright]]
  temp = temp * sign(temp )
  newnames = gsub(leftname, replacer,leftvar[hasright])
  mydataframe[,newnames]=temp
  return( mydataframe )
}



#' Convert left/right variables to an average measurement
#'
#' @param mydataframe dataframe containing relevant variables
#' @param leftvar left side variable names ie the full names of the variables to average
#' @param leftname the variable substring indicating left side
#' @param rightname the variable substring indicating right side
#' @param replacer string to replace left with in column names of output
#' @return fixed x
#' @author Avants BB
#' @export
mapLRAverageVar <- function( mydataframe, leftvar, leftname='left',rightname='right', replacer='LRAVG' ) {
  rightvar =  gsub( leftname, rightname, leftvar )
  hasright = rightvar %in% colnames(mydataframe)
  temp = mydataframe[,leftvar[hasright]] * 0.5 + mydataframe[,rightvar[hasright]] * 0.5
  newnames = gsub(leftname, replacer,leftvar[hasright])
  mydataframe[,newnames]=temp
  return( mydataframe )
}


#' Convert NA to false
#'
#' @param x a vector of bool
#' @return fixed x
#' @author Avants BB
#' @export
fs <- function( x ) {
  x[ is.na(x)]=FALSE
  x
}

#' Grep entries with a vector search parameters
#'
#' @param x a vector of search terms
#' @param desc target vector of items to be searched
#' @param intersect boolean whether to use intersection or union otherwise
#'
#' @return result of grep
#' @author Avants BB
#' @export
multigrep <- function( x, desc, intersect=FALSE ) {
  roisel = c()
  for ( xx in x ) {
    if (length(roisel)==0 | !intersect ) {
      roisel = c( roisel, grep(xx, desc) )
    } else {
      roisel = intersect( roisel, grep(xx, desc) )
    }
  }
  return(  roisel )
}



#' Extract column names with concatenated search parameters
#'
#' @param x vector of strings
#' @param demogIn the dataframe with column names to search.
#' @param exclusions the strings to exclude
#'
#' @return vector of string column names
#' @author Avants BB
#' @examples
#'
#' mydf = generateSubtyperData( 5 )
#' nms = getNamesFromDataframe( c("it","v"), mydf )
#'
#' @export
getNamesFromDataframe <- function( x, demogIn, exclusions ) {
  outnames = names(demogIn)[ grep(x[1],names(demogIn ) ) ]
  if ( length( x ) > 1 )
  for ( y in x[-1] )
    outnames = outnames[ grep(y,outnames ) ]

  if ( ! missing( exclusions ) ) {
    toexclude=grep(exclusions[1],outnames)
    if ( length(exclusions) > 1 )
      for ( zz in exclusions[-1] ) {
        toexclude = c( toexclude, grep(zz,outnames) )
      }
    if ( length( toexclude ) > 0 ) outnames = outnames[ -toexclude ]
  }
  return( outnames )
}


#' nallspdma
#'
#' Supplementary table S2 from Nalls 2019 Parkinsons disease meta-analysis of
#' genome wide association studies.
#'
#' @docType data
#' @references \url{https://pubmed.ncbi.nlm.nih.gov/31701892/}
#' @format A data frame (subset of columns described here):
#' \describe{
#'   \item{SNP}{the SNP of interest}
#'   \item{CHR}{associated chromosome}
#'   \item{BP}{base pair}
#'   \item{NearestGene}{nearest gene}
#' }
#' @usage data(nallspdma)
"nallspdma"


#' Generate example data for illustrating `subtyper`
#'
#' @param n numer of subjects
#' @param groupmeansbytime a 9-vector for each group by time mean for cognition.
#' the cognitive data is simulated from these mean values.  Group 3 changes
#' more rapidly, by default.
#'
#' @return data frame
#' @author Avants BB
#' @examples
#'
#' mydf = generateSubtyperData( 5 )
#'
#' @export
generateSubtyperData <-function( n = 100,
  groupmeansbytime = c(
      1, 3, 8, # V0 G0,G1,G2
      0.95, 3.5, 12,  # V1 G0,G1,G2
      1.05, 4.0, 20 )  # V2 G0,G1,G2
) {

  myids = LETTERS[1:24]
  groupn = round( length( myids )/3 )
  mygroup = c( rep("G0",groupn), rep("G1",groupn), rep("G2",groupn) )
  sampler = sample( 1:length(myids), n, replace=TRUE )
  mydf = data.frame(
    Id = myids[ sampler ],
    DX = mygroup[ sampler ],
    quality = rnorm( n ),
    RandomBasisProjection01 = rnorm( n ),
    RandomBasisProjection02 = rnorm( n ),
    RandomBasisProjection03 = rnorm( n ),
    RandomBasisProjection04 = rnorm( n ),
    RandomBasisProjection05 = rnorm( n ),
    RandomBasisProjection06 = rnorm( n ) )

  # now distribute visits over each subject
  visit = rep( NA, n )
  vizzes = c("V0","V1","V2")
  for ( u in myids ) {
    losel = mydf$Id == u
    loseln = sum( losel )
    loviz = sample( vizzes, loseln, replace=TRUE )
    visit[ losel ] = loviz
  }
  mydf$visit = visit

  # now distribute repeats over each subject
  repeati = rep( NA, n )
  for ( u in myids ) {
    for ( v in vizzes ) {
      losel = mydf$Id == u & mydf$visit == v
      loseln = sum( losel )
      lorep = paste0("Rep",formatC(1:loseln, width = 3, format = "d", flag = "0"))
      repeati[ losel ] = lorep
    }
  }
  mydf$repeater = repeati

  cognition = rep( NA, n )
  ct = 1
  for ( v in vizzes ) {
    for ( k in unique( mygroup ) ) {
      losel = mydf$DX == k & mydf$visit == v
      cognition[ losel ] = rnorm( sum(losel), mean = groupmeansbytime[ct] )
      ct = ct + 1
    }
  }
  mydf$cognition = cognition
  return( mydf )
}


#' Plot longitudinal trajectories across subtypes
#'
#' Uses `ggplot2` to visualize differences between subtypes over time.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param idvar variable name for unique subject identifier column
#' @param measurement the numeric variable to plot over time stratified by subtype
#' @param subtype the categorical (possibly ordinal) subtype variable
#' @param vizname the name of the grouped time variable (e.g. years change rounded to nearest quarter year)
#' @param whiskervar character either ci or se
#' @param extra additional naming variable
#' @param xlab additional naming variable
#' @param ylab additional naming variable
#'
#' @return ggplot
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' summ = plotSubtypeChange( mydf, "Id", "cognition", "DX", "visit" )
#' @export
#' @importFrom imbalance mwmote oversample racog rwo
#' @importFrom dCUR CUR
#' @importFrom MASS ginv
#' @importFrom dplyr arrange
#' @importFrom mclust Mclust predict.Mclust densityMclust
#' @importFrom visreg visreg
#' @importFrom dplyr sym bind_rows
#' @importFrom caret createDataPartition 
#' @importFrom ModTools OverSample UnderSample
#' @importFrom pheatmap pheatmap
#' @importFrom gprofiler2 gost gsnpense
#' @importFrom NMF nmf basismap
#' @importFrom stats lm predict qt rnorm var na.omit kmeans
#' @importFrom DDoutlier  LOOP  LOF  INFLO  RDOS  KDEOS  LDF  KNN_AGG  KNN_IN  KNN_SUM  RKOF
#' @importFrom ggplot2 aes ylim guides theme_bw scale_colour_hue geom_errorbar position_dodge element_text geom_line geom_point ggplot guide_legend
#' @importFrom ggplot2 xlab ylab theme rel geom_violin geom_boxplot scale_colour_manual
#' @importFrom plyr ddply rename
#' @importFrom grDevices dev.off pdf
plotSubtypeChange <-function( mxdfin,
                           idvar,
                           measurement,
                           subtype,
                           vizname='timer', # should be roundable to integer
                           whiskervar = c('ci','se'),
                           extra = "",
                           xlab = '',
                           ylab = ''
                          )
    {
      stopifnot( idvar %in% colnames(mxdfin) )
      stopifnot( measurement %in% colnames(mxdfin) )
      stopifnot( subtype %in% colnames(mxdfin) )
      stopifnot( vizname %in% colnames(mxdfin) )

      if ( ! all( whiskervar %in% c("ci","se") ) )
        stop("Must choose ci or se as a string for whiskervar")

      mysd <- function( x,na.rm ) sqrt(var(x,na.rm=na.rm))

      summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                            conf.interval=.95, .drop=TRUE) {
          # New version of length which can handle NA's: if na.rm==T, don't count them
          length2 <- function (x, na.rm=FALSE) {
              if (na.rm) sum(!is.na(x))
              else       length(x)
          }

          # This does the summary. For each group's data frame, return a vector with
          # N, mean, and sd
          datac <- plyr::ddply(data, groupvars, .drop=.drop,
            .fun = function(xx, col) {
              c(N    = length2(xx[[col]], na.rm=TRUE),
                mean = mean   (xx[[col]], na.rm=TRUE),
                sd   = mysd     (xx[[col]], na.rm=TRUE)
              )
            },
            measurevar
          )
          # Rename the "mean" column
          datac <- plyr::rename(datac, c("mean" = measurevar))

          datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

          # Confidence interval multiplier for standard error
          # Calculate t-statistic for confidence interval:
          # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
          ciMult <- qt(conf.interval/2 + .5, datac$N-1)
          datac$ci <- datac$se * ciMult

          return(datac)
      }


    lmxdf = mxdfin
    lmxdf$timer = lmxdf[ , vizname ]
# When the population standard deviation is known, the formula for a confidence interval (CI) for a population mean is x̄ ± z* σ/√n, where x̄ is the sample mean, σ is the population standard deviation, n is the sample size, and z* represents the appropriate z*-value from the standard normal distribution for your desired confidence level.

    tgcWithin <- summarySE( lmxdf, measurevar=measurement,
            groupvars=c('timer',subtype), na.rm=TRUE, conf.interval=.95 )

    colnames( tgcWithin )[1:5] = c('timer',"Quant","N","Changer","Changer_norm")
    tgcWithin$ymin = tgcWithin[,"Changer"] - tgcWithin[,whiskervar[1]]
    tgcWithin$ymax = tgcWithin[,"Changer"] + tgcWithin[,whiskervar[1]]
    yran = range( c( tgcWithin$ymin,tgcWithin$ymax ) )
    pd <- ggplot2::position_dodge( 0.05 ) # move them .05 to the left and right
    ww = 0.5
    namer=paste( subtype, ' vs ', measurement , extra )
    return( ggplot( tgcWithin,
      aes(x=timer, y=Changer,
        group = Quant, color = Quant  )) + #, shape = Quant
      geom_errorbar(aes(ymin=ymin, ymax=ymax), colour="black", width=ww, position=pd) +
      geom_line(lwd = 1, show.legend = FALSE ) +
      geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
      xlab( xlab ) + ylab( paste(ylab, measurement," +/-", whiskervar[1] )) +
      ylim(yran[1], yran[2]) +
      theme(strip.text = element_text(size = rel(1.1))) +
      theme(legend.position = 'top') +
      guides(fill=guide_legend(title.position="top")) +
      scale_colour_hue(name=subtype, l=40) + # Use darker colors, lightness=40
        theme_bw() +
        theme(text = element_text(size=20),
            axis.text.x = element_text(angle=0, vjust=1),
            legend.direction = "horizontal", legend.position = "top" )
      )

#    return( tgcWithin )
    }



#' Get subjects within a given set of visits
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param visitvar variable name for the visit column
#' @param visitselect the desired visits
#'
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydf = filterByVisit( mydf,  "visit", c("V0","V1"))
#' @export
filterByVisit  <-function(
  mxdfin,
  visitvar,
  visitselect ) {

return( mxdfin[ mxdfin[,visitvar] %in% visitselect, ] )

}




#' Get subjects and timepoints with the best quality
#'
#' This filters data that has both repeats (e.g. multiple images on the same day)
#' and (optionally) longitudinal data.  The quality measure should be a "higher
#' is better" measurement.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param idvar variable name for unique subject identifier column
#' @param visitvar variable name for the visit or date column
#' @param qualityvar variable name for the quality column; higher values should map to
#' higher quality data.
#'
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydfhq = highestQualityRepeat( mydf, "Id", "visit", "quality")
#' @export
highestQualityRepeat  <-function(
  mxdfin,
  idvar,
  visitvar,
  qualityvar ) {

  if ( ! ( visitvar %in% names( mxdfin ) ) ) stop("visitvar not in dataframe")
  if ( ! ( idvar %in% names( mxdfin ) ) ) stop("idvar not in dataframe")
  if ( ! ( qualityvar %in% names( mxdfin ) ) ) stop("qualityvar not in dataframe")
  vizzes = unique( mxdfin[,visitvar] )
  uids = unique( mxdfin[,idvar] )
  useit = rep( FALSE, nrow( mxdfin ) )
  for ( u in uids ) {
    losel = mxdfin[,idvar] == u
    vizzesloc = unique( mxdfin[ losel, visitvar ] )
    for ( v in vizzesloc ) {
      losel = mxdfin[,idvar] == u & mxdfin[,visitvar] == v
      mysnr = mxdfin[losel,qualityvar]
      myw = which( losel )
      if ( length( myw ) > 1 ) {
        if ( any( !is.na(mysnr) )  )
          useit[ myw[ which.max(mysnr) ] ] = TRUE
        } else useit[ myw ] = TRUE
      }
    }
  return( mxdfin[ useit, ] )
}


#' Reject subjects and timepoints with lowest quality
#'
#' This filters data that has both repeats (e.g. multiple images on the same day)
#' and (optionally) longitudinal data.  The quality measure should be a "higher
#' is better" measurement.  This function can be called recursively until it
#' converges.  It is less aggressive than \code{highestQualityRepeat}.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param idvar variable name for unique subject identifier column
#' @param visitvar variable name for the visit or date column
#' @param qualityvar variable name for the quality column; higher values should map to
#' higher quality data.
#' @param verbose boolean
#'
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydfhq = rejectLowestQualityRepeat( mydf, "Id", "visit", "quality")
#' @export
rejectLowestQualityRepeat <-function(
  mxdfin,
  idvar,
  visitvar,
  qualityvar,
  verbose = FALSE ) {

  if ( ! ( visitvar %in% names( mxdfin ) ) ) stop("visitvar not in dataframe")
  if ( ! ( idvar %in% names( mxdfin ) ) ) stop("idvar not in dataframe")
  if ( ! ( qualityvar %in% names( mxdfin ) ) ) stop("qualityvar not in dataframe")
  vizzes = sort(unique( mxdfin[,visitvar] ))
  uids = sort(unique( mxdfin[,idvar] ))
  useit = rep( TRUE, nrow( mxdfin ) )
  for ( u in uids ) {
    losel = mxdfin[,idvar] == u
    vizzesloc = unique( mxdfin[ losel, visitvar ] )
    for ( v in vizzesloc ) {
      losel = mxdfin[,idvar] == u & mxdfin[,visitvar] == v
      mysnr = mxdfin[losel,qualityvar]
      myw = which( losel )
      if ( length( myw ) > 1 )
        if ( any( !is.na(mysnr) )  )
          useit[ myw[ which.min(mysnr) ] ] = FALSE
      }
    }
  if ( verbose ) print( table( useit ) )
  return( mxdfin[ useit, ] )
}



#' Average repeat data
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param idvar variable name for unique subject identifier column
#' @param visitvar variable name for the visit or date column
#'
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydfhq = averageRepeats( mydf, "Id", "visit")
#' @export
averageRepeats  <-function(
  mxdfin,
  idvar,
  visitvar ) {

  if ( ! ( visitvar %in% names( mxdfin ) ) ) stop("visitvar not in dataframe")
  if ( ! ( idvar %in% names( mxdfin ) ) ) stop("idvar not in dataframe")
  vizzes = unique( mxdfin[,visitvar] )
  uids = unique( mxdfin[,idvar] )
  useit = rep( FALSE, nrow( mxdfin ) )
  num_cols <- unlist(lapply(mxdfin, is.numeric))
  for ( u in uids ) {
    losel = mxdfin[,idvar] == u
    vizzesloc = unique( mxdfin[ losel, visitvar ] )
    for ( v in vizzesloc ) {
      losel = mxdfin[,idvar] == u & mxdfin[,visitvar] == v
      myw = which( losel )
      if ( length( myw ) > 1 ) {
          useit[ myw[ 1 ] ] = TRUE
          mxdfin[ myw[ 1 ] , num_cols ] = colMeans( mxdfin[ myw, num_cols ], na.rm=TRUE )
        } else useit[ myw ] = TRUE
      }
    }
  return( mxdfin[ useit, ] )
}



#' Outlierness scoring
#'
#' Produce several complementary measurements of outlierness
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param measureColumns vector string of column names to use in outlier detection
#' @param calck optional integer for knn
#' @param outlierfunctions vector of strings naming outlier functions to report
#' @param verbose boolean
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mydf[8,rbfnames] = mydf[8,rbfnames] * 4.0
#' mydfol = outlierness( mydf, rbfnames )
#' olnames = names(mydf)[grep("OL_",names(mydf))]
#' mydf[6,olnames]
#' mydf[8,olnames]
#' @export
#' @importFrom stats quantile
#' @importFrom utils tail
outlierness <- function( mxdfin, measureColumns, calck,
  outlierfunctions = c("LOOP", "LOF",  "INFLO", "RDOS", "KDEOS",
    "LDF", "KNN_AGG", "KNN_IN", "KNN_SUM", "RKOF" ),
  verbose=FALSE ) {
  olnames = c("OL_LOOP", "OL_LOF", "OL_NOF", "OL_INFLO", "OL_RDOS", "OL_KDEOS",
    "OL_LDF", "OL_KNN_AGG", "OL_KNN_IN", "OL_KNN_SUM", "OL_RKOF" )
  oldata = scale( mxdfin[,measureColumns], TRUE, TRUE )
  if ( missing( calck ) )
    calck = DDoutlier::NAN( oldata )$r # this + LOOP is good
  for ( myolf in outlierfunctions ) {
    if ( myolf == "LOOP" ) {
      if ( verbose ) cat( "calc: LOOP ... ")
      mxdfin$OL_LOOP = DDoutlier::LOOP(  oldata, k=calck )
    } else if ( myolf == "LOF"  ) {
      if ( verbose ) cat( "calc: LOF ... ")
      mxdfin$OL_LOF = DDoutlier::LOF(  oldata, k=calck )
    } else if ( myolf == "INFLO"  ) {
      if ( verbose ) cat( "calc: INFLO ... ")
      mxdfin$OL_INFLO = ( DDoutlier::INFLO( oldata, k=calck ) )
    } else if ( myolf == "NOF" ) {
      if ( verbose ) cat( "calc: NOF ... ")
      mxdfin$OL_NOF = DDoutlier::NOF( oldata )$NOF
    } else if ( myolf == "KDEOS" ) {
      if ( verbose ) cat( "calc: KDEOS ... ")
      mxdfin$OL_KDEOS = DDoutlier::KDEOS( oldata, k_min = round( calck * 0.5 ), k_max = calck + 1 )
    } else if ( myolf == "LDF" ) {
      if ( verbose ) cat( "calc: LDF ... ")
      mxdfin$OL_LDF = DDoutlier::LDF( oldata, k = calck  )$LDF
    } else if ( myolf == "KNN_AGG" ) {
      if ( verbose ) cat( "calc: KNN_AGG ... ")
      mxdfin$OL_KNN_AGG = DDoutlier::KNN_AGG( oldata, k_min = round( calck * 0.5 ), k_max = calck + 1  )
    } else if ( myolf == "KNN_IN" ) {
      if ( verbose ) cat( "calc: KNN_IN ... ")
      mxdfin$OL_KNN_IN = -1.0 * DDoutlier::KNN_IN( oldata, k = calck )
    } else if ( myolf == "KNN_SUM" ) {
      if ( verbose ) cat( "calc: KNN_SUM ... ")
      mxdfin$OL_KNN_SUM = DDoutlier::KNN_SUM( oldata, k = calck  )
    } else if ( myolf == "RKOF" ) {
      if ( verbose ) cat( "calc: RKOF ... ")
      mxdfin$OL_RKOF = DDoutlier::RKOF( oldata, k = calck  )
    }
  }
  return( mxdfin )
}




#' Carefully replace a column name with a new one
#'
#' This will replace a column name with a new column name - it will not make
#' any changes to the underlying data.
#'
#' @param dataIn Input data frame with the old name within
#' @param oldName string column name to replace
#' @param newName the replacement name
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mydf = replaceName( mydf, rbfnames[1], 'mileyCyrus' )
#' @export
replaceName <- function( dataIn, oldName, newName ) {
  if ( oldName %in% names( dataIn ) )
    names(dataIn)[ names(dataIn) == oldName ] = newName
  return( dataIn )
}


#' Find test-retest data within a dataframe
#'
#' We assume there is a quality or related measure that will let us compute
#' distances between different time points and that these distances relate to
#' how similar two images are.  Note that the measurement may not relate to
#' quality at all but should map images into a similar metric space based on
#' appearance or some other objective quality that leads to the same value
#' when images are the same and changes continuously as images differ.
#'
#' @param dataIn Input data frame with the old name within
#' @param qualitycolumns the name(s) of quality-related column measurements
#' @param subjectID the unique subject id column name
#' @param visitID the column name that gives the visit / date; typically we look
#' for data the exists on the same visit.
#' @param whichVisit optional visit to use for the analysis
#' @param measureToRepeat the measurement to assemble for computing trt stats (e.g. ICC)
#' @param uniqueID optional column to contatenate to trt dataframe
#' @param covariates optional additional column name(s) to add to dataframe
#' @param nozeroes optional boolean - do not allow zero distance time
#' @return data frame with test-retest friendly organization ... trt0 and trt1 show the row indices of the test retest data
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mytrt = assembleTestRetest( mydf, rbfnames[1], 'Id', 'visit' )
#' @export
assembleTestRetest <- function(
  dataIn,
  qualitycolumns,
  subjectID,
  visitID,
  whichVisit,
  measureToRepeat,
  uniqueID, 
  covariates,
  nozeroes=FALSE ) {

  usubs = sort( unique( dataIn[,subjectID] ) )
  trtdf = data.frame( )
  for ( u in usubs ) {
    usel = dataIn[,subjectID] == u
    if ( ! missing( whichVisit ) ) {
      usel = dataIn[,subjectID] == u & dataIn[,visitID] == whichVisit
    }
    usel[ is.na( usel ) ] = FALSE
    if (sum(usel) > 1 ) {
      rowinds = which( usel )
      mydist = as.matrix( dist( dataIn[usel,qualitycolumns] ) )
      diag(mydist) = Inf
      if ( nozeroes ) mydist[ mydist == 0 ] = Inf
      mymin = min(mydist,na.rm=TRUE)
      bestind = which( mydist == mymin, arr.ind = T)[1,]
      n = nrow(trtdf) + 1
      trtdf[ n , subjectID ] = u
      if ( ! missing( whichVisit ) )
        trtdf[ n , visitID ] = whichVisit
      trtdf[ n , "distance" ] = mymin
      trtdf[ n , "trt0" ] = rowinds[ bestind[1] ]
      trtdf[ n , "trt1" ] = rowinds[ bestind[2] ]
      if ( ! missing( measureToRepeat )) {
        trtdf[ n , paste0(measureToRepeat,0) ] = dataIn[ trtdf[ n , "trt0" ], measureToRepeat ]
        trtdf[ n , paste0(measureToRepeat,1) ] = dataIn[ trtdf[ n , "trt1" ], measureToRepeat ]
      }
      if ( ! missing( uniqueID )) {
        trtdf[ n , paste0(uniqueID,0) ] = dataIn[ trtdf[ n , "trt0" ], uniqueID ]
        trtdf[ n , paste0(uniqueID,1) ] = dataIn[ trtdf[ n , "trt1" ], uniqueID ]
      }
      if ( ! missing( covariates ) ) {
        for (covariate in covariates) {
          trtdf[ n , paste0(covariate,"_0") ] = dataIn[ trtdf[ n , "trt0" ], covariate ]
          trtdf[ n , paste0(covariate,"_1") ] = dataIn[ trtdf[ n , "trt1" ], covariate ]
        }
      } 
    }
  }
  return( trtdf )
}


#' Uses outlierness scoring to filter a data frame
#'
#' This will produce an inclusion function based on scores genetrated by the
#' \code{outlinerness} function within this package.  We assume the data frame
#' has values of the type \code{RandBasisProj*} or \code{RandomBasisProj*} within its column names.
#'
#' @param dataIn Input data frame
#' @param bestolvar the outlier variable to prioritize
#' @param quantileThreshold the quantile to use for thresholing outlierness values
#' @param extraFiltering will perform brain-based application-specific filtering
#' @param calck optional parameter to pass to outlierness
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydfol = filterForGoodData( mydf, calck = 8 )
#' @export
#' @importFrom stats quantile
#' @importFrom utils tail
filterForGoodData <- function( dataIn,
    bestolvar = "OL_RKOF",
    quantileThreshold=0.99,
    extraFiltering = FALSE,
    calck = 256 ) {
  rbj = "RandBasisProj"
  if ( length( grep(rbj,names(dataIn)) ) == 0 ) {
    rbj = "RandomBasisProj"
    if ( length( grep(rbj,names(dataIn)) ) == 0 ) {
      stop("No RandBasisProj or RandomBasisProj column names. Cannot proceed.")
    }
  }
  epshvbrain_adjusted = 1.0 - quantileThreshold
  goodsel = rep( TRUE, nrow( dataIn ) )
  if ( ( "dkt_parc_wmSNR_dkt_parc_wmSNR" %in% names( dataIn ) ) & extraFiltering ) {
    qval = quantile( dataIn$dkt_parc_wmSNR_dkt_parc_wmSNR, epshvbrain_adjusted, na.rm=T  )
    goodsel = goodsel & dataIn$dkt_parc_wmSNR_dkt_parc_wmSNR > qval
  }
  if ( ( "GM_vol_tissues" %in% names( dataIn ) ) & extraFiltering ) {
    qval = quantile( dataIn$GM_vol_tissues, c(epshvbrain_adjusted, 1.0 - epshvbrain_adjusted), na.rm=T  )
    goodsel = goodsel & dataIn$GM_vol_tissues > qval[1]
  }

  rbvars = sort( names(dataIn)[grep(rbj,names(dataIn))] )
  goodsel = goodsel[  !is.na( dataIn[,rbvars[1]]) ]
  dataIn = dataIn[ !is.na( dataIn[,rbvars[1]]),  ]
  if ( length( grep("Pos",rbvars) ) > 0 )
    rbvars = rbvars[ -grep("Pos",rbvars)]
  dataIn = outlierness( dataIn, rbvars, calck=calck, verbose = TRUE)
  olvar = names(dataIn)[grep("OL_",names(dataIn)) ]
  goodsel = goodsel &
    dataIn[,bestolvar] < quantile( dataIn[,bestolvar], quantileThreshold, na.rm=T) &
    !is.na(   dataIn[,bestolvar])
  olvarinds = 1:length( olvar )
  olvarinds = olvarinds[ ! ( olvarinds %in% which( olvar == bestolvar ) ) ]
  for ( k in olvarinds ) {
    if ( ! is.na( dataIn[1,olvar[k]] ) ) {
      nextsel = dataIn[,olvar[k]] < quantile( dataIn[,olvar[k]], quantileThreshold, na.rm=T )
      goodsel = goodsel & nextsel
      }
    }
  goodsel[ is.na( goodsel )] = FALSE
  return( dataIn[ goodsel ,  ] )
}



#' fill baseline column
#'
#' build a column of data that maps baseline values to every row. baseline values
#' starting points or reference points for each subject, e.g. a value of a measurement
#' at time zero.  the method will produce a mean value if multiple entries match
#' the subjectID and visitID conditions.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param columnName string defining valid columns in the data frame.
#' @param subjectID the unique subject id column name
#' @param visitID names of the column that defines the variable defining
#' baseline-ness. e.g. \code{isBaseline}
#' @param baselineVisitValue the value defining baseline e.g. \code{TRUE}
#' @param baselineExt string appended to column name defining the baseline variable
#' @param deltaExt string appended to column name defining the change variable
#' @param fast boolean; will only return subjects with baseline values; if there 
#' are several baseline entries, these will be averaged. only works with numeric data.
#' @param verbose boolean
#' @return data frame with new columns
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' mydf = fillBaselineColumn( mydf, "RandomBasisProjection01", "Id", "visit", "V0" )
#' @export
fillBaselineColumn <- function(
   mxdfin,
   columnName,
   subjectID,
   visitID,
   baselineVisitValue,
   baselineExt = "_BL",
   deltaExt = "_delta",
   fast = TRUE,
   verbose=TRUE
) {

  myMode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
    }
  mostfreq <- function( x, na.rm=TRUE ) {
    if ( is.numeric( x ) ) return( mean(x,na.rm=na.rm ) )
    return( myMode( x ) )
  }
  if ( ! ( subjectID %in% names( mxdfin ) ) )
    stop("subjectID not in data frame's columns")
  if ( ! all( columnName %in% names( mxdfin ) ) )
    stop("columnName not in data frame's columns")
  if ( ! ( visitID %in% names( mxdfin ) ) )
    stop("visitID not in data frame's columns")
  if ( ! ( baselineVisitValue %in% mxdfin[,visitID] ) )
    stop("baselineVisitValue not in data frame's visitID")
  newcolname = paste0( columnName, baselineExt )
  newcolnamed = paste0( columnName, deltaExt )
  if ( fast ) { # use columnful merge (was use_data.table)
    # get the sid freq
    mxdfin[,subjectID]=as.character(mxdfin[,subjectID])
    mxdfin = mxdfin[ order( mxdfin[,subjectID] ), ]
    sidfreq = table( mxdfin[,subjectID])
    bldf = mxdfin[ fs(mxdfin[,visitID] == baselineVisitValue),  c(subjectID,columnName) ]
    sidtblbl = table( bldf[,subjectID]) 
    maxbln = max( sidtblbl )
    if ( maxbln > 1 ) {
      onlyones = names( sidtblbl )[ sidtblbl == 1 ]
      morethanones = names( sidtblbl )[ sidtblbl > 1 ]
      bldf1 = bldf[ bldf[,subjectID] %in% onlyones, ]
      sel2 = bldf[,subjectID] %in% morethanones
      if ( verbose )
        print( paste( "aggregate ", columnName, " in ", sum(sel2), "subjects with > 1 row ... these subjects have as many as", maxbln, "entries" ) )
      multisubs=bldf[sel2,subjectID]
      bldf2 = aggregate( bldf[sel2,c(subjectID,columnName)], list(multisubs), mostfreq, na.rm=T)
      bldf2[,subjectID]=unique( multisubs )
      bldf = base::rbind( bldf1[,c(subjectID,columnName)], bldf2[,c(subjectID,columnName)] )
      bldf = bldf[ order( bldf[,subjectID] ), ]
      if ( verbose )
        print("aggregation done")
      }
    inmcols = colnames( bldf ) %in% columnName
    colnames( bldf )[ inmcols ] = newcolname
    mxdfin = mxdfin[ mxdfin[,subjectID] %in% bldf[,subjectID], ]
    sidfreq = sidfreq[ names(sidfreq) %in% bldf[,subjectID] ]
    sidfreq = sidfreq[ as.character(bldf[,subjectID]) ]
    rownames(bldf)=bldf[,subjectID]
    if ( verbose ) {
      print( "begin repeat")
      print( as.integer(sidfreq) )
    }
    bldf = bldf[rep.int( bldf[,subjectID] , as.integer(sidfreq)), ]
    if ( verbose ) {
      print("end repeat/begin fill")
      print(dim(bldf))
      print(dim(mxdfin))
    }
    if ( verbose ) {
      print("selection done")
      print(dim(mxdfin))
      print(dim(bldf))
    }
    if ( identical( as.character(mxdfin[,subjectID]), as.character(bldf[,subjectID]) ) ) {
      mxdfin[,newcolname]=bldf[,newcolname]
    } else {
      print("Subject IDs are not identical")
      print(head(mxdfin[,subjectID]))
      print(head(bldf[,subjectID]))
      print(table( mxdfin[,subjectID]!=bldf[,subjectID]  ))
      stop("Subject IDs are not identical")
    }
    if ( verbose )
      print("end fill")
#    filldf = dplyr::bind_rows( mxdfin, bldf, .id=subjectID )
#    filldf = data.table(mxdfin)[bldf, on = subjectID, allow.cartesian=FALSE] # data.table
    if ( ! is.na( deltaExt ) ) {
      mxdfin[,newcolnamed] = mxdfin[,columnName] - mxdfin[,newcolname]
    }
    return( list(mxdfin, newcolname, newcolnamed ) )
  }
  mxdfin[,newcolname]=NA
  mxdfin[,newcolnamed]=NA
  visitidisnumeric = is.numeric(mxdfin[,visitID])
  usubs = unique( mxdfin[,subjectID] )
  nsubs = length( usubs )
  ct=0
  for ( u in usubs ) {
    if ( ct %% 20 == 0 ) cat( paste0(round(ct/nsubs*100),"%.") )
    ct=ct+1
    losel = fs( mxdfin[,subjectID] == u )
    lomxdfin = mxdfin[ losel ,  ]
    selbase = fs( lomxdfin[,visitID] == baselineVisitValue )
    if ( sum(selbase) == 0 & visitidisnumeric ) { # take next best value
      minval = min( lomxdfin[,visitID],na.rm=T )
      selbase = fs(lomxdfin[,visitID] == minval)
      }
    selbase[ is.na(selbase) ] = FALSE
    for ( jj in 1:length( columnName ) ) {
      columnNameLoc = columnName[jj]
      newcolnameLoc = newcolname[jj]
      isFactor = class( lomxdfin[,columnNameLoc] ) != "numeric"
      baseval = NA
      if ( sum( selbase ) > 0  & !isFactor ) {
        baseval = mean( lomxdfin[ selbase, columnNameLoc ],  na.rm=T )
      } else if ( sum( selbase ) > 0  & isFactor ) {
        baseval = median( lomxdfin[ selbase, columnNameLoc ],  na.rm=T )
      }
      mxdfin[ losel , newcolnameLoc ] = baseval
    }
  }
  if ( ! is.na( deltaExt ) )
    mxdfin[,newcolnamed] = mxdfin[,columnName] - mxdfin[,newcolname]
  return( list(mxdfin, newcolname, newcolnamed ) )
}


#' Covariate adjustment
#'
#' Adjust a training vector value by nuisance variables eg field strength etc.
#' One may want to use a specific sub-group for this, e.g. controls.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param adjustmentFormula string defining a valid formula to be used in \code{lm}.
#' @param groupVariable names of the column that defines the group to use for training.
#' @param group string defining a subgroup on which to train
#' @return data frame with adjusted measurement variable as defined in formula
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' myform = "cognition ~ RandomBasisProjection01 "
#' mydf = adjustByCovariates(mydf,myform,"DX","G0")
#' @export
adjustByCovariates  <- function(
  mxdfin,
  adjustmentFormula,
  groupVariable,
  group
) {
  outcomevar = gsub( " ", "", unlist(strsplit( adjustmentFormula, "~" ))[[1]] )
  if ( ! missing( group ) & ! missing( groupVariable ) ) {
    if ( ! (groupVariable %in% names( mxdfin ) ) ) stop("group name is wrong")
    gsel = mxdfin[,groupVariable] %in% group
    if ( sum(gsel) < 5 ) stop("too few subjects in subgroup for training")
    subdf = mxdfin[ gsel, ]
  } else subdf = mxdfin
  for ( xx in all.vars( as.formula( adjustmentFormula ) ) ) {
#    print( xx )
 #   print( table( is.na( mxdfin[,x] ) ) )
  }
  baseform = paste( outcomevar, " ~ 1" )
  imodel = lm( baseform, data=subdf ) # just the intercept
  ctlmodel = lm( adjustmentFormula, data=subdf )
#  print( table(subdf$istrain, subdf$studyName ) )
#  print( table(mxdfin$istrain, mxdfin$studyName ) )
  predintercept = predict( imodel, newdata = mxdfin  )
  predvol = predict( ctlmodel, newdata = mxdfin ) - predintercept
  adjustedoutcome = paste0( outcomevar, "_adjusted" )
  mxdfin[ , adjustedoutcome ] = mxdfin[,outcomevar] - predvol 
  return( mxdfin )
}


#' Covariate adjustment for many variables based on a single variable
#'
#' Adjust a training vector value by nuisance variables eg field strength etc.
#' One may want to use a specific sub-group for this, e.g. controls.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param adjustmentFormula string defining a valid formula to be used in \code{lm}.
#' @param columnstoadjust a vector of strings
#' @param groupVariable names of the column that defines the group to use for training.
#' @param group string defining a subgroup on which to train
#' @return data frame with adjusted measurement variable as defined in formula
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' myform = "cognition ~ RandomBasisProjection01 "
#' mydf = adjustByCovariates(mydf,myform,"DX","G0")
#' @export
adjustByCovariatesUni  <- function(
  mxdfin,
  adjustmentFormula,
  columnstoadjust,
  groupVariable,
  group
) {
  ##############################################
  cstrained <- function( x, istrain ) {
    xout = rep(NA,length(x))
    scaledTrainData = scale(x[istrain])
    xout[istrain] = scaledTrainData
    xout[!istrain] = scale(x[!istrain], 
        center=attr(scaledTrainData, "scaled:center"), 
        scale=attr(scaledTrainData, "scaled:scale"))
    return( xout )
  }
  outcomevar = gsub( " ", "", unlist(strsplit( adjustmentFormula, "~" ))[[1]] )
  if ( ! (groupVariable %in% names( mxdfin ) ) ) stop("group name is wrong")
  gsel = mxdfin[,groupVariable] %in% group
  if ( sum(gsel) < 5 ) stop("too few subjects in subgroup for training")
  subdf = mxdfin
  subdf$adjustByCovariatesUniTempVar = cstrained( 
    subdf[,outcomevar], gsel )
  adjustmentFormula2 = gsub(outcomevar, "adjustByCovariatesUniTempVar", adjustmentFormula )
  ctlmodel = lm( adjustmentFormula2, data=subdf[gsel,] )
  predvol = predict( ctlmodel, newdata = subdf )
  for ( zz in columnstoadjust ) {
    adjustedoutcome = paste0( zz, "_adjusted" )
    myscld = cstrained( mxdfin[,zz], gsel )
    mxdfin[ , adjustedoutcome ] = myscld - predvol
    }
  return( mxdfin )
}


#' Train subtype for univariate data
#'
#' This is the training module for subtype definition based on a vector.
#' Produces a data frame denoting the measurement variable, its cutpoints and
#' the names of the subtypes.
#' One may want to use a specific sub-group for this, e.g. patients.
#'
#' @param mxdfin Input data frame
#' @param measureColumn string defining a valid column to use for subtyping.  E.g.
#' the output of the covariate adjustment done by \code{adjustByCovariates}.
#' @param subtypes string vector naming the subtypes
#' @param quantiles numeric vector defining the quantiles that will yield subgroups
#' @param subtypename string naming the subtype variable
#' @param groupVariable names of the column that defines the group to use for training.
#' @param group string defining a subgroup on which to train
#' @return data frame defining the subtypes and cutpoints
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' qdf = trainSubtypeUni( mydf, "cognition", c("C0","C1","C2"), c(0.33,0.66) )
#' @export
trainSubtypeUni  <- function(
  mxdfin,
  measureColumn,
  subtypes,
  quantiles = c(0.25,0.75),
  subtypename = "subtype",
  groupVariable,
  group
) {

  if ( length( subtypes ) != ( length( quantiles ) + 1 ) )
    stop( "length( subtypes ) != ( length( quantiles ) + 1 )" )
  if ( ! missing( group ) & ! missing( groupVariable ) ) {
    if ( ! (groupVariable %in% names( mxdfin ) ) ) stop("group name is wrong")
    gsel = mxdfin[,groupVariable] %in% group
    if ( sum(gsel) < 5 ) stop("too few subjects in subgroup for training")
    subdf = mxdfin[ gsel, ]
  } else {
    subdf = mxdfin
    group=NA
  }

  nrow = length( quantiles ) + 1
  stdf = data.frame(
    subtypes = subtypes,
    measurement = rep( measureColumn, nrow ),
    quantiles = c( 0, quantiles ),
    quantileValues = rep( NA, nrow ),
    group = rep( group, nrow ) )
  names( stdf )[1] = subtypename
  stdf$quantileValues = c( -Inf, quantile( subdf[,measureColumn],quantiles, na.rm=TRUE ) )
  return( stdf )
}


#' Predict subtype for univariate data
#'
#' This is the inference module for subtype definition based on a vector.
#' After training, we can predict subtype very easily in new data based on
#' the data frame produced by training \code{trainSubtypeUni}.  If one passes
#' the visitName to the function then we will define subtype from the baseline
#' value alone.
#'
#' @param mxdfin Input data frame
#' @param subtypeDataFrame data frame defining the subtypes and cutpoints
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @param rename boolean will rename levels to user-provided names
#' @return data frame with attached subtypess
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mydf = outlierness( mydf, rbfnames )
#' mydf = highestQualityRepeat( mydf, "Id", "visit", "OL_KNN_SUM")
#' qdf = trainSubtypeUni( mydf, "cognition", c("C0","C1","C2"), c(0.33,0.66) )
#' pdf = predictSubtypeUni( mydf, qdf, "Id" )
#' @export
#' @importFrom Hmisc cut2
predictSubtypeUni  <- function(
  mxdfin,
  subtypeDataFrame,
  idvar,
  visitName,
  baselineVisit,
  rename = TRUE ) {

  msr = as.character( subtypeDataFrame[, "measurement"][1] )
  thesubtypes = subtypeDataFrame[,1]
  mxdfin[,names(subtypeDataFrame)[1]] = NA
  quantsV = as.numeric( subtypeDataFrame[,"quantileValues"][-1] )
  quantsP = as.numeric( subtypeDataFrame[,"quantiles"][-1] )
  # by default, just run across each row
  mxdfin[,names(subtypeDataFrame)[1]] = Hmisc::cut2( mxdfin[,msr], cuts = quantsV )
  theselevs = sort( na.omit( unique(  mxdfin[,names(subtypeDataFrame)[1]] ) ) )
  if ( rename ) {
    mxdfin[,names(subtypeDataFrame)[1]] = as.character( mxdfin[,names(subtypeDataFrame)[1]] )
    for ( k in 1:length( theselevs ) ) {
      losel = mxdfin[,names(subtypeDataFrame)[1]] == theselevs[k]
      mxdfin[,names(subtypeDataFrame)[1]][ losel ] = thesubtypes[k]
    }
    mxdfin[,names(subtypeDataFrame)[1]] <- factor(mxdfin[,names(subtypeDataFrame)[1]], levels=thesubtypes )
  } else {
    mxdfin[,names(subtypeDataFrame)[1]] <- factor(mxdfin[,names(subtypeDataFrame)[1]], levels=theselevs )
  }
  if ( missing( visitName ) | missing( baselineVisit ) )
    return( mxdfin )
  if ( ! ( baselineVisit %in% unique( mxdfin[,visitName] ) ) )
    stop( "! ( baselineVisit %in% unique( mxdfin[,visitName] ) )" )
  uids = unique( mxdfin[,idvar] )
  for ( u in uids ) {
    usel = mxdfin[,idvar] == u
    losel0 = mxdfin[,idvar] == u & mxdfin[,visitName] == baselineVisit & !is.na(mxdfin[,msr])
    losel0[ is.na( losel0 ) ] = FALSE
    if ( sum( losel0 ) > 0 ) {
      mytbl = table( mxdfin[losel0,names(subtypeDataFrame)[1]] )
      mostcommon = names( mytbl )[which.max(mytbl)]
      mxdfin[usel,names(subtypeDataFrame)[1]] = mostcommon
    }
  }
  return( mxdfin )

}





#' Train subtype for multivariate data
#'
#' This is the training module for subtype definition based on a matrix.
#' Currently only supports clustering based on gaussian mixtures
#' via the ClusterR package.
#' One may want to use a specific sub-group for this, e.g. patients.
#'
#' @param mxdfin Input data frame
#' @param measureColumns vector defining the data columns to be used for clustering.
#' Note that these methods may be sensitive to scaling so the user may want to
#' scale columns accordingly.
#' @param method string GMM or kmeans or medoid
#' @param desiredk number of subtypes
#' @param maxk maximum number of subtypes
#' @param groupVariable names of the column that defines the group to use for training.
#' @param group string defining a subgroup on which to train
#' @param frobNormThresh fractional value less than 1 indicating the amount of change
#' in the reconstruction error (measured by frobenius norm) from the previous
#' iteration 1 - F_cur / F_prev that will determine the optimal number of clusters.
#' For GMM clustering.
#' @param trainTestRatio Training testing split for finding optimal number
#' of clusters. For GMM clustering.  If zero, then will not split data. Otherwise,
#' will compute reconstruction error in test data only.
#' @param distance_metric see medoid methods in ClusterR
#' @param flexweights optional weights
#' @param flexgroup optional group
#' @param groupFun optional function name to use in group-guided clustering e.g. minSumClusters
#' @return the clustering object
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' gmmcl = trainSubtypeClusterMulti( mydf, rbfnames, maxk=4 )
#' @export
#' @importFrom mlr3cluster as_task_clust
#' @importFrom mlr3 mlr_learners
#' @importFrom ClusterR Cluster_Medoids predict_GMM Clara_Medoids GMM Optimal_Clusters_GMM KMeans_rcpp Optimal_Clusters_KMeans Optimal_Clusters_Medoids predict_Medoids
trainSubtypeClusterMulti  <- function(
  mxdfin,
  measureColumns,
  method = 'kmeans',
  desiredk,
  maxk,
  groupVariable,
  group,
  frobNormThresh = 0.01,
  trainTestRatio = 0,
  distance_metric=NULL,
  flexweights=NULL,
  flexgroup=NULL,
  groupFun = NULL
) {

  clustermeth = c("clara","fanny","pam")
  ktypes = c( clustermeth, "kmeans", 'kmeansflex', "GMM", "mclust", "pamCluster", 
    "kmeansflex","kmedians",  "angle",  "ejaccard", "flexcorr","flexcanb","flexmink", "cckmeans", "VarSelCluster",
    "hardcl","neuralgas", "hierarchicalCluster", "jaccard", mlr_learners$keys("clust") )

  if ( ! (method %in% ktypes ) ) {
    allmeth = paste( ktypes, collapse=' | ')
    message(paste("your chosen method:",method, "is not available."))
    message(paste( "try one of:" , allmeth ) )
    return( ktypes )
  }
  ismlr3=FALSE
  if ( length(grep( "clust[.]", method )) > 0 ) ismlr3=TRUE

  .env <- environment() ## identify the environment of cv.step
  if ( ! missing( group ) & ! missing( groupVariable ) ) {
    if ( ! (groupVariable %in% names( mxdfin ) ) ) stop("group name is wrong")
    gsel = mxdfin[,groupVariable] %in% group
    if ( sum(gsel) < 5 ) stop("too few subjects in subgroup for training")
    subdf = mxdfin[ gsel, measureColumns ]
  } else {
    subdf = mxdfin[ , measureColumns ]
    group=NA
  }

  if ( ismlr3 ) {
    task = as_task_clust(subdf)
    learner = mlr_learners$get(method)
    myparams = learner$param_set$ids()
    plist = list()
    if ( "k" %in% myparams  ) plist[["k"]]=desiredk
    if ( "num_clusters" %in% myparams  ) plist[["num_clusters"]]=desiredk
    if ( "centers" %in% myparams  ) plist[["centers"]]=desiredk
    learner$param_set$values = plist
    myp = learner$train(task)
    predict(myp,newdata=subdf)
    return( myp )
  }
  subdf = data.matrix( subdf )

  if ( method == "GMM" ) {
    if ( missing( desiredk ) ) {
      if ( missing( maxk ) ) maxk = 10
#      opt_gmm = ClusterR::Optimal_Clusters_GMM(subdf, max_clusters = maxk, criterion = "BIC",
#        dist_mode = "maha_dist", seed_mode = "random_subset",
#        km_iter = 10, em_iter = 10, var_floor = 1e-10,
#        plot_data = FALSE )
#      desiredk = which.min( opt_gmm )
        opt_gmm = rep( NA, maxk )
        opt_gmm_delta = Inf
        ii = 1
        if ( trainTestRatio == 0 ) {
          isTrain = rep( TRUE, nrow(subdf) )
          isTest = rep( TRUE, nrow(subdf) )
        } else {
          isTrainNum = sample( 1:nrow( subdf ), round( trainTestRatio * nrow( subdf ) ) )
          isTrain = rep( FALSE, nrow(subdf) )
          isTrain[ isTrainNum ] = TRUE
          isTest = !isTrain
        }
        while ( opt_gmm_delta > frobNormThresh & ii < maxk ) {
          gmm = ClusterR::GMM(subdf[isTrain,], gaussian_comps = ii,
              dist_mode = "maha_dist", seed_mode = "random_subset",
              km_iter = 10, em_iter = 5, verbose = FALSE)
          gmmclp = predictSubtypeClusterMulti( subdf[isTest,], measureColumns, gmm,
              clustername = "predictSubtypeClusterMultiClustersX" )
          if ( ii == 1 ) {
            myform = as.formula( "data.matrix(subdf[isTest,]) ~  1" , env = .env)
          } else {
            myform = as.formula( "data.matrix(subdf[isTest,]) ~  factor(gmmclp$predictSubtypeClusterMultiClustersX)" , env = .env)
          }
          mdl = lm( myform )
          opt_gmm[ii] = norm( data.matrix(predict(mdl,newdata=gmmclp)) - data.matrix(subdf[isTest,]), "F")
          if ( ii > 1 ) opt_gmm_delta = 1.0 - opt_gmm[ ii ]/opt_gmm[ ii - 1 ]
          ii = ii + 1
        }
        desiredk = ii
    }
    gmm = ClusterR::GMM( subdf,
        gaussian_comps = desiredk,
        dist_mode = "maha_dist",
        seed_mode = "random_subset",
        km_iter = 10,
        em_iter = 5, verbose = FALSE )
    return( gmm )
  }
  if ( method == "kmeans" ) {
    if ( missing( desiredk ) ) {
      if ( missing( maxk ) ) maxk = 10
      opt = ClusterR::Optimal_Clusters_KMeans( subdf, 
        max_clusters = maxk, plot_clusters = FALSE,
        criterion = 'distortion_fK', fK_threshold = 0.85,
        initializer = 'optimal_init', tol_optimal_init = 0.2)
      desiredk = which.min(opt)
      }
    km_rc = # stats::kmeans( subdf, desiredk )
      ClusterR::KMeans_rcpp(subdf,
        clusters = desiredk, num_init = 5, max_iters = 100,
        initializer = 'optimal_init', verbose = F, fuzzy=TRUE)
    return( km_rc )
  }
  if ( method == "VarSelCluster" ) { 
    # res_with <- VarSelCluster(x, 2, nbcores = 2, initModel=40, crit.varsel = "BIC")
    return( VarSelLCM::VarSelCluster(  subdf, desiredk, nbcores = 4, initModel=50, crit.varsel = "AIC", vbleSelec=FALSE ) )
  }
  if ( method == "mclust" ) { 
    return( mclust::Mclust(  subdf, desiredk ) )
  }
  if ( method == "medoid" ) {
    cl_X = data.matrix( subdf )
    for (i in 1:ncol(subdf)) { cl_X[, i] = as.numeric(cl_X[, i]) }
    if ( missing( desiredk ) ) {
      opt = Optimal_Clusters_Medoids(
          cl_X,
          max_clusters=maxk,
          distance_metric=distance_metric,
          criterion = "dissimilarity",
          clara_samples = 0,
          clara_sample_size = 0,
          minkowski_p = 1,
          swap_phase = TRUE,
          threads = 1,
          verbose = FALSE,
          plot_clusters = FALSE,
          seed = 1
        )
      return(opt)
    }
    cl_f = Cluster_Medoids(cl_X, clusters = desiredk, 
      distance_metric = distance_metric, minkowski_p = 1,
       threads = 1,
       swap_phase = TRUE,
       fuzzy = TRUE,
       verbose = FALSE )
    return( cl_f )
  }
  if ( method == "pamCluster" ) return( pamCluster(subdf,desiredk) )
  if ( method == "nmfCluster" ) return( nmfCluster(subdf,desiredk) )
  if ( method == "hierarchicalCluster" ) {
    if ( is.null(distance_metric) ) distance_metric='euclidean'
    return( hierarchicalCluster(subdf,desiredk,distmethod=distance_metric) )
  }
  if ( method == "EMCluster" ) return( EMCluster(subdf,desiredk) )
  if ( method == "FuzzyCluster" ) return( FuzzyCluster(subdf,desiredk) )
  flexmeth = c("kmeansflex", "flexkmeans", "kmedians", 
    "angle", "jaccard", "ejaccard","bootclust",
    "hardcl","neuralgas", "flexcorr","cckmeans",'flexmink','flexcanb')
  if ( method %in% clustermeth ) {
    if ( method == "clara" ) return( clara(subdf, desiredk ))
    if ( method == "fanny" ) return( fanny(subdf, desiredk ))
    if ( method == "pam" ) return( pam(subdf, desiredk ))
  }
  if ( method %in% flexmeth ) {
    initk = ClusterR::KMeans_rcpp(subdf,
        clusters = desiredk, num_init = 5, max_iters = 100,
        initializer = 'optimal_init', verbose = F, fuzzy=TRUE)$clusters
    if ( method %in% c('flexcorr','flexcanb','flexmink' ) ) {
      if ( method == 'flexcorr' ) ejacFam <- flexclust::kccaFamily(dist=distCor,cent=centMean)
      if ( method == 'flexcanb' ) ejacFam <- flexclust::kccaFamily(dist=distCanberra,cent=centMean)
      if ( method == 'flexmink' ) ejacFam <- flexclust::kccaFamily(dist=distMinkowski,cent=centMean)
      mycl = flexclust::kcca(
          subdf,
          k = initk,
          weights=flexweights, group=flexgroup,
        family = ejacFam )
      return( mycl )
    }
    if ( method %in% c("hardcl","neuralgas","cckmeans"))
      if ( method == 'cckmeans' ) {
        return(
          cclust(subdf, k=initk, weights=flexweights, group=flexgroup, method='kmeans' ) )
      } else return(
        cclust(subdf, k=initk, weights=flexweights, group=flexgroup, method=method ) )
    if ( method == "bootclust ")
      return( 
        bootFlexclust(subdf, k=2:desiredk, nboot=100, FUN=cclust) )
    if ( method == "kmeansflex" | method == "flexkmeans" ) method='kmeans'
      if ( !is.null(groupFun ) ) {
        #  c( "minSumClusters", "majorityClusters", "differentClusters" )
        myfam = kccaFamily(method, trim=0.005, groupFun =groupFun)
        } else {
        myfam = kccaFamily(method, trim=0.005)
        }
    return( 
      flexclust::kcca(
          subdf,
          k = initk,
          weights=flexweights, group=flexgroup,
        family = myfam ) )
  }

}



#' Bug fix to predict function in the VarSelLCM package
#'
#' @param object varselobject
#' @param newdata newdata for prediction
#' @export
VarSelLCMproba.post <- function(object, newdata){

  logprob <- matrix(object@param@pi, nrow(newdata), object@model@g, byrow=TRUE)
  for (nom in colnames(newdata)){
    xnotna <- newdata[,which(colnames(newdata)==nom)]
    where <- which(!is.na(xnotna))
    xnotna <- xnotna[where]
    if (nom %in% rownames(object@param@paramContinuous@mu)){
      who <- which(nom == rownames(object@param@paramContinuous@mu))
      for (k in 1:object@model@g) logprob[where,k] <- logprob[where,k] + dnorm(xnotna, object@param@paramContinuous@mu[who,k], object@param@paramContinuous@sd[who,k], log=TRUE)
    }else if (nom %in% rownames(object@param@paramInteger@lambda)){
      who <- which(nom == rownames(object@param@paramInteger@lambda))
      for (k in 1:object@model@g) logprob[where,k] <- logprob[where,k] + dpois(xnotna, object@param@paramInteger@lambda[who,k], log=TRUE)
    }else if (nom %in% names(object@param@paramCategorical@alpha))
      who <- which(nom ==  names(object@param@paramCategorical@alpha))
      if ( length(object@param@paramCategorical@alpha ) > 0 )
        for (k in 1:object@model@g){
          for (h in 1:ncol(object@param@paramCategorical@alpha[[who]]))
            logprob[where,k] <- logprob[where,k] + log(object@param@paramCategorical@alpha[[who]][k,h] ** (xnotna == colnames(object@param@paramCategorical@alpha[[who]])[h]))
      }
  }
  prob <- exp(logprob - apply(logprob, 1, max))
  prob/rowSums(prob)
}

#' Reorder a subtype variable based on an external reference vector
#'
#' @param mxdfin Input data frame
#' @param clustername column name for the identified clusters
#' @param reorderingVariable reorder the cluster names based on this variable
#' @return the new dataframe and (if present) membership variables
#' @author Avants BB
#' @export
reorderingDataframe <- function( mxdfin, clustername, reorderingVariable ) {
  xdf=aggregate( mxdfin[,reorderingVariable], list(mxdfin[,clustername]), mean, na.rm=TRUE )
  neword = order(xdf[,'x'],decreasing=FALSE)
  newname = as.character( xdf[neword,'Group.1'] )
  newvarval = xdf[neword,'x']
  xdf = cbind( xdf, newname, newvarval )
  names(xdf)=c('originalname','x', 'newname', paste0("ord_x"))
  return(xdf)
}

#' Predict subtype from multivariate data
#'
#' This is the inference module for subtype definition based on a matrix.
#' Currently only supports clustering based on gaussian mixtures
#' via the ClusterR package.
#'
#' @param mxdfin Input data frame
#' @param measureColumns vector defining the data columns to be used for clustering.
#' @param clusteringObject a clustering object to predict clusters
#' @param clustername column name for the identified clusters
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @param reorderingDataframe reorder the cluster names based on this dataframe mapping of original to new variable names
#' @param distance_metric see medoid methods in ClusterR
#' @return the clusters attached to the data frame; also returns membership probabilities
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' gmmcl = trainSubtypeClusterMulti( mydf, rbfnames, maxk=4 )
#' gmmclp = predictSubtypeClusterMulti( mydf, rbfnames, gmmcl )
#' @export
predictSubtypeClusterMulti  <- function(
  mxdfin,
  measureColumns,
  clusteringObject,
  clustername = 'GMMClusters',
  idvar,
  visitName,
  baselineVisit,
  reorderingDataframe,
  distance_metric = 'pearson_correlation'
) {
  myclustclass = class(clusteringObject)
  if ( length(myclustclass) == 1 )
    myclustclass=c(myclustclass,'other')
  subdf = mxdfin[ , measureColumns ]
  subdf = data.matrix( subdf )
  clustermeth = c("clara","fanny","pam")
  if ( myclustclass[1] %in% clustermeth ) {
    stop(paste("Prediction is not defined for ",myclustclass[1]))
  }
  if ( myclustclass[1] == "VSLCMresults" ) {
    knownnames = clusteringObject@data@var.names
    mypred = VarSelLCMproba.post( clusteringObject, subdf )
    classesare = apply( mypred, which.max, MARGIN=1 )
    mxdfin = cbind( mxdfin, factor( paste0(clustername,classesare ) ))
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
    cluster_memberships = data.frame(mypred)
    colnames(cluster_memberships) = paste0(clustername,"_mem_",1:ncol(cluster_memberships))
    mxdfin = cbind( mxdfin, cluster_memberships ) 
  } else if ( myclustclass[2] == "Gaussian Mixture Models" ) {
    pr = ClusterR::predict_GMM( subdf, clusteringObject$centroids,
      clusteringObject$covariance_matrices, clusteringObject$weights )
    mxdfin = cbind( mxdfin, factor( paste0(clustername,pr$cluster_labels ) ))
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
    cluster_memberships = data.frame(pr$cluster_proba)
    colnames(cluster_memberships) = paste0(clustername,"_mem_",1:nrow( clusteringObject$centroids))
    mxdfin = cbind( mxdfin, cluster_memberships )
  } else if (  myclustclass[2] == "k-means clustering" ) {
    # compute distance of every subject to each centroid
#    clusteringObject = stats::kmeans( subdf, desiredk )
    cluster_labels = rep( NA, nrow( subdf ) )
    cluster_memberships = matrix( nrow=nrow( subdf ), ncol=nrow( clusteringObject$centroids) )
    for ( kk in 1:nrow( subdf ) ) {
      dd = rep( NA, nrow( clusteringObject$centroids) )
      for ( jj in 1:nrow( clusteringObject$centroids) ) {
        dd[jj]=mean( ( as.numeric(subdf[kk,])-clusteringObject$centroids[jj,])^2 )
        cluster_memberships[kk,jj] = dd[jj]
        }
      cluster_labels[kk] = which.min(dd)
    }
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
    cluster_memberships = data.frame(cluster_memberships)
    colnames(cluster_memberships) = paste0(clustername,"_mem_",1:nrow( clusteringObject$centroids))
    mxdfin = cbind( mxdfin, cluster_memberships )
  } else if ( myclustclass[2] == "cluster medoids silhouette" ) { 
    cl_X = subdf
    for (i in 1:ncol(subdf)) { cl_X[, i] = as.numeric(cl_X[, i]) }
    pr = ClusterR::predict_Medoids( cl_X, clusteringObject$medoids,
      distance_metric =  distance_metric, fuzzy=TRUE )
    return( pr )
    cluster_labels = rep( NA, nrow( subdf ) )
    nclust = nrow(clusteringObject$medoids)
    cluster_memberships = matrix( nrow=nrow( subdf ), ncol=nclust )
    for ( kk in 1:nrow( subdf ) ) {
      dd = rep( NA, nclust )
      for ( jj in 1:nclust ) {
        dd[jj]=mean( ( as.numeric(subdf[kk,])-clusteringObject$centroids[jj,])^2 )
        cluster_memberships[kk,jj] = dd[jj]
        }
      cluster_labels[kk] = which.min(dd)
    }
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
    cluster_memberships = data.frame(cluster_memberships)
    colnames(cluster_memberships) = paste0(clustername,"_mem_",1:nrow( clusteringObject$centroids))
    mxdfin = cbind( mxdfin, cluster_memberships )
  } else if ( myclustclass[1] == "kcca" ) {
    cluster_labels = flexclust::predict( clusteringObject, newdata=subdf )
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
  } else if ( myclustclass[1] == "Mclust" ) {
    mypr = mclust::predict.Mclust( clusteringObject, newdata=subdf )
    cluster_labels = mypr$classification
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
    cluster_memberships = data.frame(mypr$z)
    colnames(cluster_memberships) = paste0(clustername,"_mem_",1:ncol(mypr$z))
    mxdfin = cbind( mxdfin, cluster_memberships )
  } else if ( myclustclass[1] %in% c("FuzzyCluster","pamCluster","EMCluster","hierarchicalCluster","nmfCluster") ) {
    mypr = predict( clusteringObject, newdata=subdf )
    cluster_labels =  paste0(clustername,as.character(mypr$classification))
    mxdfin = cbind( mxdfin, factor( cluster_labels ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
  } else if ( length(grep("Learner",myclustclass[1]))==1  ) {
    cluster_labels = predict( clusteringObject, newdata=data.frame(subdf) )
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
  } else stop("Unknown class of clustering object.")

  if ( ! missing( reorderingDataframe ) ) {
    # identify the mean value of the reo variable per class
    newclustername = rep(NA,nrow(mxdfin))
    for ( zz in 1:nrow(reorderingDataframe) ) {
      newclustername[  as.character(mxdfin[,clustername]) == as.character(reorderingDataframe[zz,'originalname']) ]=reorderingDataframe[zz,'newname']
    }
    mxdfin[,clustername]=as.character(mxdfin[,clustername])
#    mxdfin[,clustername]=factor(as.character(newclustername),levels=newclustername)
    mxdfin[,clustername]=newclustername
  }

  if ( missing( visitName ) | missing( baselineVisit ) )
    return( data.frame( mxdfin ) )

  msr = measureColumns[1]
  if ( ! ( baselineVisit %in% unique( mxdfin[,visitName] ) ) )
    stop( "! ( baselineVisit %in% unique( mxdfin[,visitName] ) )" )
  uids = unique( mxdfin[,idvar] )
  for ( u in uids ) {
    usel = mxdfin[,idvar] == u
    losel0 = mxdfin[,idvar] == u & mxdfin[,visitName] == baselineVisit & !is.na(mxdfin[,msr])
    losel0[ is.na( losel0 ) ] = FALSE
    if ( sum( losel0 ) > 0 ) {
      mytbl = table( mxdfin[losel0,clustername] )
      mostcommon = names( mytbl )[which.max(mytbl)]
      mxdfin[usel,clustername] = mostcommon
    }
  }
  return( data.frame( mxdfin ) )
}





#' Matrix factorization biclustering
#'
#' Cluster both rows and columns.  Return the joint cluster membership.
#'
#' @param mxdfin Input data frame
#' @param measureColumns vector defining the data columns to be used for clustering.
#' @param k the number of clusters
#' @param verbose boolean
#' @return a matrix with label identities; 0 indicates not in a cluster.
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mybic = biclusterMatrixFactorization( mydf, rbfnames, k = 2 )
#' @importFrom magrittr %>%
#' @importFrom ggpubr ggboxplot ggdotplot
#' @importFrom mlr3 lrn msr as_task_classif as_learner
#' @importFrom mlr3pipelines po selector_type %>>%
#' @importFrom fastICA fastICA
#' @importFrom data.table as.data.table
#' @importFrom mclust Mclust predict.Mclust mclustBIC
#' @importFrom fpc pamk
#' @importFrom flexclust kccaFamily kcca bootFlexclust cclust distCor centMean centAngle centMedian centOptim distJaccard distCanberra distAngle distEuclidean distManhattan distMax distMinkowski
#' @importFrom Evacluster pamCluster nmfCluster kmeansCluster hierarchicalCluster FuzzyCluster EMCluster 
#' @importFrom VarSelLCM VarSelCluster 
#' @importFrom Evacluster predict.pamCluster predict.nmfCluster predict.kmeansCluster predict.hierarchicalCluster predict.FuzzyCluster predict.EMCluster 
#' @export
biclusterMatrixFactorization  <- function(
  mxdfin,
  measureColumns,
  k = 2,
  verbose = FALSE
) {

  fixmat <- function( x ) {
    if ( is.null( rownames( x ) ) )
      rownames( x )=as.character( 1:nrow(x) )
    if ( is.null( colnames( x ) ) )
      colnames( x )=as.character( 1:ncol(x) )
    xt = t( x )
    zz=apply( xt, FUN=var, MARGIN=2)
    x = t( xt[, zz != 0 ] )
    zz=apply( x, FUN=var, MARGIN=2)
    myrm = rowMeans( x, na.rm=T )
    for ( rr in which( zz == 0 ) )
      x[,rr]=myrm
    return( x - min( x, na.rm=TRUE ) )
  }

  mybcmat = fixmat( data.matrix( mxdfin[,measureColumns]  ) )
  biclustmat = matrix( 0, nrow = nrow( mybcmat ), ncol = ncol( mybcmat ) )
  if ( is.null( rownames( mybcmat ) ) )
    rownames( mybcmat )=as.character( 1:nrow(mybcmat) )
  colnames( biclustmat ) = colnames( mybcmat )
  rownames( biclustmat ) = rownames( mybcmat )
  mynmf = NMF::nmf( mybcmat, k, method="Frobenius", seed='ica' )
#  w <- NMF::basis( mynmf )
#  h <- NMF::coef( mynmf )
  maxw = NMF::predict( mynmf, 'columns')
  maxh = NMF::predict( mynmf, 'rows')
  for ( j in 1:k ) {
    mycol = names( maxw )[ maxw == j ]
    myrow = names( maxh )[ maxh == j ]
    if ( verbose ) {
      print( k )
      print( mycol )
      print( myrow )
      }
    biclustmat[ rownames( biclustmat ) %in% myrow ,  colnames( biclustmat ) %in% mycol  ] = j
    }
 return( list( rowClusters=maxh, colClusters=maxw, jointClusters=biclustmat ) )
}






#' subtype feature importance
#'
#' Associate predictors (or features) with subtypes; these could be diagnoses
#' or cluster assignments.  Will use regression to return data frames intended
#' to visualize the most important relationships between features and types.
#' There are two approaches - worth using both.  Can be combined with
#' boostrapping to give distributional visualizations.
#'
#' @param dataframein Input dataframe with all relevant data
#' @param subtypeLabels Input subtype assignments.
#' @param featureMatrix matrix/dataframe defining the data columns as features.
#' @param associationType either predictor features from subtypes or predict
#' subtypes from features.  will produce related but complementary results. in
#' some cases, depending on subtypes/degrees of freedom, only one will be appropriate.
#' the third option (subjects) reports rownames of the dataframe that best fit the 
#' related subtype.
#' @param covariates optional string of covariates
#' @param transform optional effect_size
#' @param significance_level to threshold effects
#' @param visualize boolean
#' @return dataframes for visualization that show feature to subtype importance e.g. via \code{pheatmap}
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' fimp = featureImportanceForSubtypes( mydf, mydf$DX, mydf[,rbfnames], "subtypes2features" )
#' @importFrom fastICA fastICA
#' @importFrom mclust quantileMclust
#' @importFrom coca coca
#' @importFrom cluster fanny pam clara
#' @importFrom cluster fanny pam clara
#' @importFrom effectsize t_to_d z_to_d
#' @export
featureImportanceForSubtypes <- function(
    dataframein,
    subtypeLabels,   # cluster
    featureMatrix,   # featureColumns
    associationType = c( "features2subtypes", "subtypes2features", "subjects" ),
    covariates = "1",
    transform = 'effect_sizes',
    significance_level = 0.001,
    visualize = FALSE ) {
  stopifnot( associationType %in% c( "features2subtypes", "subtypes2features", "subjects" ) )
  form_pred<-function (object, ...)
  {
      if ( is.character(object)) object = as.formula(object)
      if (inherits(object, "formula")) {
          object <- terms(object)
      }
      y_index <- attr(object, "response")
      if (y_index != 0) {
          object[[2]] <- NULL
          object <- terms(object)
      }
      all.vars(object, ...)
  }
  if ( !is.factor(subtypeLabels)  )
    subtypeLabels=factor(subtypeLabels)
  if ( length( subtypeLabels ) != nrow( featureMatrix ) )
    stop("length( subtypeLabels ) != nrow( featureMatrix )")
  uniqClusts = levels( subtypeLabels )
  mync = length( uniqClusts )
  # stores feature importance
  clustzdescribe = data.frame( matrix( nrow = mync, ncol = ncol( featureMatrix ) ) )
  clustsigdescribe = data.frame( matrix( nrow = mync, ncol = ncol( featureMatrix ) ) )
  colnames( clustsigdescribe ) = colnames( featureMatrix )
  rownames( clustsigdescribe ) = uniqClusts
  colnames(clustzdescribe) = colnames( featureMatrix )
  clustmat = matrix( 0, nrow = length( subtypeLabels ), ncol=mync )
  colnames(clustmat) = as.character( uniqClusts )
  for ( j in 1:mync ) {
    losel = subtypeLabels == uniqClusts[j]
    if ( sum(losel) > 0 )
      clustmat[ losel , j] = 1
    }
  if ( associationType[1] == "subjects" ) {
    bestSubjectList = list()
    for ( j in 1:mync ) {
        mygt = gt( clustmat[,j],data.matrix(featureMatrix), 
          permutations=round(1.0/significance_level), standardize=TRUE, 
          model='logistic' )
        mysub = data.frame( subjects( mygt, what="z-score", sort=TRUE, cluster=FALSE ) )
        mysub[,'zscore']=( mysub[,'Residual'] ) / mysub[,'Std.dev']
        bestSubjectList[[j]]=mysub
        }
    return( bestSubjectList )
  } else if ( associationType[1] == "features2subtypes" ) {
    # converts cluster labels to one-hot coding
    for ( j in 1:mync ) {
      # return( list(clustmat,featureMatrix))
      myrform = paste(" temp ~ ", covariates )
      featureMatrixResid=featureMatrix
      for ( kk in 1:ncol(featureMatrix) ) {
        dataframein$temp = featureMatrix[,kk]
        featureMatrixResid[,kk] = residuals( lm( myrform, data=dataframein ) )
      }
      mygt = gt( clustmat[,j],data.matrix(featureMatrixResid), 
        permutations=round(1.0/significance_level), standardize=TRUE, 
        model='logistic'  )
      mycov = covariates( mygt, what="z-score", zoom=TRUE, cluster=FALSE, plot=FALSE )
      mycov = extract( mycov )
      mycoffs = data.frame( cbind( mycov@result, mycov@extra ) )
      mycoffssub=mycoffs[ mycoffs$direction == "assoc. with clustmat[, j] = 1" & 
        mycoffs$p.value <= 0.05, ]
      myz = mycoffs[,"Statistic"]
      myzsub = mycoffssub[,"Statistic"]
      if ( transform == 'effect_sizes' ) {
        myz = as.numeric( effectsize::z_to_d( myz, nrow(featureMatrix) )$d )
        myzsub = as.numeric( effectsize::z_to_d( myzsub, nrow(featureMatrix))$d  )
        }
      clustzdescribe[j,rownames(mycoffs)]=myz
      clustsigdescribe[j,rownames(mycoffssub)]=myzsub
    }
    if ( visualize ) pheatmap::pheatmap(abs(clustzdescribe),cluster_rows=F,cluster_cols=F)
    # get the max for each column
    clustsigdescribemax = clustsigdescribe * 0
    for ( j in 1:ncol(clustzdescribe) ) {
      wmax = which.max(abs(clustzdescribe[,j]))
      clustsigdescribemax[ wmax, j ] = clustzdescribe[wmax,j]
    }
    if ( visualize ) pheatmap::pheatmap((clustsigdescribemax),cluster_rows=F,cluster_cols=F)
    return( list(
      subtypeFeatureZScores = clustzdescribe,
      subtypeFeatureZScoresSignificant = clustsigdescribe,
      subtypeFeatureZScoresMax = clustsigdescribemax
    ) )
  } else {
    myform = as.formula(paste( "featureMatrix[,j] ~ clustmat[,k] + ", covariates ))
    if ( covariates != "1")
      locvars = form_pred( myform )[-c(1:2)]
    for ( j in 1:ncol(featureMatrix) ) {
      for( k in 1:mync ) {
        localdf = data.frame(
          feat=featureMatrix[,j],
          clust=clustmat[,k])
        if ( covariates != "1") {
          temper = dataframein[,locvars]
          if (length(locvars)==1) {
            temper=data.frame( temper )
            colnames(temper)=locvars
            }
          localdf=cbind(localdf, temper)
        }
        c1reg = lm( as.formula(myform), data=localdf )
        mycoffs = coefficients(summary(c1reg))
        if ( nrow(mycoffs ) > 1 ) {
          sigthresh = as.numeric( mycoffs[2,"Pr(>|t|)"] <= significance_level )
          myz = mycoffs[2,"t value"]
          if ( transform == 'effect_sizes' ) 
            myz = data.frame(effectsize::t_to_d( myz, nrow(localdf) ))[1,1]
          clustzdescribe[k,j]=myz
          clustsigdescribe[k,j]=myz * sigthresh
        } else {
          clustzdescribe[k,j]=NA
          clustsigdescribe[k,j]=NA
        }
      }
    }
    # get the max for each column
    clustsigdescribemax = clustsigdescribe * 0
    for ( j in 1:nrow(clustzdescribe) ) {
      wmax = which.max(abs(clustzdescribe[j,]))
      clustsigdescribemax[ j, wmax ] = clustzdescribe[j, wmax]
    }
    if ( visualize ) pheatmap::pheatmap((clustsigdescribemax),cluster_rows=F,cluster_cols=F)
    rownames(clustzdescribe)=rownames(clustsigdescribemax)
    return( list(
      subtypeFeatureTScores = clustzdescribe,
      subtypeFeatureTScoresSignificant = clustsigdescribe,
      subtypeFeatureTScoresMax = clustsigdescribemax
    ) )

  }
}



#' select top features via linear regression on an outcome
#'
#' @param dataframein Input dataframe with all relevant data
#' @param subtypeLabels Input subtype assignments.
#' @param featureMatrix matrix/dataframe defining the data columns as features.
#' @param covariates optional string of covariates
#' @param n_features select this many features per level
#' @param associationType either predictor features from subtypes or predict
#' subtypes from features.  will produce related but complementary results. in
#' some cases, depending on subtypes/degrees of freedom, only one will be appropriate.
#' @return vector of feature names
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' fimp = regressionBasedFeatureSelection( mydf, mydf$DX, mydf[,rbfnames], 
#'    associationType ="subtypes2features" )
#' @export
regressionBasedFeatureSelection <- function( 
    dataframein, subtypeLabels, featureMatrix, covariates="1", n_features=25, associationType ="features2subtypes" ) {
    fimp = featureImportanceForSubtypes(
        dataframein,
        subtypeLabels, 
        featureMatrix, 
        associationType[1], 
        covariates=covariates )
    if ( "subtypeFeatureTScoresSignificant" %in% names(fimp) ) {
      fimpname = "subtypeFeatureTScoresSignificant"
    } else if ( "subtypeFeatureZScoresSignificant" %in% names(fimp) ) {
      fimpname = "subtypeFeatureZScoresSignificant"
    }
    thefeats = c()
    myfimp = fimp[[fimpname]]
    for ( k in 1:nrow( myfimp )) {
        mm = abs(myfimp[k,])
        myor = order( as.numeric(mm[1,]), decreasing=TRUE )
        thefeats = c( thefeats, colnames( mm )[ head(myor,n_features) ] )
        }
    return( unique( thefeats ) )
    }





#' subtype plotting
#'
#' Assemble a set of standard plots looking at subtype results to support comparing
#' across a hierarchy of types both cross-sectionally and longitudinally.
#'
#' @param inputDataFrame Input complete data frame
#' @param variableToVisualize string naming the variable to display across subtypes
#' @param hierarchyOfSubtypes string vector of subtypes with increasing degrees of specificity
#' @param idvar variable name for unique subject identifier column
#' @param vizname the name of the grouped time variable (e.g. years change rounded to nearest quarter year)
#' @param whiskervar character either ci or se
#' @param consistentSubset display longitudinal data only from subjects that are consistently present at all visits
#' @param manualColors a list of user defined manual colors; the length of this list
#' should match the length of hierarchyOfSubtypes and colors should be named according
#' to the levels therein.  each entry in the list should be a string vector of color names.
#' @param outputPrefix filename prefix for the stored pdf plots; if missing, just plot to display
#' @param width the width of the graphics region in inches.
#' @param height the height of the graphics region in inches.
#' @return the output is a set of plots saved at the outputPrefix location
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 1000 )
#' qdf = trainSubtypeUni( mydf, "cognition", c("C0","C1","C2"), c(0.33,0.66) )
#' qdf = predictSubtypeUni( mydf, qdf, "Id" )
#' \dontrun{
#' hierarchicalSubtypePlots( qdf, "cognition", c("DX", "subtype" ),
#'  "Id", "visit", outputPrefix='/tmp/X' )
#' }
#' @importFrom globaltest gt covariates extract subjects
#' @importFrom ggstatsplot ggbetweenstats
#' @importFrom wesanderson wes_palette wes_palettes
#' @importFrom ggthemes theme_tufte
#' @importFrom magrittr %>%
#' @importFrom ggbeeswarm geom_beeswarm
#' @importFrom dplyr sym left_join mutate summarize n group_by
#' @importFrom ggplot2 labs element_text theme geom_smooth geom_violin geom_boxplot geom_dotplot ggtitle stat_smooth
#' @importFrom ggplot2 facet_grid facet_wrap label_both vars
#' @export
hierarchicalSubtypePlots <- function(
    inputDataFrame,
    variableToVisualize,
    hierarchyOfSubtypes,
    idvar,
    vizname,
    whiskervar=c('ci','se'),
    consistentSubset = FALSE,
    manualColors,
    outputPrefix,
    width=12, height=8 ) {
  if ( ! ( any(hierarchyOfSubtypes %in% names(inputDataFrame) ) ) )
    stop("some hierarchyOfSubtypes variables do not exist in the data frame.")
  if ( ! ( any(variableToVisualize %in% names(inputDataFrame) ) ) )
    stop(paste("the variableToVisualize", variableToVisualize, "does not exist in the data frame."))
  if ( ! missing( manualColors ) )
    stopifnot( length(manualColors) == length( hierarchyOfSubtypes ) )
  figs = c()
  ct = 1
  if ( missing( outputPrefix ) ) visualize = TRUE else visualize = FALSE

  # first level in hierarchy
  for (  k in 1:length( hierarchyOfSubtypes ) ) {
    losel = !is.na(inputDataFrame[,hierarchyOfSubtypes[k]])
    losel[ is.na(losel) ] = FALSE
    xplot0 <- ggstatsplot::ggbetweenstats(
          data = inputDataFrame[losel,],
          x = !!dplyr::sym( hierarchyOfSubtypes[k] ),
          y = !!dplyr::sym( variableToVisualize ),
          plot.type = "boxviolin", # type of plot
          bf.message = FALSE,
          results.subtitle = FALSE,
#          ggtheme = ggthemes::theme_tufte(), # ggplot2::theme_minimal(), # a different theme
          package = "wesanderson",
          palette = "Royal2"
        ) + ggplot2::labs( x = hierarchyOfSubtypes[k], y = variableToVisualize,
        title = paste("Subtype:",hierarchyOfSubtypes[k], "vs", variableToVisualize ) ) +
        ggplot2::theme(text = ggplot2::element_text(size = 20) ) +
        ggplot2::theme(plot.title = ggplot2::element_text(size=22))
    if ( ! missing( manualColors ) )
      xplot0 <- xplot0 + scale_colour_manual(values = manualColors[[k]] )
    if ( visualize ) print( xplot0 ) else {
      poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize, ".pdf" )
      outfn = paste0( outputPrefix, "_", poster )
      pdf( outfn, width=width, height=height )
      print( xplot0 )
      dev.off()
      figs[ct] = outfn
      ct = ct + 1
      }
    }

  # 2nd level in hierarchy - do level 1 vs other levels
  # e.g. if 1st level is DX, then plot subtypes within each DX
  if ( is.factor( inputDataFrame[,hierarchyOfSubtypes[1]] ) ) {
    unqDX = levels( inputDataFrame[,hierarchyOfSubtypes[1]] )
  } else unqDX = sort( unique( inputDataFrame[,hierarchyOfSubtypes[1]] ) )
  if ( length( hierarchyOfSubtypes ) > 1 ) {
    for (  k in 2:length( hierarchyOfSubtypes ) ) {
      for ( j in 1:length( unqDX ) ) {
        losel = inputDataFrame[,hierarchyOfSubtypes[1]] == unqDX[j] &
          !is.na(inputDataFrame[,hierarchyOfSubtypes[k]])
        losel[ is.na(losel) ] = FALSE
        xplot1 <- ggstatsplot::ggbetweenstats(
              data = inputDataFrame[ losel, ],
              x = !!dplyr::sym( hierarchyOfSubtypes[k] ),
              y = !!dplyr::sym( variableToVisualize ),
              plot.type = "boxviolin", # type of plot
              bf.message = FALSE,
              results.subtitle = FALSE,
#              ggtheme = ggthemes::theme_tufte(), # ggplot2::theme_minimal(), # a different theme
              package = "wesanderson",
              palette = "Royal2"
            ) + ggplot2::labs( x = hierarchyOfSubtypes[k], y = variableToVisualize,
            title = paste("Subtype:",hierarchyOfSubtypes[k], "vs", variableToVisualize, "@", unqDX[j] ) ) +
            ggplot2::theme(text = ggplot2::element_text(size = 20) ) +
            ggplot2::theme(plot.title = ggplot2::element_text(size=22))
        if ( ! missing( manualColors ) )
          xplot1 <- xplot1 + scale_colour_manual(values = manualColors[[k]] )
        if ( visualize ) print( xplot1 ) else {
          poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize, "_at_", toString(unqDX[j]),".pdf" )
          outfn = paste0( outputPrefix, "_", poster )
          pdf( outfn, width=width, height=height )
          print( xplot1 )
          dev.off()
          figs[ct] = outfn
          ct = ct + 1
          }
        }
      }
  }

  # for facetting
  wrap_by <- function(...) {
        facet_wrap(vars(...), labeller = label_both)
      }

  # now do something similar with longitudinal data
  if ( ! missing( vizname ) ) {
    uviz = sort( unique(  inputDataFrame[,vizname] ) )
    if ( length( uviz ) == 1 ) {
#      message("length( uviz ) == 1")
#      return( figs )
    }
    mysubs = unique( inputDataFrame[,idvar] )
    if ( consistentSubset ) {
      mysubs = unique( inputDataFrame[ inputDataFrame[,vizname] == uviz[1],idvar] )
      for ( jj in 2:length(uviz) ) {
        mysubs = intersect( mysubs,
          unique( inputDataFrame[ inputDataFrame[,vizname] == uviz[jj],idvar] ) )
        }
#      print( paste( "length(mysubs) = ", length(mysubs)) )
      }
    # do first level plots
    for (  k in 1:length( hierarchyOfSubtypes ) ) {
      myxlab = paste( hierarchyOfSubtypes[k], "vs", variableToVisualize,
        "longitudinal" )
      losel = !is.na(inputDataFrame[,hierarchyOfSubtypes[k]]) &
        inputDataFrame[,idvar] %in% mysubs
      losel[ is.na(losel) ] = FALSE
      lplot0 = plotSubtypeChange(
        inputDataFrame[losel,],
        idvar = idvar,
        measurement = variableToVisualize,
        subtype = hierarchyOfSubtypes[k],
        vizname = vizname, xlab = myxlab,
        whiskervar = whiskervar[1] )
      if ( ! missing( manualColors ) )
        lplot0 <- lplot0 + scale_colour_manual(values = manualColors[[k]] )
      if ( visualize ) print( lplot0 ) else {
        poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize, "_longitudinal.pdf" )
        outfn = paste0( outputPrefix, "_", poster )
        pdf( outfn, width=width, height=height )
        print( lplot0 )
        dev.off()
        figs[ct] = outfn
        ct = ct + 1
        }



      sym2 = sym(hierarchyOfSubtypes[k])
      ttt = inputDataFrame[losel,]
      viznameB=paste0(vizname,'_c')
      ttt[,viznameB] = as.numeric( as.factor( ttt[,vizname] ) )
      lplot1 =  ggplot( ttt,
            aes(
              y=!!sym(variableToVisualize), x=!!sym(viznameB),
                fill = (!!sym2), color =  (!!sym2)
            ) ) + wrap_by(!!sym2) +
                theme_bw() +
                geom_beeswarm() +
                geom_smooth(method = "lm", alpha = .15, aes(fill = (!!sym2)) ) +
                ggtitle( myxlab ) +
                theme(text = element_text(size=20),legend.position='top')
        if ( ! missing( manualColors ) )
          lplot1 <- lplot1 + scale_colour_manual(values = manualColors[[k]] )
        if ( visualize ) print( lplot1 ) else {
          poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize,
            "_longitudinal_swarmplot.pdf" )
          outfn = paste0( outputPrefix, "_", poster )
          pdf( outfn, width=width, height=height )
          print( lplot1 )
          dev.off()
          figs[ct] = outfn
          ct = ct + 1
          }
      }


    if ( length( hierarchyOfSubtypes ) > 1 ) {
      for (  k in 2:length( hierarchyOfSubtypes ) ) {
        for ( j in 1:length( unqDX ) ) {
          myxlab = paste( hierarchyOfSubtypes[k], "vs", variableToVisualize,
            "at", toString(unqDX[j]),
            "longitudinal" )
          losel = inputDataFrame[,hierarchyOfSubtypes[1]] == unqDX[j] &
            !is.na(inputDataFrame[,hierarchyOfSubtypes[k]]) &
              inputDataFrame[,idvar] %in% mysubs
          losel[ is.na(losel) ] = FALSE
          lplot1 = plotSubtypeChange(
            inputDataFrame[losel,],
            idvar = idvar,
            measurement = variableToVisualize,
            subtype = hierarchyOfSubtypes[k],
            vizname = vizname,
            whiskervar = whiskervar[1],
            xlab = myxlab )
          if ( ! missing( manualColors ) )
            lplot1 <- lplot1 + scale_colour_manual(values = manualColors[[k]] )
          if ( visualize ) print( lplot1 ) else {
            poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize,
              "_at_", toString(unqDX[j]),
              "_longitudinal.pdf" )
            outfn = paste0( outputPrefix, "_", poster )
            pdf( outfn, width=width, height=height )
            print( lplot1 )
            dev.off()
            figs[ct] = outfn
            ct = ct + 1
            }

          sym2 = sym(hierarchyOfSubtypes[k])
          sample_size = inputDataFrame[losel,] %>% group_by(!!sym(vizname)) %>% summarize(num=n())
          ttt = inputDataFrame[losel,]
          viznameB=paste0(vizname,'_c')
          ttt[,viznameB] = as.numeric( as.factor( ttt[,vizname] ) )
          lplot1 =  ggplot( ttt,
              aes(
                y=!!sym(variableToVisualize), x=!!sym(viznameB),
                  fill = (!!sym2), color =  (!!sym2)
#                  group=interaction(!!sym(viznameB), !!sym2)
              ) ) + wrap_by(!!sym2) +
                  theme_bw() +
                  geom_beeswarm() + # stat_smooth(method='lm')
                  geom_smooth(method = "lm", alpha = .15, aes(fill = (!!sym2)) ) +
                  ggtitle( myxlab ) +
                  theme(text = element_text(size=20),legend.position='top')
          if ( ! missing( manualColors ) )
            lplot1 <- lplot1 + scale_colour_manual(values = manualColors[[k]] )
          if ( visualize ) print( lplot1 ) else {
            poster = paste0( hierarchyOfSubtypes[k], "_vs_", variableToVisualize,
              "_at_", toString(unqDX[j]),
              "_longitudinal_swarmplot.pdf" )
            outfn = paste0( outputPrefix, "_", poster )
            pdf( outfn, width=width, height=height )
            print( lplot1 )
            dev.off()
            figs[ct] = outfn
            ct = ct + 1
            }
          }
        }
      }


    }
  return( figs )

}




#' genetic variants data frame from plink data
#'
#' @param rootFileName root for pgen psam and pvar files
#' @param targetSNPs snps to extract (optional - if absent, will get all)
#' @param verbose boolean
#' @return dataframes with both variants and subject ids
#' @author Avants BB
#' @examples
#' # mydf = plinkVariantsDataFrame( fn, c( 'rs6469804', 'rs6859' ) )
#' @importFrom pgenlibr NewPvar NewPgen ReadList
#' @importFrom data.table fread
#' @importFrom gaston as.matrix read.bed.matrix
#' @export
plinkVariantsDataFrame <- function( rootFileName, targetSNPs,  verbose=FALSE ) {
  type=getExt(rootFileName)
  stopifnot( type %in% c("pgen","bed") )
  if ( type == 'pgen' ) {
    rootFileName2 = tools::file_path_sans_ext( rootFileName )
    f.pvar = paste0(rootFileName2, '.pvar')
    f.pgen = paste0(rootFileName2, '.pgen')
    f.sam = paste0(rootFileName2, '.psam')
    subjectIDs = fread( f.sam )
    myvariants = fread( f.pvar, select = 3)
    if ( missing(targetSNPs) ) i = 1:length(myvariants$ID)
    if ( !missing(targetSNPs) ) i  <- which( myvariants$ID %in% targetSNPs)
    pvar <- pgenlibr::NewPvar(f.pvar)
    pgen <- pgenlibr::NewPgen(f.pgen, pvar=pvar) #,sample_subset = c(1,2,3,4) )
    myderk = pgenlibr::ReadList( pgen, i, meanimpute=F )
    snpdf = data.frame( myderk )
    colnames( snpdf ) = myvariants$ID[i]
    outdf = data.frame( subjectIDs = subjectIDs$IID, snpdf )
    return( outdf )
  } else {
    gwas=gaston::read.bed.matrix( rootFileName )
    ww=which( gwas@snps$id %in% targetSNPs )
    if ( length( ww ) == 0 ) {
      mymsg = paste("No target SNPs are present ")
      print(mymsg)
      message(mymsg)
      return(NA)
      }
    y=gwas[,ww]
    snpnames=y@snps$id
    y=gaston::as.matrix( y )
    gmydf=data.frame(subjectIDs=rownames(y))
    gmydf=cbind(gmydf,y)
    return(gmydf)
  }
}




#' three way interaction plot from raw data
#'
#' @param indf dataframe with relevant variables
#' @param xvar variable for x-axis of plots
#' @param yvar variable for y-axis of plots
#' @param colorvar variable by which to color plots
#' @param anat continuous variable by which to split plots
#' @param anatshow optional character name for continuous variable to show on plots - can be up to vector length 3 - one for each panel of the plot
#' @param ggpalette optional palette
#' @param showpoints boolean
#' @return the plot
#' @author Avants BB
#' @examples
#' # FIXME
#' @importFrom ggpubr ggscatter set_palette
#' @importFrom gridExtra grid.arrange
#' @export
threewayinteraction <- function( indf, xvar, yvar, colorvar, anat, anatshow, ggpalette='jco', showpoints=FALSE ) {
  glist = list()
  if ( ! missing(anatshow) ) {
    while ( length(anatshow) <  3 )
      anatshow=c(anatshow,"")
  }
  if ( missing(anatshow) ) {
    anatshow=gsub("T1Hier_","",anat)
    anatshow=rep(anatshow,3)
  }
  indf[,colorvar]=factor(indf[,colorvar])
  indf$snapfact=factor(indf[,colorvar])
  ylimmer = range( indf[,yvar])
  glist[[length(glist)+1]]= ggscatter(indf, x = xvar, y = yvar, color=colorvar,   size=3.45,  point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ggtitle(paste(anatshow[1])) + ylim( ylimmer ) #+ theme(legend.position = "none")

  if ( is.numeric(indf[,anat]) ) {
    medsplit = median( indf[,anat], na.rm=T )
    hisel = indf[,anat] > medsplit
    loclev=loclev2=anat
  } else {
    loclev = (unique(indf[,anat])[1])
    hisel=indf[,anat]==loclev
    loclev2=paste0("!",loclev)
  }
  glist[[length(glist)+1]]=ggscatter(indf[hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ggtitle(paste('+',anatshow[2],loclev)) + theme(legend.position = "none") + ylim( ylimmer )

  
  glist[[length(glist)+1]]=ggscatter(indf[!hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE  ) + theme(text = element_text(size=12))+ ggtitle(paste('-',anatshow[3],loclev2))+ theme(legend.position = "none") + ylim(  ylimmer )

  grid.arrange(grobs=glist,ncol=3)

}




#' Get the file extension from a file-name. from `reader` package.
#'
#' @param fn filename(s) (with full path is ok too)
#' @return returns the (usually) 3 character file extension of a filename
#' @export 
#' @author Nicholas Cooper \email{nick.cooper@@cimr.cam.ac.uk}
getExt <- function (fn) 
{
    if (length(fn) < 1) {
        warning("fn had length of zero")
        return(fn)
    }
    if (all(is.na(fn)) | !is.character(fn)) {
        stop("fn should not be NA and should be of type character()")
    }
    strip.file.frags <- function(X) {
        file.segs <- strsplit(X, ".", fixed = TRUE)[[1]]
        lss <- length(file.segs)
        if (lss > 1) {
            out <- paste(file.segs[lss])
        }
        else {
            out <- ""
        }
        return(out)
    }
    return(sapply(fn, strip.file.frags))
}





#' Match a pair of vector distributions based on quantiles 
#'
#' @param vecToBeTransformed input vector to be transformed to match the reference
#' @param vecReference reference vector
#' @param quantiles a vector of quantile points to match
#' @param polynomialOrder integer greater than or equal to one
#' @param truncate boolean
#' @return the transformed vector
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' newvec = matchVectorDistributionByQuantiles( mydf[,rbfnames[1]], mydf[,rbfnames[1]] )
#' @export
matchVectorDistributionByQuantiles  <- function(
  vecToBeTransformed,
  vecReference,
  quantiles = 1:9/10.0, 
  polynomialOrder = 1,
  truncate = TRUE ) {
    myq1 = quantile( vecToBeTransformed, quantiles, na.rm=T )
    myq2 = quantile( vecReference, quantiles, na.rm=T )
    mytd = data.frame( myq2=myq2, myq1=myq1 )
    mdl = lm( myq2 ~ stats::poly(myq1,polynomialOrder), data=mytd )
    mynd=data.frame( myq1 = vecToBeTransformed )
    outvec = predict( mdl, newdata=mynd  )
    if ( truncate ) {
      outvec[ outvec < min(vecReference,na.rm=T)]=min(vecReference,na.rm=T)
      outvec[ outvec > max(vecReference,na.rm=T)]=max(vecReference,na.rm=T)
    }
    return( outvec )
  }






#' quantify by quantiles (quantSquared)
#' 
#' transform input vector to a quantile representation relative to a reference. 
#'
#' @param vecToBeTransformed Input vector should have same entries as reference;
#' should be a row of a data frame
#' @param matrixReferenceDistribution data frame defining the reference data;  
#' ideally this will not contain missing data. 
#' @return the quantile transformed vector
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' zz=data.frame(t(colMeans(mydf[,rbfnames])))
#' newvec = quantSquared( zz, mydf )
#' newvec2 = quantSquared( mydf[1,rbfnames], mydf )
#' @export
quantSquared  <- function(
  vecToBeTransformed,
  matrixReferenceDistribution ) {
    mycnt = colnames( vecToBeTransformed )
    newvec = rep( NA, length( vecToBeTransformed ) )
    names(newvec)=mycnt
    for ( x in mycnt ) {
      refvec = matrixReferenceDistribution[,x]
      nna = sum( !is.na( refvec ) )
      boolvec =  as.numeric(vecToBeTransformed[x]) > as.numeric(refvec)
      newvec[x] = sum( boolvec , na.rm=TRUE )/ nna
    }
    return( newvec )
  }



#' prplot
#' 
#' partial residual regression plot using ggpubr for display 
#' and visreg for partial residual calculation 
#'
#' @param mdl input fitted model
#' @param xvariable a model term
#' @param byvariable a model term interacting with the xvariable
#' @param titlestring string for title
#' @param ystring string for title
#' @param addpoints continuous value greater than zero
#' @param palette string
#' @param colorvar string
#' @return the quantile transformed vector
#' @author Avants BB
#' @export
prplot  <- function(
  mdl, xvariable, byvariable, titlestring='', ystring='', addpoints=0, palette='npg', colorvar=''
   ) {
  addthepoints=FALSE
  if ( addpoints > 0 ) addthepoints=TRUE
  if ( ! missing( byvariable ) ) {
    vv=visreg::visreg( mdl, xvariable, by=byvariable, plot=FALSE)
    if ( is.factor(vv$res[,xvariable] ) | is.character(vv$res[,xvariable]) ) {
      return( ggdotplot(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    conf.int=T,
                    point=addthepoints,
                    fill=colorvar, facet.by=byvariable,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
    } else return( ggscatter(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    point=addthepoints, add='reg.line', conf.int=T,
                    color=colorvar, facet.by=byvariable,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
  }
  if ( missing( byvariable ) ) {
    vv=visreg::visreg( mdl, xvariable, plot=FALSE)
     if ( is.factor(vv$res[,xvariable] ) | is.character(vv$res[,xvariable]) ) {
      return( ggboxplot(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    conf.int=T,
                    point=addthepoints,
                    fill=colorvar,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
    } else return( ggscatter(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, 
                    point=addthepoints, add='reg.line', conf.int=T,
                    color=colorvar, palette=palette,
                    cor.coef=TRUE ) +
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
    }
  }


#' balanced sampling of a variable
#' 
#' resample a dataframe to counteract variable imbalance
#' 
#' @param x data frame to extend
#' @param variable column name
#' @param method one of permissible methods listed in function (will print if wrong method passed)
#' 
#' @return new data frame
#' @author Avants BB
#' @export
balancedDataframe <- function( x, variable, method ) {
    valbal = c("none", "rwo", "racog", "mwmote", "over", "under" )
    if ( ! ( method %in% valbal ) )
        stop(paste("method must be one of ", paste(valbal, collapse=" / ")))
    if ( method == 'none' ) return( x )
    temptbl = table( x[,variable] )
    myinstsize = max( temptbl ) - min( temptbl )
    if ( method == 'over') {
        return( ModTools::OverSample( x, variable ) )
    }
    if ( method == 'under') {
        return( ModTools::UnderSample( x, variable ) )
    }
    if ( method == 'mwmote')
        balDF = mwmote( dataset = x, numInstances = myinstsize, classAttr = variable)
    if ( method == 'rwo')
        balDF = rwo( dataset = x, numInstances = myinstsize, classAttr = variable)
    if ( method == 'racog')
        balDF = racog( dataset = x, numInstances = myinstsize, classAttr = variable)
    x = rbind( x, balDF )
    return(x)
    }



#' balanced sampling of a multi-class variable
#' 
#' resample a dataframe to counteract variable imbalance
#' 
#' @param x data frame to extend
#' @param variable column name
#' @param method one of permissible methods listed in function (will print if wrong method passed)
#' @param minimum_ratio minimum acceptable ratio of smallest to largest group n
#' 
#' @return new data frame
#' @author Avants BB
#' @export
balanceDataMultiClass <- function( x, variable, method, minimum_ratio=0.99 ) {
  if ( method == 'none' ) return( x )
  if ( method == 'under') {
      minmaxdifftbl = table( x[,variable] )
      mintbl = min( minmaxdifftbl )
      diffis = min( minmaxdifftbl ) / max( minmaxdifftbl )
      while( diffis < minimum_ratio ) {
          minmaxdifftbl = table( x[,variable] )
          maxtbl = max( minmaxdifftbl )
          mintbl = min( minmaxdifftbl )
          diffis = mintbl/maxtbl
          if ( diffis < minimum_ratio ) {
              maxname = names( minmaxdifftbl[ minmaxdifftbl == max(minmaxdifftbl) ] )
              samplethismany = round( mintbl + ( maxtbl - mintbl ) * minimum_ratio ) - 1
              indsformax = sample( which( x[,variable] == maxname ), samplethismany, replace=FALSE )
              otherinds = which( x[,variable] != maxname )
              x=x[c(otherinds,indsformax),]
              }
          }
      } else if ( method == 'over') {
          minmaxdifftbl = table( x[,variable] )
          mintbl = min( minmaxdifftbl )
          diffis = min( minmaxdifftbl )/max( minmaxdifftbl )
          while( diffis < minimum_ratio ) {
              minmaxdifftbl = table( x[,variable] )
              maxtbl = max( minmaxdifftbl )
              diffis = min( minmaxdifftbl )/max( minmaxdifftbl )
              if ( diffis < minimum_ratio ) {
                  minname = names( minmaxdifftbl[ minmaxdifftbl == min(minmaxdifftbl) ] )
                  samplethismany = round( maxtbl * minimum_ratio )+1
                  indsformax = sample( which( x[,variable] == minname ), samplethismany, replace=TRUE )
                  otherinds = which( x[,variable] != minname )
                  x=x[c(otherinds,indsformax),]
                  }
              }
      } else if ( method == 'mwmote' ) {
        minmaxdifftbl = table( x[,variable] )
        diffis = min( minmaxdifftbl )/max( minmaxdifftbl )
        while ( diffis < minimum_ratio ) {    
            minmaxdifftbl = table( x[,variable] )
            diffis = min( minmaxdifftbl )/max( minmaxdifftbl )
            if ( diffis < minimum_ratio ) {
                samplethismany = round( (max( minmaxdifftbl ) - min( minmaxdifftbl )) * minimum_ratio ) + 1
                temp=imbalance::mwmote( x, samplethismany, classAttr = variable )
                x=rbind( x, temp)
            }
        }
    }
  return( x )
}

#' balanced data partition
#' 
#' caret-based data train-test data partition function compatible with mlr3
#' 
#' @param x target vector to split to train test
#' @param perc percentage ( 0 to 1 ) to place in training partition
#' @param subjectIDs unique IDs per subject; aids with repeated measurement data
#' by ensuring subjects exist in uniquely in either train or test
#' 
#' @return vector of strings
#' @author Avants BB
#' @export
dataPartition <- function( x, perc, subjectIDs=NULL ) {
  if ( is.null( subjectIDs ) ) {
    train=caret::createDataPartition( x, p=perc )$Resample1
    test = c(1:length(x))
    test = test[ ! ( test %in% train) ]
    return( list( train=train, test=test ) )
  } else {
    if ( ! is.numeric( x ) ) {
      xuse = as.numeric( as.factor(x) )
    } else xuse = x
    xdf=data.frame( x=xuse, sid=subjectIDs )
    xdf=aggregate( x ~ sid, data=xdf, median )
    train=caret::createDataPartition( xdf$x, p=perc )$Resample1
    trainIDs = xdf$sid[train]
    testIDs = xdf$sid[-train]
    test = which( subjectIDs %in% testIDs )
    train = which( subjectIDs %in% trainIDs )
    return( list( train=train, test=test ) )
  }
}


#' mlr3classifiers
#' 
#' good classification learners from mlr3 (in my experience)
#' 
#' @param twoclass boolean
#' @param all boolean
#' 
#' @return vector of strings
#' @author Avants BB
#' @export
mlr3classifiers <- function( twoclass=TRUE, all=FALSE ) {
    if ( all ) {
      mymsg = 'Use \n mlr3extralearners::list_mlr3learners(select = c("id", "required_packages")) '
      message(mymsg)
      # myl = mlr3extralearners::list_mlr3learners(select = c("id", "required_packages"))
      # myl = myl[ -grep("RWeka",myl$required_packages), ]
      # myl2 = as.data.table(mlr_learners)
      # myl2 = myl2[ -grep("RWeka",myl2$packages), ]
      return( mymsg )
    }
    mylearners =  paste0( "classif.", 
        c("gbm","glmnet","kknn", 'ranger', 'ksvm',# 'mob',
        # "lssvm", 
        "rpart", "xgboost", #"e1071",
         "randomForest", "fnn",
        # 'liblinear', 
        'naive_bayes' ) )
    twoclassers=c('classif.imbalanced_rfsrc')
    if ( twoclass ) mylearners = c( mylearners, twoclassers )
    return( mylearners )
}

#' mlr3classifier
#' 
#' build a mlr3 classification model
#'
#' @param dfin dataframe input
#' @param tcols columns for the prediction task - first is the target outcome
#' @param learnerName which mlr3 learner to instantiate
#' @param partrate partition ratio for the training 0.8 equals 80 percent train 20 test
#' @param dup_size integer for over/under/smote sampling
#' @param balancing string over, under, smote, rwo, mwmote, racog, none are the options
#' @param subjectIDs unique IDs per subject; aids with repeated measurement data
#' by ensuring subjects exist in uniquely in either train or test
#' @param verbose boolean
#' @return dataframe with task, learner, accuracy and balanced accuracy
#' @author Avants BB
#' @export
mlr3classification <- function( dfin, tcols, learnerName, partrate=0.80, dup_size=0, balancing="smote", subjectIDs=NULL, verbose=TRUE ) {
    tarzan = tcols[1]
    valbal = c("none","over","under","smote", "rwo", "racog", "mwmote" )
    if ( ! ( balancing %in% valbal ) )
        stop(paste("balancing must be one of ", paste(valbal, collapse=" / ")))
    if ( sum( table( dfin[,tarzan] ) > 0 ) == 2 ) {
        twoclass=TRUE 
        dfin[,tarzan]=as.factor(as.character(dfin[,tarzan]))
        } else twoclass=FALSE
    selectviz = subtyper::fs( !is.na( dfin[,tarzan])  )
    trainppmisplit = dfin[ selectviz, tcols ]
    split = dataPartition(trainppmisplit[,tarzan], partrate, subjectIDs )
    task_penguins = as_task_classif( formula(paste(tarzan, " ~ .")), data = trainppmisplit)

    if ( balancing %in% c( "rwo", "racog", "mwmote" ) & dup_size > 1 ) {
      temptbl = table( trainppmisplit[,tarzan] )
      myinstsize = min( temptbl ) * ( dup_size - 1 )
      if ( balancing == 'mwmote')
        balDF = mwmote( dataset = trainppmisplit, numInstances = myinstsize, classAttr = tarzan)
      if ( balancing == 'rwo')
        balDF = rwo( dataset = trainppmisplit, numInstances = myinstsize, classAttr = tarzan)
      if ( balancing == 'racog')
        balDF = racog( dataset = trainppmisplit, numInstances = myinstsize, classAttr = tarzan)
      newrownum=(nrow(trainppmisplit)+1)
      trainppmisplit = rbind( trainppmisplit, balDF )
      split$train = c( split$train, newrownum:nrow(trainppmisplit) )
    }

    # SMOTE enriches the minority class with synthetic data
    if ( balancing == 'smote' )
        gr_smote =
            po("colapply", id = "int_to_num",
                applicator = as.numeric, 
                affect_columns = selector_type("integer")) %>>%
            po("smote", dup_size = dup_size) %>>%
            po("colapply", id = "num_to_int",
                applicator = function(x) as.integer(round(x, 0L)), 
                affect_columns = selector_type("numeric"))
            # enrich minority class by factor (dup_size + 1)

    # oversample majority class (relative to majority class)
    if ( balancing == 'over' )
        gr_smote = po("classbalancing",
            id = "oversample", adjust = "minor",
            reference = "minor", shuffle = FALSE, ratio = dup_size)
            # enrich minority class by factor 'ratio'

    if ( balancing == 'under' )
        gr_smote = po("classbalancing",
            id = "undersample", adjust = "major",
            reference = "major", shuffle = FALSE, ratio = 1 / dup_size )

    jj = learnerName
    learner = lrn( jj )
    if ( dup_size > 0 & jj != "classif.imbalanced_rfsrc" & 
      balancing %in% c('smote','over','under') ) {
        learner = as_learner(gr_smote %>>% learner)
      } 
    learner$train(task_penguins, split$train)
    prediction = learner$predict( task_penguins, split$test )
    measure = msr("classif.bacc")
    bacc=prediction$score(measure)
    measure = msr("classif.acc")
    acc=prediction$score(measure)
       
    return( list( task=task_penguins, split=split, learner=learner, acc=acc, bacc=bacc ) )
}


#' mlr3classifiercv
#' 
#' cross-validate a mlr3 classification model
#'
#' @param dfin dataframe input
#' @param tcols columns for the prediction task - first is the target outcome
#' @param nrepeats number of subsampling-driven train test runs
#' @param partrate partition ratio for the training 0.8 equals 80 percent train 20 test
#' @param dup_size integer for over/under/smote sampling
#' @param balancing string over, under, smote, none are the options
#' @param mylearners the mlr3 learners over which to search ; defaults to those that support multiclass
#' @param subjectIDs unique IDs per subject; aids with repeated measurement data
#' by ensuring subjects exist in uniquely in either train or test
#' @param verbose boolean
#' @return dataframe quantifying performance
#' @author Avants BB
#' @export
mlr3classifiercv <- function( dfin, tcols, nrepeats=10, partrate=0.80, dup_size=0, balancing="smote", mylearners = mlr3classifiers(), subjectIDs=NULL, verbose=TRUE ) {
    tarzan = tcols[1]
    valbal = c("none","over","under","smote", "rwo", "racog", "mwmote" )
    if ( ! ( balancing %in% valbal ) )
        stop(paste("balancing must be one of ", paste(valbal, collapse=" / ")))
    if ( sum( table( dfin[,tarzan] ) > 0 ) == 2 ) {
        twoclass=TRUE 
        dfin[,tarzan]=as.factor(as.character(dfin[,tarzan]))
        } else twoclass=FALSE
    selectviz = subtyper::fs( !is.na( dfin[,tarzan])  )
    trainppmi = dfin[ selectviz, tcols ]
    print( table( trainppmi[,tarzan] ))
    clsdf = data.frame()
    for ( r in 1:nrepeats ) {
        if ( verbose ) print(paste("repeat",r))
        split = dataPartition(trainppmi[,tarzan], partrate, subjectIDs )
        trainppmisplit = trainppmi
        if ( balancing %in% c( "rwo", "racog", "mwmote" )  ) {
          temptbl = table( trainppmisplit[,tarzan] )
          myinstsize = max( temptbl ) - min( temptbl )
          if ( balancing == 'mwmote')
            balDF = mwmote( dataset = trainppmisplit[split$train,], numInstances = myinstsize, classAttr = tarzan)
          if ( balancing == 'rwo')
            balDF = rwo( dataset = trainppmisplit[split$train,], numInstances = myinstsize, classAttr = tarzan)
          if ( balancing == 'racog')
            balDF = racog( dataset = trainppmisplit[split$train,], numInstances = myinstsize, classAttr = tarzan)
          newrownum=(nrow(trainppmisplit)+1)
          print( paste("adding",nrow(balDF),"extrapolated data via", balancing) )
          trainppmisplit = rbind( trainppmisplit, balDF )
          split$train = c( split$train, newrownum:nrow(trainppmisplit) )
        }

        task_penguins = as_task_classif( formula(paste(tarzan, " ~ .")), data = trainppmisplit)

        # SMOTE enriches the minority class with synthetic data
        if ( balancing == 'smote' )
        gr_smote =
            po("colapply", id = "int_to_num",
                applicator = as.numeric, 
                affect_columns = selector_type("integer")) %>>%
            po("smote", dup_size = dup_size) %>>%
            po("colapply", id = "num_to_int",
                applicator = function(x) as.integer(round(x, 0L)), 
                affect_columns = selector_type("numeric"))
            # enrich minority class by factor (dup_size + 1)

        # oversample majority class (relative to majority class)
        if ( balancing == 'over' )
        gr_smote = po("classbalancing",
            id = "oversample", adjust = "minor",
            reference = "minor", shuffle = FALSE, ratio = dup_size)
            # enrich minority class by factor 'ratio'

        if ( balancing == 'under' )
        gr_smote = po("classbalancing",
            id = "undersample", adjust = "major",
            reference = "major", shuffle = FALSE, ratio = 1 / dup_size )

        for ( jj in mylearners ) {
            if ( verbose ) print(paste("mylearners",jj))
            learner = lrn( jj )
            if ( dup_size > 0 & jj != "classif.imbalanced_rfsrc" & 
                balancing %in% c('smote','over','under') )
                learner = as_learner(gr_smote %>>% learner)
            learner$train(task_penguins, split$train)
            prediction = learner$predict( task_penguins, split$test )
            if ( verbose )
              print( prediction$confusion )
            measure = msr("classif.bacc")
            bacc=prediction$score(measure)
            measure = msr("classif.acc")
            acc=prediction$score(measure)
            n=nrow(clsdf)+1
            clsdf[n,'mdl']=jj
            clsdf[n,'bacc']=bacc
            clsdf[n,'acc']=acc
            clsdf[n,'r']=r
            }
        }
       
    return( clsdf )
}


#' dcurvarsel
#' 
#' variable selection with CUR - using default parameters from package example
#'
#' @param curdf dataframe input
#' @param variablenames columns for the task
#' @param fraction scalar between zero and one
#' @return vector of selected variables
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' vtosel = getNamesFromDataframe("Random",mydf)
#' # vimp = dcurvarsel( mydf, vtosel, 0.5 )
#' @export
dcurvarsel <- function( curdf, variablenames, fraction ) {
  curEnv=environment()
  myk = min(c(20, dim(curdf[,variablenames])))-1
  loccur = dCUR::CUR
  environment(loccur) <- curEnv
  pdf(file="/dev/null")
  result = loccur( 
              data=curdf, 
              variables=variablenames,
              k=myk, rows = 1, columns = .2, standardize = TRUE,
              cur_method = "mixture" )
  temp=capture.output(dev.off())
  return( head( result$leverage_columns_sorted$var_names, 
    round(length(variablenames)*fraction) ) )
  }


#' clearcolname
#' 
#' remove a column from a dataframe
#'
#' @param mydf dataframe input
#' @param mycolname columns for the task
#' @return trimmed dataframe
#' @author Avants BB
#' @export
clearcolname = function( mydf, mycolname ) {
    if ( mycolname %in% colnames( mydf ) ) {
        mydf=mydf[,-grep(mycolname,colnames( mydf))]
    }
    return( mydf )
}



#' consensusSubtypingTrain
#' 
#' apply several clustering methods and append results to a dataframe;
#' usually will include one subject per row
#'
#' @param dataToClust dataframe input that contains the relevant variables (may have others as well)
#' @param featureNames names to use in the clustering
#' @param clustVec names of the clustering methods to use
#' @param ktrain the number of clusters
#' @param reorderingVariable the name of the column to use to reorder the cluster names
#' @param mvcl character prefix for the new cluster column names
#' @param ksearch the cluster number(s) for clustering of the methods by concordance
#' @param verbose boolean
#' @return a list with newdata: dataframe with new variables attached; models contains the trained models; reorderers contains a dataframe for reordering cluster levels; quality measurements and clustering based on cluster concordance (adjusted rand index) are also returned.
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @export
consensusSubtypingTrain = function( dataToClust, featureNames, clustVec, ktrain, reorderingVariable, mvcl='MVST', ksearch=3, verbose=FALSE ) {

    # from mclust - avoid import for this 1 function
    adjustedrandindex = function( x, y ) { 
        x <- as.vector(x)
        y <- as.vector(y)
        if (length(x) != length(y)) 
            stop("arguments must be vectors of the same length")
        tab <- table(x, y)
        if (all(dim(tab) == c(1, 1))) 
            return(1)
        a <- sum(choose(tab, 2))
        b <- sum(choose(rowSums(tab), 2)) - a
        c <- sum(choose(colSums(tab), 2)) - a
        d <- choose(sum(tab), 2) - a - b - c
        ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b + 
            a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
        return(ARI)
      }

    if ( missing( reorderingVariable ) ) reorderingVariable = featureNames[1]
    stopifnot( reorderingVariable %in% colnames(dataToClust) )
    stopifnot( all( featureNames %in% colnames(dataToClust) ) )
    clustmodels = list()
    reoModels = list()
    newclustnames = paste0(mvcl,"_",clustVec)
    successfulclust = rep( FALSE, length(clustVec) )
    names(successfulclust)=clustVec
    for ( myclust in clustVec ) {
        if ( verbose ) print(paste("Train:",mvcl,myclust))
        mvclLocal = paste0(mvcl,"_",myclust)
        clustmodels[[ myclust ]] =
            trainSubtypeClusterMulti( dataToClust, 
                featureNames, myclust, desiredk=ktrain )
        turkeyshoot = predictSubtypeClusterMulti( dataToClust, 
            featureNames, clustmodels[[ myclust ]], mvclLocal )
        if ( length(table(turkeyshoot[,mvclLocal])) == ktrain ) {
            reodf = reorderingDataframe( turkeyshoot, mvclLocal, reorderingVariable )
            reodfCheck = reodf
            ct = 0
            if ( verbose )
              print(paste("st:",mvclLocal))
            while ( ! all( reodfCheck$x == reodfCheck$ord_x ) & ct < 3 ) {
                temp = predictSubtypeClusterMulti( dataToClust, 
                    featureNames, clustmodels[[ myclust ]], mvclLocal,reorderingDataframe=reodf )
                reodfCheck = reorderingDataframe( temp, mvclLocal, reorderingVariable )
                if (  ! all( reodfCheck$x == reodfCheck$ord_x )  ) reodf=reodfCheck
                ct = ct + 1
                }
            reodf$method=myclust
            reoModels[[myclust]]=reodf
            dataToClust = clearcolname(dataToClust, mvclLocal )
            dataToClust = predictSubtypeClusterMulti( dataToClust, 
                    featureNames, clustmodels[[ myclust ]], mvclLocal,reorderingDataframe=reodf )
            successfulclust[myclust]=TRUE
            }
        }


  clustVecGood=names(successfulclust[successfulclust])
  nn = length(clustVec[successfulclust])
  mystat = matrix(nrow=nn,ncol=nn)
  colnames(mystat)=rownames(mystat)=clustVecGood
  for ( j in newclustnames ) {
    if ( j %in% colnames(dataToClust) )
    for ( k in newclustnames ) {
        if ( k != j & k %in% colnames(dataToClust)  ) {
            cl1 = dataToClust[,j]
            cl2 = dataToClust[,k]
            myari = adjustedrandindex(cl1,cl2)
            # print(paste("j",j,"k",k,"ARI",myari))
            mystat[gsub(paste0(mvcl,"_"),"",j),gsub(paste0(mvcl,"_"),"",k)]=myari
#            else if ( agreestat == 'conf') {
#                cm <- caret::confusionMatrix(factor(cl1), reference = factor(cl2))
#                mystat[j,k]=cm$overall[1]
#            } else if ( agreestat == 'kap' ) mystat[j,k]=psych::cohen.kappa(cbind(cl1,cl2))$kappa
#            mystat[j,k]=chisq.test(cl1,cl2,simulate.p.value = TRUE)$statistic
#            print(paste("chi statistic=",chisq.test(cl1,cl2)$statistic ))
            }
          }
        }
    mymeans = apply( mystat, FUN=mean, MARGIN=2, na.rm=T)
    mymaxes = apply( mystat, FUN=max, MARGIN=2, na.rm=T)
    mostdisagreeable = names(mymeans[fs(mymeans==min(mymeans,na.rm=T))])
    mostagreeable = names(mymeans[fs(mymeans==max(mymeans,na.rm=T))])
    if ( verbose ) {
      message(paste("The most agreeable method is:",mostagreeable))
      print(paste("The most agreeable method is:",mostagreeable))
      message(paste("The most disagreeable method is:",mostdisagreeable))
      print(paste("The most disagreeable method is:",mostdisagreeable))
    }
    mystat2=mystat
    diag(mystat2)=NA
    if ( verbose ) pheatmap::pheatmap( mystat2 )
    if ( any( is.na( ksearch ) ) ) ksearch = 2:round(nn/2)
    mypam = fpc::pamk(mystat2,krange=ksearch,criterion="asw", usepam=TRUE,
        scaling=FALSE, alpha=0.001, critout=FALSE, ns=10 )
    myClusterCluster = mypam$pamobject$clustering
    myClusterScores = rowMeans(mypam$pamobject$medoids,na.rm=T)
    clustind = which( myClusterScores == 
        sort(myClusterScores)[length(myClusterScores) ] )# 2nd best
    moreagreeable = names( myClusterCluster[ myClusterCluster == clustind ] )
    if ( verbose ) {
      message("ClusterCluster")
      print("ClusterCluster")
      print(myClusterCluster)
      print("ClusterSimilarityScores")
      print( myClusterScores )
      message("highestARIMethods")
      print("highestARIMethods")
      print(moreagreeable)
      } 

    return( list( 
            newdata = dataToClust,
            models = clustmodels,
            reorderers = reoModels,
            highestARIMethods = moreagreeable,
            ClusterClusters = myClusterCluster,
            ClusterSimilarityScores = myClusterScores
            )
        )

    }

#' consensusSubtypingPredict
#' 
#' apply several clustering methods and append results to a dataframe;
#' assumes trained models already exist.
#'
#' @param dataToClust dataframe input that contains the relevant variables (may have others as well)
#' @param featureNames names to use in the clustering
#' @param clustVec names of the clustering methods to use
#' @param clustmodels the trained models
#' @param reorderers the reordering data frames
#' @param mvcl character prefix for the new cluster column names
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @return new dataframe with new variables attached
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @export
consensusSubtypingPredict = function( dataToClust, featureNames, clustVec, clustmodels, 
  reorderers, mvcl, idvar, visitName, baselineVisit  ) {
    stopifnot( all(featureNames %in% colnames(dataToClust)))
    namestoclear = getNamesFromDataframe(mvcl,dataToClust)
    for ( nm in namestoclear )
        dataToClust = clearcolname(dataToClust, nm )
    for ( myclust in clustVec ) {
        mvclLocal = paste0(mvcl,"_",myclust)
        if ( myclust %in% names(reorderers) & myclust %in% names(clustmodels) ) {
          if ( ! missing(idvar) & ! missing(visitName) & ! missing( baselineVisit ) ) {
            stopifnot( all(c(idvar,visitName) %in% colnames(dataToClust)))
            stopifnot( any( baselineVisit %in% dataToClust[,visitName] ) )
            dataToClust = predictSubtypeClusterMulti( dataToClust, 
                    featureNames, clustmodels[[ myclust ]], mvclLocal, 
                    idvar, visitName, baselineVisit, 
                    reorderingDataframe=reorderers[[myclust]] )
          } else {
            dataToClust = predictSubtypeClusterMulti( dataToClust, 
                    featureNames, clustmodels[[ myclust ]], mvclLocal, 
                    reorderingDataframe=reorderers[[myclust]] )
          }
        }
      }
    return( dataToClust )
    }


#' consensusSubtypingCOCA
#' 
#' apply consensus clustering give several clustering solutions
#'
#' @param dataToClust dataframe input that contains the relevant variables (may have others as well)
#' @param targetk the desired number of classes
#' @param cocanames names of columns to use for the consensus
#' @param newclustername the column name for the consensus clustering
#' @param reorderingVariable the name of the column to use to reorder the cluster names
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @param consensusmethod either kmeans or hclust
#' @param verbose boolean
#' @return new dataframe with new variables attached
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @importFrom caret dummyVars contr.ltfr
#' @export
consensusSubtypingCOCA = function( dataToClust, targetk, cocanames, newclustername, reorderingVariable, idvar, visitName, baselineVisit, consensusmethod='kmeans',verbose=TRUE ) {
    # assume we already ran consensuscluster
    if ( !missing(idvar) )
      stopifnot( idvar %in% colnames(dataToClust) )
    if ( !missing(visitName) )
      stopifnot( visitName %in% colnames(dataToClust) )
    if ( !missing(reorderingVariable) )
      stopifnot( reorderingVariable %in% colnames(dataToClust) )
    if ( !missing(baselineVisit) )
      stopifnot( baselineVisit %in% dataToClust[,visitName] )
    stopifnot( all( cocanames %in% colnames(dataToClust) ) )
    usebaseline=FALSE
    if ( !missing(idvar) & !missing(visitName) & !missing(baselineVisit) )
      usebaseline=TRUE
    if ( length(cocanames) == 1 ) {
      dataToClust[,newclustername]=dataToClust[,cocanames]
      message("COCA irrelevant - only a single clustering result is available.")
      return( dataToClust )
    }
    dataToClust = clearcolname(dataToClust, newclustername )
    isbl = rep(TRUE,nrow(dataToClust))
    if ( usebaseline )
      isbl = dataToClust[,visitName] == baselineVisit
    cocoform = paste("~", paste( cocanames, collapse="+" ))
    if ( verbose ) {
        for ( x in cocanames ) {
            print(table(dataToClust[isbl,x] ))
        }
    }
    dmy = dummyVars(cocoform, data = dataToClust[,cocanames])
    dmytx = data.frame(predict(dmy, newdata = dataToClust[isbl,cocanames]))
    cocatx = coca::coca(dmytx, K = targetk, B=1000, maxIterKM=5000, ccClMethod=consensusmethod )
#    cocatx = coca::coca(dmytx, maxK = 6, B=5000 )
#    coca = coca::coca( dmytx, maxK = 10, hclustMethod = "average")
    cocatxlab = as.numeric( cocatx$clusterLabels )
    if ( verbose )
      print( table( cocatxlab ) )
    dataToClust[,newclustername]=NA
    dataToClust[isbl,newclustername]=cocatxlab
    if ( usebaseline ) {
        temp = fillBaselineColumn( dataToClust,
            newclustername, 
            idvar, visitName, baselineVisit, 
            fast=TRUE, verbose=FALSE )[[1]]
        dataToClust[rownames(temp),newclustername]=temp[,paste0(newclustername,'_BL')]
    }
    if ( !missing(reorderingVariable) ) {
      # now reorder 
      xdf=aggregate( dataToClust[isbl,reorderingVariable], 
        list(dataToClust[isbl,newclustername]),  
        mean, na.rm=TRUE )
      newreo=order(xdf[,2])
      olabels = dataToClust[,newclustername]
      placeholder = olabels
      for ( jj in 1:nrow(xdf) ) placeholder[ olabels == newreo[jj] ] = xdf[jj,1]
      dataToClust[,newclustername] = placeholder
      xdf=aggregate( dataToClust[isbl,reorderingVariable], list(dataToClust[isbl,newclustername]), mean, na.rm=TRUE )
      if ( verbose ) print(xdf)
      }
    dataToClust[,newclustername] = paste0(newclustername, dataToClust[,newclustername] )
    return( dataToClust )
    }
