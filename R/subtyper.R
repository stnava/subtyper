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



#' Convert left/right variables to a measure of asymmetry
#'
#' @param mydataframe dataframe containing relevant variables
#' @param leftvar left variable names
#' @param replacer string to replace left with in column names of output
#' @return fixed x
#' @author Avants BB
#' @export
mapAsymVar <- function( mydataframe, leftvar, replacer='Asym' ) {
  rightvar =  gsub( "left", "right", leftvar )
  hasright = rightvar %in% colnames(mydataframe)
  temp = mydataframe[,leftvar[hasright]] - mydataframe[,rightvar[hasright]]
  temp = temp * sign(temp )
  newnames = gsub("left", replacer,leftvar[hasright])
  colnames(temp)=newnames
  return( temp )
}



#' Convert left/right variables to an average measurement
#'
#' @param mydataframe dataframe containing relevant variables
#' @param leftvar left variable names
#' @param replacer string to replace left with in column names of output
#' @return fixed x
#' @author Avants BB
#' @export
mapLRAverageVar <- function( mydataframe, leftvar, replacer='LRAVG' ) {
  rightvar =  gsub( "left", "right", leftvar )
  hasright = rightvar %in% colnames(mydataframe)
  temp = mydataframe[,leftvar[hasright]] * 0.5 + mydataframe[,rightvar[hasright]] * 0.5
  newnames = gsub("left", replacer,leftvar[hasright])
  colnames(temp)=newnames
  return( temp )
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
#' mydfhq = averageRepeats( mydf, "Id", "visit", "quality")
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
   deltaExt = "_delta"
) {
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
  mxdfin[,newcolname]=NA
  mxdfin[,newcolnamed]=NA
  visitidisnumeric = class(mxdfin[,visitID]) == "numeric"
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
      }
      if ( sum( selbase ) > 0  & isFactor ) {
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
  baseform = paste( outcomevar, " ~ 1" )
  imodel = lm( baseform, data=subdf ) # just the intercept
  ctlmodel = lm( adjustmentFormula, data=subdf )
  predvol = predict( ctlmodel, newdata = mxdfin ) - predict( imodel, newdata = mxdfin  )
  adjustedoutcome = paste0( outcomevar, "_adjusted" )
  mxdfin[ , adjustedoutcome ] = mxdfin[,outcomevar] - predvol
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

  ktypes = c( "kmeans", 'kmeansflex', "GMM", "mclust", "pamCluster", 
    "kmeansflex","kmedians",  "angle",  "ejaccard", "flexcorr",
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
  subdf = data.matrix( subdf )

  if ( ismlr3 ) {
    task = as_task_clust(subdf)
    learner = mlr_learners$get(method)
    learner$param_set$values = list(centers = 3L)
    myp = learner$train(task)
    predict(myp,newdata=subdf)
    return( myp )
  }

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
    "hardcl","neuralgas", "flexcorr")
  if ( method %in% flexmeth ) {
    initk = ClusterR::KMeans_rcpp(subdf,
        clusters = desiredk, num_init = 5, max_iters = 100,
        initializer = 'optimal_init', verbose = F, fuzzy=TRUE)$clusters
    if ( method == 'flexcorr') {
      ejacFam <- flexclust::kccaFamily(dist=distCor,cent=centMean)
      mycl = flexclust::kcca(
          subdf,
          k = initk,
          weights=flexweights, group=flexgroup,
        family = ejacFam )
      return( mycl )
    }
    if ( method %in% c("hardcl","neuralgas"))
      return( 
        cclust(subdf, k=initk, weights=flexweights, group=flexgroup ) )
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
  distance_metric = 'pearson_correlation'
) {
  myclustclass = class(clusteringObject)
  if ( length(myclustclass) == 1 )
    myclustclass=c(myclustclass,'other')
  subdf = mxdfin[ , measureColumns ]
  subdf = data.matrix( subdf )
  if ( myclustclass[2] == "Gaussian Mixture Models" ) {
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
    cluster_labels = mypr$classification
    mxdfin = cbind( mxdfin, factor( paste0(clustername,cluster_labels) ) )
    colnames( mxdfin )[ ncol( mxdfin ) ] = clustername
  } else stop("Unknown class of clustering object.")
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
#' @importFrom fastICA fastICA
#' @importFrom mclust Mclust predict.Mclust mclustBIC
#' @importFrom fpc pamk
#' @importFrom flexclust kccaFamily kcca bootFlexclust cclust
#' @importFrom Evacluster pamCluster nmfCluster kmeansCluster hierarchicalCluster FuzzyCluster EMCluster 
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
#' some cases, depending on subtypes/degrees of freedom, only one will work.
#' @param covariates optional string of covariates
#' @param transform optional effect_size
#' @param significance_level to threshold effects
#' @param visualize boolean
#' @return dataframes for visualization that show feature to subtype importance e.g. via \code{pheatmap}
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' fimp = featureImportanceForSubtypes( mydf$DX, mydf[,rbfnames] )
#' fimp = featureImportanceForSubtypes( mydf$DX, mydf[,rbfnames], "subtypes2features" )
#' @importFrom fastICA fastICA
#' @importFrom effectsize t_to_d z_to_d
#' @export
featureImportanceForSubtypes <- function(
    dataframein,
    subtypeLabels,   # cluster
    featureMatrix,   # featureColumns
    associationType = c( "features2subtypes", "subtypes2features" ),
    covariates = "1",
    transform = 'effect_sizes',
    significance_level = 0.001,
    visualize = FALSE ) {

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
  if ( class(subtypeLabels) != 'factor' )
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
  clustmat = matrix( 0, nrow = length( subtypeLabels ), ncol=mync )
  colnames(clustmat) = as.character( uniqClusts )
  for ( j in 1:mync ) {
    losel = subtypeLabels == uniqClusts[j]
    if ( sum(losel) > 0 )
      clustmat[ losel , j] = 1
    }
  if ( associationType[1] == "features2subtypes" ) {
    # converts cluster labels to one-hot coding
    for ( j in 1:mync ) {
      c1reg = glm( factor( clustmat[,j]) ~ data.matrix(featureMatrix), family='binomial' )
      mycoffs = coefficients(summary(c1reg))
      sigthresh = rep( 0, ncol(featureMatrix))
      sigthresh[ mycoffs[-1,"Pr(>|z|)"] <= significance_level ] = 1
      myz = mycoffs[-1,"z value"]
      if ( transform == 'effect_sizes' ) 
        myz = as.numeric( effectsize::z_to_d( myz, nrow(featureMatrix) ) )
      clustzdescribe[j,]=myz
      clustsigdescribe[j,]=myz * sigthresh
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
#' @return vector of feature names
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' fimp = regressionBasedFeatureSelection( mydf$DX, mydf[,rbfnames] )
#' @export
regressionBasedFeatureSelection <- function( 
    dataframein, subtypeLabels, featureMatrix, covariates="1", n_features=25 ) {
    fimp = featureImportanceForSubtypes(
        dataframein,
        subtypeLabels, 
        featureMatrix, 
        "subtypes2features", 
        covariates=covariates )
    thefeats = c()
    myfimp = fimp$subtypeFeatureTScoresSignificant
    for ( k in 1:nrow( fimp$subtypeFeatureTScores )) {
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
#' @importFrom gaston as.matrix
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
    library(gaston)
    gwas=read.bed.matrix( rootFileName )
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
    mydf=data.frame(sid=rownames(y))
    mydf=cbind(mydf,y)
    return(mydf)
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
  glist[[length(glist)+1]]= ggscatter(indf, x = xvar, y = yvar, color=colorvar,   size=3.45,  point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ggtitle(paste(anatshow[1])) #+ theme(legend.position = "none")

  if ( is.numeric(indf[,anat]) ) {
    medsplit = median( indf[,anat], na.rm=T )
    hisel = indf[,anat] > medsplit
    loclev=loclev2=anat
  } else {
    loclev = (unique(indf[,anat])[1])
    hisel=indf[,anat]==loclev
    loclev2=paste0("!",loclev)
  }
  glist[[length(glist)+1]]=ggscatter(indf[hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ggtitle(paste('+',anatshow[2],loclev)) + theme(legend.position = "none")

  
  glist[[length(glist)+1]]=ggscatter(indf[!hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE  ) + theme(text = element_text(size=12))+ ggtitle(paste('-',anatshow[3],loclev2))+ theme(legend.position = "none")

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