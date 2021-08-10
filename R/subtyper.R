# Here's a good place to put your top-level package documentation

.onLoad <- function (lib, pkgname="subtyper") {
    ## Put stuff here you want to run when your package is loaded
    invisible()
}





#' Generate example data for illustrating `subtyper`
#'
#' @param n numer of subjects
#'
#' @return data frame
#' @author Avants BB
#' @examples
#'
#' mydf = generateSubtyperData( 5 )
#'
#' @export
generateSubtyperData <-function( n = 100 ) {

  myids = c("A","B","C","D","E", "F", "G", "H", "I", "J", "K", "L" )
  mygroup = c( rep("G0",4), rep("G1",4), rep("G2",4) )
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
  mymeans = c( 1, 3, 8, 1, 4, 10, 1, 5, 16 )
  ct = 1
  for ( v in vizzes ) {
    for ( k in unique( mygroup ) ) {
      losel = mydf$DX == k & mydf$visit == v
      cognition[ losel ] = rnorm( sum(losel), mean = mymeans[ct] )
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
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' summ = plotSubtypeChange( mydf, "Id", "cognition", "DX", "visit" )
#' @export
#' @importFrom stats lm predict qt rnorm var
#' @importFrom DDoutlier  LOOP  LOF  INFLO  RDOS  KDEOS  LDF  KNN_AGG  KNN_IN  KNN_SUM  RKOF
#' @importFrom ggplot2 aes ylim guides theme_bw scale_colour_hue geom_errorbar position_dodge element_text geom_line geom_point ggplot guide_legend
#' @importFrom ggplot2 xlab ylab theme rel
#' @importFrom plyr ddply rename
plotSubtypeChange <-function( mxdfin,
                           idvar,
                           measurement,
                           subtype,
                           vizname='timer', # should be roundable to integer
                           whiskervar = 'ci',
                           extra = "",
                           xlab = '',
                           ylab = ''
                          )
    {

      if ( ! ( whiskervar %in% c("ci","se") ) )
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
          datac <- ddply(data, groupvars, .drop=.drop,
            .fun = function(xx, col) {
              c(N    = length2(xx[[col]], na.rm=TRUE),
                mean = mean   (xx[[col]], na.rm=TRUE),
                sd   = mysd     (xx[[col]], na.rm=TRUE)
              )
            },
            measurevar
          )
          # Rename the "mean" column
          datac <- rename(datac, c("mean" = measurevar))

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
    tgcWithin$ymin = tgcWithin[,"Changer"] - tgcWithin[,whiskervar]
    tgcWithin$ymax = tgcWithin[,"Changer"] + tgcWithin[,whiskervar]
    yran = range( c( tgcWithin$ymin,tgcWithin$ymax ) )
    pd <- position_dodge( 0.05 ) # move them .05 to the left and right
    ww = 0.5
    namer=paste( subtype, ' vs ', measurement , extra )
    print( ggplot( tgcWithin,
      aes(x=timer, y=Changer, group = Quant, color = Quant , shape = Quant )) +
      geom_errorbar(aes(ymin=ymin, ymax=ymax), colour="black", width=ww, position=pd) +
      geom_line(lwd = 1, show.legend = FALSE ) +
      geom_point(position=pd, size=3, shape=21, fill="white") + # 21 is filled circle
      xlab( xlab ) + ylab( paste(ylab, measurement," +/-", whiskervar )) +
      ylim(yran[1], yran[2]) +
      theme(strip.text = element_text(size = rel(1.1))) +
      theme(legend.position = 'top') +
      guides(fill=guide_legend(title.position="top")) +
      scale_colour_hue(name=subtype, l=40) + # Use darker colors, lightness=40
#        ggtitle( namer ) + scale_y_continuous(breaks=0:20*4) +  # Set tick every 4
        theme_bw() +
        theme(text = element_text(size=14),
            axis.text.x = element_text(angle=0, vjust=1),
            legend.direction = "horizontal", legend.position = "top" )
      )
    return( tgcWithin )
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
#' and (optionally) longitudinal data.  The quality measure should be a "higher"
#' is better" measurement.
#'
#' @param mxdfin Input data frame with repeated measurements and a grouped time variable
#' @param idvar variable name for unique subject identifier column
#' @param visitvar variable name for the visit column
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
    for ( v in vizzes ) {
      losel = mxdfin[,idvar] == u & mxdfin[,visitvar] == v
      mysnr = mxdfin[losel,qualityvar]
      myw = which( losel )
      if ( any( !is.na(mysnr) )  )
        useit[ myw[ which.max(mysnr) ] ] = TRUE
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
#' @param report boolean return vector string of column names for outlier detection results
#' @return data frame
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' rbfnames = names(mydf)[grep("Random",names(mydf))]
#' mydfol = outlierness( mydf, rbfnames )
#' @export
outlierness <- function( mxdfin, measureColumns, calck, report=FALSE ) {
  olnames = c("OL_LOOP", "OL_LOF", "OL_NOF", "OL_INFLO", "OL_RDOS", "OL_KDEOS",
    "OL_LDF", "OL_KNN_AGG", "OL_KNN_IN", "OL_KNN_SUM", "OL_RKOF" )
  olnames = c("OL_LOOP", "OL_LOF",  "OL_INFLO", "OL_RDOS", "OL_KDEOS",
    "OL_LDF", "OL_KNN_AGG", "OL_KNN_IN", "OL_KNN_SUM", "OL_RKOF" )
  if ( report ) return( olnames )
  oldata = scale( mxdfin[,measureColumns], TRUE, TRUE )
  if ( missing( calck ) )
    calck = DDoutlier::NAN( oldata )$r # this + LOOP is good
  mxdfin$OL_LOOP = DDoutlier::LOOP(  oldata, k=calck )
  mxdfin$OL_LOF = DDoutlier::LOF(  oldata, k=calck )
  mxdfin$OL_INFLO = ( DDoutlier::INFLO( oldata, k=calck ) )
#  mxdfin$OL_NOF = capture.output( DDoutlier::NOF( oldata )$NOF)
  mxdfin$OL_KDEOS = DDoutlier::KDEOS( oldata, k_min = round( calck * 0.5 ), k_max = round( calck * 2 ) )
  mxdfin$OL_LDF = DDoutlier::LDF( oldata, k = calck  )$LDF
  mxdfin$OL_KNN_AGG = DDoutlier::KNN_AGG( oldata, k_min = round( calck * 0.5 ), k_max = round( calck * 2 ) )
  mxdfin$OL_KNN_IN = -1.0 * DDoutlier::KNN_IN( oldata, k = calck )
  mxdfin$OL_KNN_SUM = DDoutlier::KNN_SUM( oldata, k = calck  )
  mxdfin$OL_RDOS = DDoutlier::RDOS( oldata, k = calck  )
  mxdfin$OL_RKOF = DDoutlier::RKOF( oldata, k = calck  )
  return( mxdfin )
}

# filterByQuality  <- function( ) {}  # strategies to identify a high quality subset of data



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
    gsel = mxdfin[,groupVariable] == group
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
    gsel = mxdfin[,groupVariable] == group
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
  baselineVisit ) {

  msr = as.character( subtypeDataFrame[, "measurement"][1] )
  thesubtypes = subtypeDataFrame[,1]
  defaultST = tail( thesubtypes, 1 )
  mxdfin[,names(subtypeDataFrame)[1]] = defaultST
  quantsV = as.numeric( subtypeDataFrame[,"quantileValues"][-1] )
  quantsP = as.numeric( subtypeDataFrame[,"quantiles"][-1] )
  # by default, just run across each row
  mxdfin[,names(subtypeDataFrame)[1]] = Hmisc::cut2( mxdfin[,msr], cuts = quantsV )
  mxdfin[,names(subtypeDataFrame)[1]] = as.character( mxdfin[,names(subtypeDataFrame)[1]] )
  theselevs = unique(  mxdfin[,names(subtypeDataFrame)[1]] )
  for ( k in 1:length( theselevs ) ) {
    losel = mxdfin[,names(subtypeDataFrame)[1]] == theselevs[k]
    mxdfin[,names(subtypeDataFrame)[1]][ losel ] = thesubtypes[k]
  }
  mxdfin[,names(subtypeDataFrame)[1]] <- factor(mxdfin[,names(subtypeDataFrame)[1]], levels=thesubtypes )
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
