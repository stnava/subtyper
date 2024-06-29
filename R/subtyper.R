# Here's a good place to put your top-level package documentation

.onLoad <- function (lib, pkgname="subtyper") {
    ## Put stuff here you want to run when your package is loaded
    invisible()
}


#' Filter Columns by Percentage of NA Values
#'
#' @description
#' This function filters out columns from a dataframe that have a percentage of
#' NA (missing) values greater or equal to a user-defined threshold. The percentage
#' is calculated with respect to the total number of rows in the dataframe.
#'
#' @param df A dataframe from which columns will be filtered.
#' @param percentThreshold A numeric value between 0 and 100 indicating the percentage
#'        threshold of NA values. Columns with NA percentages greater or equal to this
#'        value will be removed from the dataframe.
#'
#' @return A dataframe with columns filtered based on the NA percentage threshold.
#'
#' @examples
#' data(mtcars)
#' mtcars[1:5, 1:2] <- NA # Introduce NA values for demonstration
#' filtered_df <- filterNAcolumns(mtcars, 10)
#' print(filtered_df)
#'
#' @export
#'
#' @import dplyr
filterNAcolumns <- function(df, percentThreshold) {
  # Validate input
  if (!is.data.frame(df)) {
    stop("The input 'df' must be a dataframe.")
  }
  
  if (!is.numeric(percentThreshold) || percentThreshold < 0 || percentThreshold > 100) {
    stop("The 'percentThreshold' must be a numeric value between 0 and 100.")
  }
  
  # Calculate the percentage of NA values for each column
  na_percent <- sapply(df, function(col) {
    sum(is.na(col)) / nrow(df) * 100
  })
  
  # Filter columns based on threshold
  filtered_df <- df[, na_percent < (100.0-percentThreshold)]
  
  return(filtered_df)
}


#' Generate Predictors from ANTsPyMM Imaging Data
#'
#' This function generates a list of variable names to be used as predictors
#' in a model, based antspymm tabular version of imaging data.
#' It filters and processes the data to create meaningful predictors. LRAVG
#'
#' @param demog A dataframe containing demographic and imaging data.
#' @param doasym boolean
#' @param return_colnames boolean
#' @return A dataframe with processed demographic and imaging data.
#' @examples
#' # predictors <- antspymm_predictors(demographic_data)
#' @export
#'
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
antspymm_predictors <- function( demog, doasym=FALSE, return_colnames=FALSE ) {
  badcaud=getNamesFromDataframe("bn_str_ca",demog)
  badcaud=badcaud[ -grep("deep",badcaud)]
  xcl=c("hier_id",'background','SNR','evr','mask','msk','smoothing','minutes', "RandBasis",'templateL1', 'upsampl', 'paramset', 'nc_wm', 'nc_csf', 'censor','bandpass', 'outlier', 'meanBold', 'dimensionality', 'spc', 'org','andwidth',
  'unclassified', 'cleanup', 'slice', badcaud, 'dimx', 'dimy', 'dimz','dimt', 'modality' )
  if ( doasym & return_colnames ) xcl=c(xcl,'left','right',"_l_","_r_")
  t1namesbst = getNamesFromDataframe( c("T1Hier",'brainstem','vol'), demog, exclusions=c("tissues","lobes"))[-1]
  testnames=c(
          getNamesFromDataframe( "T1w_" , demog, exclusions=xcl),
          getNamesFromDataframe( "mtl" , demog, exclusions=xcl),
          getNamesFromDataframe( "cerebellum" , demog, exclusions=c(xcl,"_cerebell")),
          getNamesFromDataframe( "T1Hier_" , demog, exclusions=c("hier_id","[.]1","[.]2","[.]3",'background','tissue','dktregions','T1Hier_resnetGrade','hemisphere','lobes','SNR','evr','area',xcl)),
          t1namesbst,
          getNamesFromDataframe( "rsfMRI_fcnxpro" , demog, exclusions=c("hier_id",'background','thk','area','vol','FD','dvars','ssnr','tsnr','motion','SNR','evr','_alff','falff_sd','falff_mean',xcl)),
          getNamesFromDataframe( "perf_" , demog, exclusions=c("hier_id",'background','thk','area','vol','FD','dvars','ssnr','tsnr','motion','SNR','evr','_alff','falff_sd','falff_mean',xcl)),
          getNamesFromDataframe( "DTI_" , demog, exclusions=c("hier_id",'background','thk','area','vol','motion','FD','dvars','ssnr','tsnr','SNR','evr','cnx','relcn',xcl)) )
  testnames = unique( testnames )
  testnames = intersect( testnames, colnames(demog))
  if ( return_colnames ) return( testnames )

  if ( FALSE ) {
    testnames = c(
                getNamesFromDataframe( "Asym" , demog ),
                getNamesFromDataframe( "LRAVG" , demog ) ) %>% unique()
    testnames = testnames[ -multigrep( c("DTI_relcn_LRAVG","DTI_relcn_Asym"), testnames ) ]
    # special LR avg for falff
    falffnames = getNamesFromDataframe( c("falff"), demog, exclusions=c('mean','sd','Unk'))
  }

  tempnames=colnames(demog)
  tempnames=gsub("Right","right",tempnames)
  tempnames=gsub("Left","left",tempnames)
  colnames(demog)=tempnames

  if ( doasym )
    demog=mapAsymVar( demog, 
              testnames[ grep("_l_", testnames) ], '_l_', "_r_" )
  demog=mapLRAverageVar( demog, 
              testnames[ grep("_l_", testnames) ], '_l_', "_r_" )
  if ( doasym )
    demog=mapAsymVar( demog, 
                  testnames[ grep("left", testnames) ] )
  demog=mapLRAverageVar( demog, 
              testnames[ grep("left", testnames) ] )
  return( demog )
  }


#' Interpret ICC Value
#'
#' Interprets the Intra-Class Correlation Coefficient (ICC) value based on 
#' the criteria set by Cicchetti (1994) or Koo (2015).
#'
#' @param icc Numeric value; the ICC value to interpret.
#' @param criterion Character string; either "Cicchetti" or "Koo", 
#' indicating which set of criteria to use for interpretation. 
#' Defaults to "Cicchetti".
#'
#' @return A character string indicating the level of agreement or 
#' reliability as per the chosen criterion.
#'
#' @examples
#' interpret_icc(0.65, criterion = "Cicchetti")
#' interpret_icc(0.65, criterion = "Koo")
#'
#' @export
interpret_icc <- function(icc, criterion = c("Cicchetti", "Koo")) {
  criterion <- match.arg(criterion)
  
  # Cicchetti's criteria
  if (criterion == "Cicchetti") {
    if (icc < 0.40 & icc >= 0 ) {
      return("Slight")
    } else if (icc < 0.60 & icc >= 0) {
      return("Fair")
    } else if (icc < 0.75 & icc >= 0) {
      return("Moderate")
    } else if (icc <= 1.00 & icc >= 0 )  {
      return("Substantial")
    } else return(NA)
  } 
  # Koo's criteria
  else if (criterion == "Koo") {
    if (icc < 0.50 & icc >= 0 ) {
      return("Slight")
    } else if (icc < 0.75 & icc >= 0) {
      return("Fair")
    } else if (icc < 0.90 & icc >= 0 ) {
      return("Moderate")
    } else if (icc <= 1.00 & icc >= 0) {
      return("Substantial")
    } else return(NA)
  }
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
#' @importFrom nmfbin nmfbin
#' @importFrom MASS ginv
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
      if ( length( cols) == 1 ) {
        subset_df1=data.frame(df1[subset_indices, cols] )
        colnames(subset_df1)=cols
      } 
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

  library(stringr)
  library(purrr)
  replace_values <- function(input_string) {
    # Function to modify a number based on the specified rules
    modify_value <- function(number) {
      num <- as.numeric(number)
      if (num >= 1 && num <= 249) {
        return(as.character(num + 249))
      } else {
        return(number)
      }
    }
    
    # Extract all numbers from the string
    numbers <- str_extract_all(input_string, "\\b\\d+\\b")[[1]]
    
    # Apply the modification to the numbers
    modified_numbers <- map_chr(numbers, modify_value)
    
    # Replace old numbers with new numbers in the string
    for (i in seq_along(numbers)) {
      input_string <- str_replace(input_string, numbers[i], modified_numbers[i])
    }

    return(input_string)
  }

  rightvar =  gsub( leftname, rightname, leftvar )
#  for ( k in 1:length(rightvar) ) {
#    r=rightvar[k]
#    if ( length( grep("rsfMRI_",r) > 0 ) )
#      rightvar[k]=replace_values(r)
#  }
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
#' @param allowMissing boolean
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
  group,
  allowMissing=FALSE
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
  if ( ! allowMissing ) {
    mxdfin[ , adjustedoutcome ] = mxdfin[,outcomevar] - predvol 
  } else {
    wtoadjust = !is.na(predvol)
    mxdfin[ , adjustedoutcome ] = mxdfin[,outcomevar]
    mxdfin[ wtoadjust, adjustedoutcome ] = mxdfin[wtoadjust,outcomevar] - predvol[wtoadjust] 
  }
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

#' Scale Numeric Variables in a Data Frame
#'
#' This function scales the numeric variables in a data frame that are 
#' specified in a given equation. Non-numeric variables mentioned in the 
#' equation are ignored. Scaling is performed by centering and dividing by the 
#' standard deviation.
#'
#' @param mydf A data frame containing the variables to be scaled.
#' @param myeq An equation in the form of a string or an object of class 
#'   \code{formula} specifying the variables to be scaled. Only the variables 
#'   on the right-hand side of the equation are considered.
#' @param variables_to_exclude vector of strings not to scale
#'
#' @return A data frame with the specified numeric variables scaled. The 
#'   function modifies the input data frame in place and returns it.
#'
#' @examples
#' data(iris)
#' scaled_iris <- scale_variables_in_equation(iris, myeq = "Sepal.Length ~ Sepal.Width + Petal.Length")
#' head(scaled_iris)
#'
#' @export
scale_variables_in_equation <- function( mydf, myeq, variables_to_exclude ) {
  myterms = all.vars(as.formula(myeq))[-1]
  if ( ! missing( variables_to_exclude ) )
    myterms = myterms[ ! ( myterms %in% variables_to_exclude )]
  myterms = intersect(myterms, colnames(mydf))
  for ( x in myterms ) {
    if ( is.numeric( mydf[,x] )) {
      mydf[,x] = c(scale(mydf[,x]))
    }
  }
  return(mydf)
}



#' Determine the Variable Type Based on Its Name
#'
#' This function inspects the input character vector for specific patterns
#' indicating the type of variable (e.g., "T1", "rsfMRI", "DTI", "NM2") and
#' returns a corresponding string identifier for the first matched type.
#' If no known pattern is found, it returns `NA`.
#'
#' @param x A character vector containing names or identifiers to be checked
#'          against known patterns for variable types. It must be a character vector.
#'
#' @return A character string indicating the type of variable matched based
#'         on the predefined patterns ("T1", "rsfMRI", "DTI", "NM2DMT").
#'         Returns `NA` if no pattern matches.
#'
#' @examples
#' antspymm_vartype("This is a T1 weighted image")  # Returns "T1"
#' antspymm_vartype("Subject underwent rsfMRI")    # Returns "rsfMRI"
#' antspymm_vartype("DTI sequence")                # Returns "DTI"
#' antspymm_vartype("Analysis of NM2")             # Returns "NM2DMT"
#' antspymm_vartype("Unknown type")                # Returns NA
#'
#' @note This function only checks for the first occurrence of a pattern
#'       and does not account for multiple different patterns within the
#'       same input. The order of pattern checking is fixed and may affect
#'       the result if multiple patterns are present.
#'
#' @export
antspymm_vartype <- function(x) {
  # Validate input
  if (!is.character(x)) {
    stop("Input must be a character vector.")
  }
  
  # Define patterns and corresponding returns in a named list
  patterns <- list(
    T1 = "T1", rsfMRI = "rsfMRI", DTI = "DTI", NM2 = "NM2DMT", T2="T2Flair", 
    t1 = "T1", rsfmri = "rsfMRI", dti = "DTI", nm2 = "NM2DMT", t2="T2Flair", 
    t1 = "T1", rs = "rsfMRI", dt = "DTI", nm2 = "NM2DMT", t2="T2Flair", 
    perf='perf')
  
  # Iterate through the patterns
  for (pattern in names(patterns)) {
    if (any(grepl(pattern, x))) {
      return(patterns[[pattern]])
    }
  }
  
  # Default return if no pattern matched
  return(NA)
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
        if ( nrow( mycoffssub ) > 0 ) {
          myzsub = as.numeric( effectsize::z_to_d( myzsub, nrow(featureMatrix))$d  )
        } else myzsub=myz
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
  glist[[length(glist)+1]]= ggscatter(indf, x = xvar, y = yvar, color=colorvar,   size=3.45,  point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ylim( ylimmer )
  # ggtitle(paste(anatshow[1])) + ylim( ylimmer ) #+ theme(legend.position = "none")

  if ( is.numeric(indf[,anat]) ) {
    medsplit = median( indf[,anat], na.rm=T )
    hisel = indf[,anat] > medsplit
    loclev=loclev2=anat
  } else {
    loclev = (unique(indf[,anat])[1])
    hisel=indf[,anat]==loclev
    loclev2=paste0("!",loclev)
  }
  glist[[length(glist)+1]]=ggscatter(indf[hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE ) + theme(text = element_text(size=12))+ ggtitle(paste('+',anatshow[2])) + theme(legend.position = "none") + ylim( ylimmer )

  
  glist[[length(glist)+1]]=ggscatter(indf[!hisel,], x = xvar, y = yvar, color=colorvar,   size=3.45, point=showpoints, add = "reg.line", palette=ggpalette, conf.int=T, cor.coef=TRUE  ) + theme(text = element_text(size=12))+ ggtitle(paste('-',anatshow[3]))+ theme(legend.position = "none") + ylim(  ylimmer )

  grid.arrange(grobs=glist,ncol=3,top=anatshow[1])

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

#' Augment a Data Frame with Custom Color Column
#'
#' This function adds a new column named 'custom_color' to a given data frame.
#' The new column maps factor levels of a specified column to colors
#' from a user-defined palette.
#'
#' @param df A data frame to be augmented.
#' @param column_name The name of the column in `df` whose factor levels are to be 
#'   mapped to colors. This column should exist in `df` and should be a factor or
#'   convertible to a factor.
#' @param color_palette A vector of colors, corresponding to the factor levels of
#'   the specified column. The number of colors must match the number of unique
#'   factor levels in the column.
#'
#' @return A data frame identical to `df` but with an additional column 
#'   'custom_color', which contains the color mappings.
#'
#' @examples
#' data <- data.frame(
#'   category = c("A", "B", "C", "A", "B"),
#'   value = c(10, 20, 30, 40, 50)
#' )
#' palette <- c("red", "green", "blue")
#' augmented_data <- augment_with_custom_color(data, "category", palette)
#'
#' @importFrom ggpubr get_palette
#' @importFrom stats setNames
#' @export
augment_with_custom_color <- function(df, column_name, color_palette, color_palette_name = NULL ) {
  # Check if the column exists in the data frame
  if (!column_name %in% names(df)) {
    stop("Column not found in the data frame")
  }

  # Convert the column to a factor if it's not already
  factor_col <- factor(df[[column_name]])

  # Extract the factor levels of the column
  factor_levels <- levels(factor_col)

  if ( ! missing( color_palette_name ) ) {
    color_palette = get_palette(palette = color_palette_name, length(factor_levels))
  }

  # Check if the number of colors matches the number of factor levels
  if (length(color_palette) != length(factor_levels)) {
    stop("The number of colors must match the number of factor levels")
  }

  # Create a named vector for color mapping
  color_mapping <- setNames(color_palette, factor_levels)

  # Map the factor levels to the colors
  df$custom_color <- color_mapping[factor_col]

  return(df)
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
#' @param extradata dataframe with additional information to be rendered
#' @return the quantile transformed vector
#' @author Avants BB
#' @export
prplot  <- function(
  mdl, xvariable, byvariable, titlestring='', ystring='', addpoints=0, palette='npg', colorvar='', extradata=NULL ) {
  addthepoints=FALSE
  colorvarnotempty=TRUE
  if ( colorvar=='') {
    colorvarnotempty=FALSE
    colorvar='black'
  }
  if ( !is.null( extradata ) ) {
    extradata=extradata[ names(predict(mdl)), ]
  }
  if ( addpoints > 0 ) addthepoints=TRUE
  if ( ! missing( byvariable ) ) {
    vv=visreg::visreg( mdl, xvariable, by=byvariable, plot=FALSE)
    if ( colorvar %in% model.frame(mdl) ) {
        vv$res[,colorvar]=model.frame(mdl)[ names(predict(mdl)), colorvar ]
      } else if ( !is.null( extradata ) ) {
        if ( colorvar %in% colnames( extradata ))
          vv$res[,colorvar]=extradata[,colorvar]
      }
#    vv$res=augment_with_custom_color(vv$res, colorvar, color_palette_name=palette)
    if ( is.factor(vv$res[,xvariable] ) | is.character(vv$res[,xvariable]) ) {
      return( ggdotplot(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    conf.int=T,
                    point=addthepoints,
                    facet.by=byvariable,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
    } else {
      myaddparams = list(color = "blue", fill = "cyan")
      myaddparams = list()
      mypp = predict( mdl )
      mylims = range(mypp) * c( 1.3, 1. )
      return( ggscatter(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    point=addthepoints, add='reg.line', conf.int=T,
                    color=colorvar, facet.by=byvariable,
                    add.params = myaddparams,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) +
                    theme(legend.position = "top", legend.title = element_blank())  )
    }
  }
  if ( missing( byvariable ) ) {
    vv=visreg::visreg( mdl, xvariable, plot=FALSE)
    if ( colorvar %in% model.frame(mdl) ) {
      vv$res[,colorvar]=model.frame(mdl)[ names(predict(mdl)), colorvar ]
    } else if ( !is.null( extradata ) ) {
      if ( colorvar %in% colnames( extradata )) {
        vv$res[,colorvar]=extradata[,colorvar]
        }
    }
#    vv$res=augment_with_custom_color(vv$res, colorvar, color_palette_name=palette)
    if ( is.factor(vv$res[,xvariable] ) | is.character(vv$res[,xvariable]) ) {
      return( ggboxplot(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, palette=palette,
                    conf.int=T,
                    point=addthepoints,
                    add.params = list(color = "blue", fill = "cyan"),
                    fill=colorvar,
                    cor.coef=TRUE ) +  
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) )
    } else {
      return( ggscatter(vv$res, x = xvariable, y = 'visregRes', 
                    size=addpoints, 
                    point=addthepoints, add='reg.line', conf.int=T,
                    color=colorvar,
                    add.params = list(color = "blue", fill = "cyan"),
                    palette=palette,
                    cor.coef=TRUE ) +
                    theme(text = element_text(size=12))+ ylab(ystring) + 
                    ggtitle( titlestring ) + theme(legend.position = "top", legend.title = element_blank())  # Position legend at top
                  )
    }
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
  x[,variable]=as.character( x[,variable] )
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



#' return nuisance variable strings
#' 
#' these strings can be used with getNamesFromDataFrame or multigrep to 
#' either include or exclude nuisance variables.
#' 
#' @export
antspymm_nuisance_names <-function(x){
xcl = c("snr_","bandp","_mean","censor","smooth","outlier","motion","FD","despik","_nc_","_evr","minut","left","right","paramset","_sd","upsampling","mhdist","RandBasis","templateL1")
return( xcl )
}

#' shorter antspymm names
#' 
#' @export
shorten_pymm_names <-function(x){
    xx=tolower(x)
    xx=gsub( "sagittal_stratum_.include_inferior_longitidinal_fasciculus_and_inferior_fronto.occipital_fasciculus..","ilf_and_ifo",xx,fixed=TRUE)
    xx=gsub("_.cres..stria_terminalis_.can_not_be_resolved_with_current_resolution..","",xx,fixed=TRUE)
    xx=gsub("_",".",xx)
    xx=gsub("longitudinal.fasciculus",'l.fasc',xx,fixed=TRUE)
    xx=gsub("corona.radiata",'cor.rad',xx,fixed=TRUE)
    xx=gsub("central",'cent',xx,fixed=TRUE)
    xx=gsub("deep.cit168",'dp.',xx,fixed=TRUE)
    xx=gsub("cit168",'',xx,fixed=TRUE)
    xx=gsub(".include",'',xx,fixed=TRUE)
    xx=gsub("mtg.sn",'',xx,fixed=TRUE)
    xx=gsub("brainstem",'.bst',xx,fixed=TRUE)
    xx=gsub("rsfmri.",'rsf.',xx,fixed=TRUE)
    xx=gsub("dti.mean.fa.",'dti.fa.',xx,fixed=TRUE)
    xx=gsub("perf.cbf.mean.",'cbf.',xx,fixed=TRUE)
    xx=gsub(".jhu.icbm.labels.1mm",'',xx,fixed=TRUE)
    xx=gsub("..include.optic.radiation..",'',xx,fixed=TRUE)
    xx=gsub("..",'.',xx,fixed=TRUE)
    xx=gsub("..",'.',xx,fixed=TRUE)
    xx=gsub("cerebellar.peduncle",'cereb.ped',xx,fixed=TRUE)
    xx=gsub("anterior.limb.of.internal.capsule",'ant.int.cap',xx,fixed=TRUE)
    xx=gsub("posterior.limb.of.internal.capsule",'post.int.cap',xx,fixed=TRUE)
    xx=gsub("t1hier.",'t1.',xx,fixed=TRUE)
    xx=gsub("anterior",'ant',xx,fixed=TRUE)
    xx=gsub("posterior",'post',xx,fixed=TRUE)
    xx=gsub("inferior",'inf',xx,fixed=TRUE)
    xx=gsub("superior",'sup',xx,fixed=TRUE)
    xx=gsub("dktcortex",'.ctx',xx,fixed=TRUE)
    xx=gsub(".lravg",'',xx,fixed=TRUE)
    xx=gsub("dti.mean.fa",'dti.fa',xx,fixed=TRUE)
    xx=gsub("retrolenticular.part.of.internal","rent.int.cap",xx,fixed=TRUE)
    xx=gsub("iculus.could.be.a.part.of.ant.internal.capsule","",xx,fixed=TRUE)
    xx=gsub("iculus.could.be.a.part.of.ant.internal.capsule","",xx,fixed=TRUE)
    xx=gsub(".fronto.occipital.",".frnt.occ.",xx,fixed=TRUE)
    xx=gsub(".longitidinal.fasciculus.",".long.fasc.",xx,fixed=TRUE)
    xx=gsub(".longitidinal.fasciculus.",".long.fasc.",xx,fixed=TRUE)
    xx=gsub(".external.capsule",".ext.cap",xx,fixed=TRUE)
    xx=gsub("of.internal.capsule",".int.cap",xx,fixed=TRUE)
    xx=gsub("fornix.cres.stria.terminalis","fornix.",xx,fixed=TRUE)
    xx=gsub("capsule","",xx,fixed=TRUE)
    xx=gsub("and.inf.frnt.occ.fasciculus.","",xx,fixed=TRUE)
    xx=gsub("crossing.tract.a.part.of.mcp.","",xx,fixed=TRUE)
    xx=gsub("post.thalamic.radiation.optic.radiation","post.thalamic.radiation",xx,fixed=TRUE)
    xx=gsub("adjusted",'adj',xx,fixed=TRUE)
    xx=gsub("..",'.',xx,fixed=TRUE)
    xx=gsub("t1w.mean","t1vth",xx,fixed=TRUE)
    xx=gsub("fcnxpro129","p2",xx,fixed=TRUE)
    xx=gsub("fcnxpro134","p3",xx,fixed=TRUE)
    xx=gsub("fcnxpro122","p1",xx,fixed=TRUE)
#    for ( x in 1:length(xx) ) {
#      xx[x]=substr(xx[x],0,40)
#    }
    return(xx)
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



#' fill1col2another
#' 
#' fill missing data values into one column from another column;
#' will only fill values that are NA
#'
#' @param df dataframe input
#' @param x column to fill
#' @param y column to fill from
#' @return filled dataframe
#' @author Avants BB
#' @export
fill1col2another <- function( df, x, y ) {
    losel = is.na( df[,x] ) & !is.na(  df[,y] )
    losel[ is.na( losel )] = FALSE
    if ( sum( losel ) > 0 )
      df[ losel, x ] = df[ losel, y ]
    return( df )
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


#' consensusSubtypingPrep
#' 
#' generate input for consensus clustering give several clustering algorithms
#'
#' @param dataToTrain dataframe input that contains the relevant variables (may have others as well) on which training will be based
#' @param dataToPredict dataframe input that contains the relevant variables (may have others as well) for which prediction will be done
#' @param featureNames names to use in the clustering
#' @param clustVec names of the clustering methods to use
#' @param maxK the maximum desired number of classes
#' @param reorderingVariable the name of the column to use to reorder the cluster names
#' @param mvcl character prefix for the new cluster column names
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @param whichrank allows user to get 2nd (or 3rd) rank set of methods
#' @param ntoreturnperk number of method results per k to return
#' @param verbose boolean
#' @return new dataframe with new variables attached
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @importFrom caret dummyVars contr.ltfr
#' @export
consensusSubtypingPrep = function( dataToTrain, dataToPredict, featureNames, clustVec, maxK, 
  reorderingVariable, mvcl, idvar, visitName, baselineVisit, whichrank=0, ntoreturnperk=2, verbose=FALSE ) {
  ############################
  dataToPredictAll = dataToPredict
  for ( desiredk in 2:maxK ) {
    mvclK=paste0(mvcl,"_cstk_",desiredk)
    ## cluster method training
    cclusterobj = consensusSubtypingTrain( 
        dataToTrain, featureNames, clustVec, desiredk, 
        reorderingVariable=reorderingVariable, mvcl=mvclK, verbose=verbose )
    dataToTrain = cclusterobj[["newdata"]]
    clustmodels = cclusterobj[["models"]]
    reoModels = cclusterobj[["reorderers"]]
    clustlist = names( clustmodels )
    myClusterScores = cclusterobj[["ClusterSimilarityScores"]]
    myClusterCluster = cclusterobj[["ClusterClusters"]]
    clustind = which( myClusterScores == 
        sort(myClusterScores)[length(myClusterScores) - whichrank  ])# 2nd best
    moreagreeable = names( myClusterCluster[ myClusterCluster == clustind ] )
    if ( verbose ) {
      message("moreagreeable")
      print(moreagreeable)
      }
    moreagreeable = head( names( myClusterScores[ order(myClusterScores,decreasing=T) ] ),
      ntoreturnperk )
    if ( verbose ) {
      message("mostagreeable")
      print(moreagreeable)
      }
    # first run the clustering on all PPMI data 
    clustdata = consensusSubtypingPredict( dataToPredict, 
        featureNames, moreagreeable, clustmodels, 
        reoModels, mvclK, 
        idvar=idvar, visitName=visitName, baselineVisit=baselineVisit )
    locclustnames = getNamesFromDataframe( mvclK, clustdata, exclusions="_mem" )
    dataToPredictAll[, locclustnames] = clustdata[, locclustnames ]
    }
    return( dataToPredictAll )
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
#' @param maxK maximum number of clusters
#' @param consensusmethod either kmeans or hclust
#' @param returnonehot boolean 
#' @param binnmf integer (0 is default - no nmf) 
#' @param verbose boolean
#' @return new dataframe with new variables attached
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @importFrom caret dummyVars contr.ltfr
#' @export
consensusSubtypingCOCA = function( dataToClust, targetk, cocanames, newclustername, reorderingVariable, idvar, visitName, baselineVisit, maxK, consensusmethod='kmeans', returnonehot=FALSE, binnmf=0, verbose=TRUE ) {
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
    if ( binnmf > 1 ) {
      dmytx=nmfbin::nmfbin(data.matrix(dmytx), binnmf )$W
      if ( verbose ) print("nmfbin done")
    }
    if ( returnonehot ) return( dmytx )
    if ( ! missing( targetk ) & missing( maxK ) ) {
      cocatx = coca::coca(dmytx, K = targetk, B=1000, maxIterKM=5000, ccClMethod=consensusmethod )
    } else if ( ! missing( maxK ) ) {
      message(paste(consensusmethod,'spearman',maxK))
      cocatx = coca::coca(dmytx, K = NULL, maxK=maxK, 
#        dunn2s=TRUE,
        choiceKmethod='AUC', 
#        widestGap=TRUE, ccDistHC='spearman',
        B=1000, maxIterKM=5000, ccClMethod=consensusmethod )
    } else stop("Must set either maxK or targetk")
    if ( verbose ) message("COCA complete")
#    cocatx = coca::coca(dmytx, maxK = 6, B=5000 )
#    coca = coca::coca( dmytx, maxK = 10, hclustMethod = "average")
    cocatxlab = as.numeric( cocatx$clusterLabels )
    if ( verbose ) {
      print( table( cocatxlab ) )
      }
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


#' consensusSubtypingPAM
#' 
#' apply consensus clustering give several clustering solutions using PAM as the consensus method
#'
#' @param dataToClust dataframe input that contains the relevant variables (may have others as well)
#' @param targetk the desired number of classes
#' @param cocanames names of columns to use for the consensus
#' @param newclustername the column name for the consensus clustering
#' @param reorderingVariable the name of the column to use to reorder the cluster names
#' @param idvar variable name for unique subject identifier column
#' @param visitName the column name defining the visit variables
#' @param baselineVisit the string naming the baseline visit
#' @param maxK maximum number of clusters
#' @param criterion see fpc pamk
#' @param verbose boolean
#' @return new dataframe with new variables attached
#' @author Avants BB
#' @examples
#' mydf = generateSubtyperData( 100 )
#' @export
consensusSubtypingPAM = function( dataToClust, targetk, cocanames, newclustername, reorderingVariable, idvar, visitName, baselineVisit, maxK, criterion='asw', verbose=TRUE ) {
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
    if ( ! missing( targetk ) & missing( maxK ) ) {
      cocatx = fpc::pamk( dmytx, targetk, usepam=FALSE, criterion=criterion )
    } else if ( ! missing( maxK ) ) {
      cocatx = fpc::pamk( dmytx, 2:maxK, usepam=FALSE, criterion=criterion )
    } else stop("Must set either maxK or targetk")
    if ( verbose ) message("COCA complete")
    cocatxlab = cocatx$pamobject$clustering
    if ( verbose ) {
      print( table( cocatxlab ) )
      }
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


#' Process Clinical/Demographic and imaging Data for an ADNI Study
#'
#' This function merges two data frames based on a common patient ID and the closest date, ensuring that each row in the first data frame (`dfA`) is matched with the row from the second data frame (`dfB`) that has the closest date for the same patient ID. The final merged data frame includes all columns from both `dfA` and `dfB`, excluding the patient ID and date columns from `dfB` to avoid duplication. The function is designed to handle date columns as Date objects and includes a progress bar to indicate the matching process's progress. EXAMDATE is assumed present in dfA and date present in dfB where date is numerically formatted as YYYYMMDD.
#'
#' @param dfA The first data frame to be merged, expected to contain columns for patient ID and date.   may need \code{dfA$subjectID = dfA$PTID}.
#' @param dfB The second data frame to be merged, expected to contain columns for patient ID and date. Rows from `dfB` are matched to `dfA` based on the closest date for each patient ID.
#' @param patientidcol Character string specifying the column name in both data frames that contains the patient ID. Default is 'subjectID'.
#' @param verbose Logical indicating whether to print additional information during processing, such as the number of common IDs found and dimensions of the data frames before and after merging. Default is TRUE.  setting verbose to 2 provides more feedback.
#'
#' @return Returns a merged data frame with the same number of rows as `dfA` and includes all columns from both `dfA` and `dfB`, with `dfB` columns matched based on the closest date for each patient ID. The date columns from `dfB` are excluded to avoid duplication.
#'
#' @examples
#' # Assuming dfA and dfB are already defined and have 'subjectID' and 'date' columns
#' # merged_df <- merge_ADNI_antspymm_by_closest_date(dfA, dfB)
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#'
#' @note The function assumes that the date columns in the dfB data frame formatted as 'YYYYMMDD' and converts them to Date objects for processing. The progress bar functionality uses base R's txtProgressBar, which is displayed in the console.  dfA's EXAMDATE is of the form YYYY-MM-DD.
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
merge_ADNI_antspymm_by_closest_date <- function(dfA, dfB, patientidcol='subjectID', verbose=TRUE) {
  # Safety checks for required columns in dfA and dfB
  if (!all(c(patientidcol, 'EXAMDATE') %in% names(dfA))) {
    stop("dfA is missing one of the required columns: ", patientidcol, " or EXAMDATE ... try setting dfA$subjectID = dfA$PTID ")
  }
  if (!all(c(patientidcol, 'date') %in% names(dfB))) {
    stop("dfB is missing one of the required columns: ", patientidcol, " or date")
  }
  
  # Attempt to convert dfB's date to Date object
  tryCatch({
    dfB[,'EXAMDATE'] <- as.Date(as.character(dfB[,'date']), format="%Y%m%d")
  }, error = function(e) {
    stop("Error in converting 'date' column in dfB to Date format: ", e$message)
  })
  
  # Intersect patient IDs to ensure matching on patientidcol
  isect <- intersect(dfA[, patientidcol], dfB[, patientidcol])
  if (verbose) {
    print(paste(length(isect), "common ids"))
  }
  dfA <- dfA[dfA[, patientidcol] %in% isect, ]
  dfB <- dfB[dfB[, patientidcol] %in% isect, ]
  
  if (verbose) {
    print("Initial dfA dimensions")
    print(dim(dfA))
    print(head(dfA[, 'EXAMDATE'], 3))
    print(head(dfB[, 'EXAMDATE'], 3))
  }
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nrow(dfA), style = 3)
  
  matched_rows <- vector("list", nrow(dfA))
  date_diff_years <- numeric(nrow(dfA))
  
  for (i in 1:nrow(dfA)) {
    current_row <- dfA[i, ]
    matching_rows <- dfB[dfB[, patientidcol] == current_row[, patientidcol], ]
    
    if (nrow(matching_rows) > 0) {
      date_diffs <- abs(difftime(matching_rows[,'EXAMDATE'], current_row[,'EXAMDATE'], units = "days"))
      date_diffs <- as.numeric(date_diffs) / 365.25
      closest_row_index <- which.min(date_diffs)
      matched_rows[[i]] <- matching_rows[closest_row_index, ]
      date_diff_years[i] <- min(date_diffs)
    } else {
      matched_rows[[i]] <- setNames(as.list(rep(NA, ncol(dfB))), names(dfB))
      date_diff_years[i] <- NA
    }
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  dfA$date_diff_years <- date_diff_years
  
  matched_dfB <- do.call(rbind, matched_rows)
  merged_df <- cbind(dfA, matched_dfB[, !(names(matched_dfB) %in% c('EXAMDATE', patientidcol))])
  
  if (verbose) {
    print("Final merged_df dimensions")
    print(dim(merged_df))
  }

  return(merged_df)
}


#' Process Clinical/Demographic and imaging Data for a PPMI Study
#'
#' This function merges, cleans, and enriches clinical and demographic data from a neurological study.
#' It checks for required columns, replaces placeholder values, merges datasets, and computes new variables.
#'
#' @param demog DataFrame containing PPMI demographic information with columns like PATNO, BIRTHDT.  An example filename from LONI would be Demographics_06Feb2024.csv
#' @param ppmidemog0 DataFrame containing curated PPMI clinical data including PATNO, EVENT_ID, age_at_visit, age.  An example filename from LONI would be PPMI_Curated_Data_Cut_Public_20230612_rev.csv.
#' @param pymf DataFrame containing imaging and additional clinical data keyed by subjectID, date and filename
#' @param pymversion A character string identifying the ANTsPyMM pipeline variant.
#' @param saa DataFrame (optional) containing supplemental clinical measurements with PATNO, EVENT_ID and CSFSAA.
#' @param verbose boolean
#' @return A processed DataFrame with merged and enriched clinical and demographic information.
#' @export
#' @examples
#' \dontrun{
#'   processed_data <- process_clinical_demographic_data(demog, ppmidemog0, saa, pymf)
#' }
merge_ppmi_imaging_clinical_demographic_data <- function(demog, ppmidemog0, pymf, pymversion, saa, verbose=TRUE ) {
  # Load required libraries
  library(dplyr)
  library(forcats)
  if ( missing( saa ) ) {
    saa = ppmidemog0[,c("PATNO","EVENT_ID","CSFSAA")]
  }

  # Ensure required columns are present
  stopifnot(all(c("PATNO", "BIRTHDT") %in% names(demog)),
            all(c("PATNO", "EVENT_ID", "age_at_visit", "age") %in% names(ppmidemog0)),
            all(c("PATNO", "EVENT_ID") %in% names(saa)),
            all(c("subjectID", "projectID", "filename","date","T1Hier_resnetGrade") %in% names(pymf)))

  demog[demog=='.']=NA
  saa[saa=='.']=NA
  ppmidemog0[ppmidemog0=='.']=NA
  ppmidemog0 = ppmidemog0[ , !(colnames(ppmidemog0) %in% "CSFSAA")]
  ppmi = merge( ppmidemog0, saa, by=c("PATNO","EVENT_ID"), suffixes = c("",".y"), all.x=TRUE )
  ppmi[ppmi=='.']=NA
  ppmi$yearsbl = ppmi$age_at_visit - ppmi$age
  rm( ppmidemog0 )
  idcols = c("subjectID","date","subjectIDdate", "imageID",'filename')
  idcolsdemog = c("PATNO","visit_date","EVENT_ID")
  clin2 = pymf
  clin2$PATNO = pymf$subjectID

  udx = sort(unique(ppmi$subgroup))
  dxmapper = data.frame( ppmisubgroup=udx )
  rownames(dxmapper)=udx
  for ( gdx in c("GBA","LRRK2","PINK1","SNCA","PRKN") )
    dxmapper[ udx[grep(gdx,udx)],c('genetictype')]=gdx

  dxmapper[is.na(dxmapper$genetictype),'genetictype']='Sporadic'
  dxmapper["Healthy Control",'genetictype']='CN'
  dxmapper["GBA",'jdx']='PDGBA'         # 1
  dxmapper["GBA + Hyposmia",'jdx']='ProdromalGBA'         # 2
  dxmapper["GBA + RBD + Hyposmia",'jdx']='ProdromalGBA'   # 3
  dxmapper["Healthy Control",'jdx']='CN' # 4
  dxmapper["Hyposmia",'jdx']='ProdromalSporadic' # 5
  dxmapper["LRRK2",'jdx']='PDLRRK2'     # 6
  dxmapper["LRRK2 + GBA",'jdx']='PDLRRK2'     # 7
  dxmapper["LRRK2 + GBA + Hyposmia",'jdx']='ProdromalLRRK2'     # 8
  dxmapper["LRRK2 + Hyposmia",'jdx']='ProdromalLRRK2'     # 9
  dxmapper["PINK1",'jdx']='PDPINK1' # 10
  dxmapper["PRKN",'jdx']='PDPRKN'   # 11
  dxmapper["RBD",'jdx']='ProdromalSporadic' # 12
  dxmapper["RBD + Hyposmia",'jdx']='ProdromalSporadic' # 13
  dxmapper["SNCA",'jdx']='PDSNCA' # 14
  dxmapper["SNCA + Hyposmia",'jdx']='ProdromalSNCA' # 15
  dxmapper["Sporadic",'jdx']='PDSporadic'   # 16
  ####################################################################
  clin2ppmi=unique( clin2$PATNO[ clin2$projectID == 'PPMI' ] )
  ppmidx = ppmi[ , c("PATNO", 'COHORT', 'CONCOHORT', "subgroup")]
  ppmidx$CONCOHORT[ ppmidx$CONCOHORT == 1  ]='PD'
  ppmidx$CONCOHORT[ ppmidx$CONCOHORT == 2  ]='CN'
  ppmidx$CONCOHORT[ ppmidx$CONCOHORT == 3  ]='SWEDD'
  ppmidx$CONCOHORT[ ppmidx$CONCOHORT == 4  ]='Prodromal'
  ppmidx$COHORT[ ppmidx$COHORT == 1  ]='PD'
  ppmidx$COHORT[ ppmidx$COHORT == 2  ]='CN'
  ppmidx$COHORT[ ppmidx$COHORT == 3  ]='SWEDD'
  ppmidx$COHORT[ ppmidx$COHORT == 4  ]='Prodromal'
  paireddx = data.frame( PATNO=clin2ppmi )
  rownames(  paireddx )=clin2ppmi
  clin2$joinedDX=NA
  for ( uid in clin2ppmi ) {
          uid = as.character( uid )
          dx2=na.omit(unique( ppmidx$subgroup[ subtyper::fs(ppmidx$PATNO == uid )]))
          dx1=na.omit(as.character(unique( clin2$joinedDX[ subtyper::fs(clin2$PATNO == uid )])))
          dx1new=dxmapper[dx2,'jdx']
          seluid = subtyper::fs(ppmidx$PATNO == uid )
          if ( sum(seluid) > 0 & length(dx2) > 0 ) {
            ccdx=na.omit(unique( ppmidx$CONCOHORT[ subtyper::fs(ppmidx$PATNO == uid )]))
            ccdx2=na.omit(unique( ppmidx$COHORT[ subtyper::fs(ppmidx$PATNO == uid )]))
            if ( length( ccdx) == 0 ) ccdx=ccdx2
            stopifnot( length(ccdx)==1)
            mygroup=dxmapper[dx2,'genetictype']
            if ( dx2 == "Healthy Control" ) mygroup=""
            if ( !is.na(ccdx)) newjdx = paste0(ccdx,mygroup)
            if ( is.na(ccdx)) newjdx = paste0(ccdx2,mygroup)
            clin2$joinedDX[ subtyper::fs(clin2$PATNO == uid )]=newjdx
            paireddx[uid,c('ppmidx','ccdx', 'ccdx2', 'group', 'joinedDX')]=c(dx2,ccdx,ccdx2, mygroup, newjdx)
          }  else if ( verbose ) print(paste("missing",uid))
      }
  paireddx$isConsensus=!is.na(paireddx$ccdx)
  # table( paireddx$ppmidx, paireddx$joinedDX  )
  clin2$joinedDX[ clin2$joinedDX == 'CNGBA'] = 'CN'
  ####################################################################

  # message("HERE WE DEAL WITH STATS ASYN")
  clin2$AsynStatus=NA
  clin2$age_BL=NA
  ppmi$CSFSAA[ ppmi$CSFSAA %in% c(2)]=NA
  uids = as.character(unique( ppmi$PATNO ))
  clin2$PATNO=as.character(clin2$PATNO)
  for ( u in uids ) {
    clin2sel=clin2$PATNO == as.character(u)
    if ( sum(clin2sel) > 0 ) {
      clin2[clin2sel,'age_BL']=min(ppmi[ ppmi$PATNO == as.character(u), 'age'],na.rm=T)
      suppas = unique(ppmi[ ppmi$PATNO == as.character(u), 'CSFSAA'])
      if ( 1 %in% suppas ) {
        suppas='Positive'
      } else if ( 3 %in% suppas ) {
        suppas='PosMSA'
      } else if ( 0 %in% suppas ) {
        suppas='Negative'
      } else suppas=NA
      clin2[clin2sel,'AsynStatus']=suppas
    }
  }
#  table( is.na(clin2$age_BL))
  ####################################
  imagingdate=rep( NA, nrow( clin2 ) )
  for ( k in 1:nrow( clin2 ) ) {
    imagingdate[k] = substr( unlist(strsplit( clin2$filename[k], '-' ))[3],0,4)
    }
  clin2$imaging_year = imagingdate
  testid="102529"
  ppmi$PATNO=as.character(ppmi$PATNO)
  demog$PATNO = as.character( demog$PATNO )
  commonids = unique( intersect( ppmi$PATNO, demog$PATNO ) )
  ppmi = merge( ppmi, demog, by='PATNO', all.x=TRUE )
  rmnames=names(ppmi)[ grep("[.]y",names(ppmi))]
  ppmi=ppmi[,!(names(ppmi) %in% rmnames)]
  # names(ppmi)=gsub("[.]x","",names(ppmi))
  ppmi$birthdate=ppmi$BIRTHDT
  # convert birtdate to NRG format
  for ( k in 1:nrow(ppmi) ) {
          if ( ! is.na( ppmi$BIRTHDT[k] ) ) {
              temp=unlist(strsplit(ppmi$BIRTHDT[k], "/"))
              ppmi$birthdate[k]=as.numeric(paste0(temp[2],temp[1],15)) # add the 15th as the estimated date
          }
      }
  demog$birthdate=NA
  for ( k in 1:nrow(demog) ) {
    if ( ! is.na( demog$BIRTHDT[k] ) ) {
      temp=unlist(strsplit(demog$BIRTHDT[k], "/"))
      demog$birthdate[k]=as.numeric(paste0(temp[2],temp[1],15)) 
      # add the 15th as the estimated date
      }
  }
  # now we can compute age at imaging using imaging date in addition to BIRTHDT
  clin2$ageatimaging=NA
  clin2$EVENT_ID=NA
  clin2$clin_imaging_mismatch=TRUE
  rm(clin2add)
  for ( k in 1:nrow( clin2 ) ) {
      subjectbdate = unique( demog$birthdate[ demog$PATNO == clin2$PATNO[k] ] )
      clin2$ageatimaging[k] = as.numeric( 
              difftime( 
                  nrgDateToRDate( clin2$date[k] ), 
                  nrgDateToRDate( subjectbdate ), units='weeks' )/52.0)
      # find the closest EVENT_ID using ageviz
      agesel=which( ppmi$PATNO == clin2$PATNO[k] )
      if ( length( agesel ) == 0 | is.na( clin2$ageatimaging[k] ) ) {
          next # this is like continue
          }
      locages=ppmi[agesel,'age_at_visit']
      closest=which.min( abs(clin2$ageatimaging[k]-locages))
      wclin = agesel[closest]
      clin2$EVENT_ID[k]=ppmi[wclin,'EVENT_ID.x']
      # select clin
      if ( length(wclin) == 1 ) {
          if ( ! exists("clin2add" ) ) {
              clin2add=names(ppmi)[ !(names(ppmi) %in% names(clin2) )]
              }
          clin2[k,clin2add]=ppmi[wclin,clin2add]
          clin2$clin_imaging_mismatch[k]=FALSE
          } else {
              stop("akred - should not reach this point")
              wclin = which(ppmi$PATNO == clin2$PATNO[k] )
              # find next closest match
              if ( length(wclin) > 0 ) {
                  locclin=ppmi[wclin,]
                  closest=which.min( abs(clin2$ageatimaging[k]-locclin$age_at_visit))
                  clin2[k,clin2add]=locclin[closest,clin2add]
              }
          }
  }
  clin2$yearsbl = NA
  uids = as.character(unique( clin2$PATNO ))
  for ( u in uids ) {
    clin2sel=clin2$PATNO == as.character(u)
    if ( sum(clin2sel) > 0 ) {
      clin2[clin2sel,'age_BL']=min(clin2[clin2sel, 'ageatimaging'],na.rm=T)
    }
  }
  clin2$age_BL[ is.na( clin2$age_BL )] = clin2$age[ is.na( clin2$age_BL )]
  clin2$yearsbl = clin2$ageatimaging - clin2$age_BL
  clin2$brainVolume = rowSums( clin2[,getNamesFromDataframe( c("T1Hier","vol","hemis"), clin2 )])
  clin2$brainVolume = clin2$brainVolume/mean(clin2$brainVolume,na.rm=T)
  clin2$mrimfg[ clin2$mrimfg == ""]="Unk"
  clin2$APOE[ clin2$APOE == "" ] = "e3/e3"
  clin2$abeta = as.numeric( clin2$abeta )
  clin2$tau = as.numeric( clin2$tau )
  clin2$DXSubAsyn = NA
  ###########################################################################
  nna=!is.na( clin2$AsynStatus )
  clin2$DXSubAsyn[nna]=paste0(clin2$joinedDX[nna], clin2$AsynStatus[nna] )
  updrsnames = c( getNamesFromDataframe( 'updrs' , clin2), 
    getNamesFromDataframe( 'moca' , clin2) )
  for ( u in c(updrsnames,"duration_yrs") ) clin2[,u]=as.numeric( clin2[,u])
  clin2bl=clin2[clin2$EVENT_ID=='BL',]
  aggregate( updrs_totscore ~ DXSubAsyn, clin2bl, mean, na.rm=T )
  negcnids = clin2bl$PATNO[ subtyper::fs(clin2bl$DXSubAsyn == 'CNPositive') ]
  if ( verbose ) {
    print( table( 
      clin2bl$subgroup[ subtyper::fs(clin2bl$PATNO %in% negcnids)  ],
      clin2bl$DXSub[ subtyper::fs(clin2bl$PATNO %in% negcnids)  ] ) )
  }

  ##########################
  clin2$commonID=clin2$PATNO
  clin2$commonAge=clin2$age_BL
  clin2$commonSex="Male"
  clin2$commonSex[clin2$SEX ==0 ]="Female"
  clin2$commonEdu=as.numeric( clin2$educ )
  clin2$MOCA=clin2$moca
  clin2$studyName = clin2$projectID
  clin2$hy=as.numeric( clin2$hy )
  clin2b = fillBaselineColumn( clin2,
          c('brainVolume','hy','MOCA',updrsnames), 
          'commonID', 'yearsbl', 0, 
          fast=FALSE, verbose=F )[[1]]

  if ( verbose ) {
    print( table( 
      clin2bl$joinedDX,
      clin2bl$AsynStatus ) )
  }

  return(clin2b)
}




#' Visualize the Joint Effect of Variables in a GLM
#'
#' This function visualizes the joint effect of several variables in a generalized linear model (GLM)
#' by plotting the predicted response over a range of values for the primary predictor variable, while
#' holding other predictors at their median or mode values. It allows specifying a group for focused analysis.
#'
#' @param demogmdl Data frame containing the variables used in the GLM.
#' @param qmdl Fitted GLM object from which predictions will be generated.
#' @param x Character vector specifying the names of the predictor variables, with the first being the primary.
#' @param y The name of the response variable.
#' @param group The name of the group variable or specific group to be analyzed.
#' @param titlestring Title of the plot indicating the focus of the visualization.
#' @param groupvar Optional; the name of the variable in `demogmdl` that defines group membership. Default is 'group'.
#' @param predictorsigns Optional; a named numeric vector indicating the direction of the effect of each predictor.
#' @param jdf_simulation boolean
#' @param xrange explicitly set xrange
#' @param verbose Logical; if TRUE, additional processing information will be printed to the console.
#'
#' @return Generates a plot visualizing the predicted response and confidence intervals across the range of the primary predictor.
#'
#' @examples
#' # Assuming `data` is your dataset, `fit` is a fitted GLM, and you're interested in predictors `x1` and `x2`:
#' # visglm(data, fit, c("x1", "x2"), "y", "control", "Visualization for Control Group")
#'
#' @importFrom ggplot2 ggplot geom_line geom_point
#' @importFrom ciTools add_ci
#' @import ciTools
#' @export
visglm <- function(demogmdl, qmdl, x, y, group=NULL, titlestring='',  varstoadd=NULL, groupvar = 'group', predictorsigns=NULL, jdf_simulation=FALSE, xrange=NULL, verbose = FALSE) {  
  .env <- environment() ## identify the environment of cv.step
  if ( is.null( varstoadd ) ) {
    varstoadd=all.vars( formula( qmdl ) )
  }
  simulateJointDistribution <- function(exampleData, nSimulations, distribution='unf') {
    # Step 1: Compute covariance matrix and perform Cholesky decomposition
    covMatrix <- cov(data.matrix(exampleData))
    cholMatrix <- chol(covMatrix)
    
    # Step 2: Simulate standard normally distributed data
    n <- nrow(exampleData)  # Original number of observations
    p <- ncol(exampleData)  # Original number of variables
    if ( distribution=='normal' ) {
      simDataStdNormal <- matrix(rnorm(nSimulations * p), nrow = nSimulations, ncol = p)
    } else {
      simDataStdNormal <- scale(matrix( 1:(nSimulations * p), nrow = nSimulations, ncol = p),T,T)
    #  simDataStdNormal[ simDataStdNormal > 1.5]=1.5
    #  simDataStdNormal[ simDataStdNormal < -1.5]=-1.5
    }
    # Step 3: Apply the Cholesky matrix to simulated data to preserve correlation
    simDataCorrelated <- simDataStdNormal %*% cholMatrix
    
    # Step 4: Adjust the means of the simulated data
    meanVector <- colMeans(exampleData)  # Mean of each variable in the original data
    simDataAdjusted <- sweep(simDataCorrelated, 2, meanVector, '+')
    
    # Step 5: Convert to a data frame and assign column names
    simDataFrame <- as.data.frame(simDataAdjusted)
    names(simDataFrame) <- names(exampleData)
    
    return(simDataFrame)
  }

  verbfn <- function( x, verbose=FALSE ) {
    if ( verbose ) print(paste("visglm:", x)); 
  }
  verbfn('start',verbose)
  myMode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
    }
  mostfreq <- function( x, na.rm=TRUE ) {
    if ( is.numeric( x ) ) return( mean(x,na.rm=na.rm ) )
    return( myMode( x ) )
  }

  # Check if groupvar exists in data
  if ( !(groupvar %in% colnames(demogmdl) )) {
    errmsg = paste("The group variable", groupvar, "does not exist in the dataset.")
    stop(errmsg)
  }
  verbfn('check group',verbose)
  if (is.null(predictorsigns)) {
    predictorsigns = rep(1,length(x))
    names(predictorsigns)=x
    }
  verbfn('predictorsigns',verbose)
  xrowname=paste(x,collapse='+')
  verbfn('xrowname',verbose)
  tt=x[1]
  n2sim = 500
  timeaxis0 = seq( min(demogmdl[,tt],na.rm=T), max(demogmdl[,tt],na.rm=T),length.out=n2sim)
  if ( predictorsigns[x[1]] < 0 ) timeaxis0=rev(timeaxis0)
  verbfn('confint',verbose)
  myconf = confint(qmdl)
  verbfn('confint done',verbose)
  mycolor='magenta'
  mylev=group
  verbfn('group',verbose)
  mycolor='blue'
  if ( !is.null(group) ) {
    if ( group != 'control' ) mycolor='blue'
    if ( group=='all') {
      psel=rep(TRUE,nrow(demogmdl))
      mylev=group
    } else psel = demogmdl[,groupvar] == group
  } else psel = !is.na( demogmdl[,y] )
  atpreds = list( )
  verbfn('all.vars',verbose)
  if ( missing( varstoadd) ) {
    varstoadd = all.vars(qmdl$formula)[-1]
  }
  verbfn('varstoadd',verbose)
  if ( jdf_simulation ) {
    simmed = simulateJointDistribution( data.matrix(demogmdl[,x]), n2sim )
    colnames(simmed)=x
    }
  if ( is.null( xrange )) {
    loqhiq=range(demogmdl[,x],na.rm=T)
  } else loqhiq = xrange
  for ( zz in  varstoadd ) {
    n = length( atpreds ) + 1
    if ( verbose ) print( paste( n,"of", length(varstoadd),zz))
    if ( zz %in% x ) {
        timeaxis = seq( loqhiq[1], loqhiq[2],length.out=length(timeaxis0))
        if ( predictorsigns[zz] < 0 ) timeaxis=rev(timeaxis)
        atpreds[[ n ]] = timeaxis
        if ( jdf_simulation )
          atpreds[[ n ]] = as.numeric(simmed[,zz])
    } else if ( is.numeric( demogmdl[psel,zz] )) {
        atpreds[[ n ]] = rep( mean(demogmdl[psel,zz],na.rm=T), length(timeaxis0))
    } else {
        mostfreqvar = mostfreq( demogmdl[psel,zz] )
        atpreds[[ n ]] = rep( mostfreqvar, length(timeaxis0))
    }
    names(atpreds)[[n]] = zz
  }
  
  verbfn('x1',verbose)
  verbfn('Ypred',verbose)
  Y <- predict( qmdl, atpreds, type='response', se.fit=TRUE, env=.env )
  
  verbfn('add_ci',verbose)
  demogmdl$imaging_protocol=mostfreq(demogmdl$imaging_protocol)
  demogmdl$commonID=mostfreq(demogmdl$commonID)
  dat1 <- add_ci(demogmdl, qmdl, names = c("lpb", "upb"), alpha = 0.05, nsims = 25, allow.new.levels = TRUE)
  verbfn('Yci',verbose)
  ydelta = Y$fit - Y$se.fit
  Y$ciminus = Y$fit-1.96*Y$se.fit
  Y$ciplus = Y$fit+1.96*Y$se.fit
  verbfn('ciplus',verbose)
  myyvars = c(demogmdl[psel,y],Y$ciminus, Y$ciplus)
  myyvars = demogmdl[,y]
  verbfn('myyvars',verbose)
  rangerx = range(demogmdl[psel,x],na.rm=T)
  scl = c(0.98,1.02)
  rangerx = rangerx * scl
  rangery = range(myyvars,na.rm=T) * scl
  verbfn('ranger',verbose)
  plot(timeaxis0, Y$fit, xlab = xrowname, ylab = y, type='l',
    xlim=loqhiq,
    ylim=rangery, main=paste(group, " ", titlestring))
  verbfn('plotting',verbose)
  lines(timeaxis0, Y$fit, lwd = 2, col = mycolor)
  lines(timeaxis0, Y$ciplus, lwd = 2, col = "red", lty=2)
  lines(timeaxis0, Y$ciminus, lwd = 2, col = "red", lty=2)
  verbfn('plot',verbose)
#  myco=(coefficents(summary(qmdl)))
#  corows=rownames(myco)
#  bestx=which.min( myco[ ])
  if ( length( x ) > 1 ) {
    xpoints = rowMeans( demogmdl[psel,x], na.rm=T )
  } else xpoints = demogmdl[psel,x]
  points( xpoints, demogmdl[psel,y], col=mycolor )
  if ( !is.null( group ) )
  if ( group == 'all' ) {
    notexp=demogmdl[,groupvar]!=group
    points( (demogmdl[notexp,x]), demogmdl[notexp,y], col='magenta' )
  }
}


#' Assign Quality Control Ratings Based on Specific Criteria
#'
#' This function assigns quality control ratings to a dataset based on several
#' criteria related to measurements of a specific anatomical structure, such as the substantia nigra.
#' The ratings are determined by whether the measurements fall within certain thresholds.
#'
#' @param df A data frame containing the dataset to be rated. No default value, must be provided by the user.
#' @param z_coord_col Name of the column for the Z-coordinate measurement. Default is "NM2DMT_NM_substantianigra_z_coordinate".
#' @param volume_col Name of the column for the volume measurement. Default is "NM2DMT_NM_volume_substantianigra".
#' @param avg_col Name of the column for the average measurement. Default is "NM2DMT_NM_avg_substantianigra".
#' @param sd_col Name of the column for the standard deviation measurement. Default is "NM2DMT_NM_sd".
#' @param max_col Name of the column for the maximum value measurement. Default is "NM2DMT_NM_max".
#' @param lohi A numeric vector of length 2 specifying the lower and upper thresholds for the Z-coordinate. Defaults to c(0.3, 0.7).
#' @param volume_threshold A numeric value specifying the threshold for the volume measurement. Default is 500.
#' @param avg_threshold A numeric value specifying the threshold for the average measurement. Default is 2500.
#' @param sd_threshold A numeric value specifying the threshold for the standard deviation measurement. Default is 50.
#' @param max_threshold A numeric value specifying the threshold for the maximum value measurement. Default is 5500.
#' @param verbose  boolean
#' @return The original data frame with additional columns containing the quality control ratings. 1=pass; 0=fail.
#' @examples
#' # Example usage:
#' # df_rated <- assign_qc_ratings(df)
#' # View(df_rated)
#' @export
assign_qc_ratings_NM2DMT <- function(df, 
                              z_coord_col = "NM2DMT_NM_substantianigra_z_coordinate", 
                              volume_col = "NM2DMT_NM_volume_substantianigra", 
                              avg_col = "NM2DMT_NM_avg_substantianigra", 
                              sd_col = "NM2DMT_NM_sd", 
                              max_col = "NM2DMT_NM_max", 
                              lohi = c(0.3, 0.7), 
                              volume_threshold = 500, 
                              avg_threshold = 2500, 
                              sd_threshold = 50, 
                              max_threshold = 5500, verbose=FALSE ) {  
  df$NM_QC_Ratings_Z <- NA
  df$NM_QC_Ratings_Z[subtyper::fs(df[[z_coord_col]] > lohi[1] & df[[z_coord_col]] < lohi[2])] <- 1
  df$NM_QC_Ratings_Z[subtyper::fs(df[[z_coord_col]] <= lohi[1] | df[[z_coord_col]] >= lohi[2])] <- 0
  df$NM_QC_Ratings_SNVol <- NA
  df$NM_QC_Ratings_SNVol[subtyper::fs(df[[volume_col]] < volume_threshold)] <- 0
  df$NM_QC_Ratings_SNVol[subtyper::fs(df[[volume_col]] >= volume_threshold)] <- 1
  df$NM_QC_Ratings_AvgIntensity <- NA
  df$NM_QC_Ratings_AvgIntensity[subtyper::fs(df[[avg_col]] <= avg_threshold)] <- 1
  df$NM_QC_Ratings_AvgIntensity[subtyper::fs(df[[avg_col]] > avg_threshold)] <- 0
  df$NM_QC_Ratings_SDIntensity <- NA
  df$NM_QC_Ratings_SDIntensity[subtyper::fs(df[[sd_col]] > sd_threshold)] <- 0
  df$NM_QC_Ratings_SDIntensity[subtyper::fs(df[[sd_col]] <= sd_threshold)] <- 1
  df$NM_QC_Ratings_MaxIntensity <- NA
  df$NM_QC_Ratings_MaxIntensity[subtyper::fs(df[[max_col]] > max_threshold)] <- 0
  df$NM_QC_Ratings_MaxIntensity[subtyper::fs(df[[max_col]] <= max_threshold)] <- 1
  

  nmqccols=c("NM_QC_Ratings_Z","NM_QC_Ratings_SNVol","NM_QC_Ratings_AvgIntensity","NM_QC_Ratings_SDIntensity","NM_QC_Ratings_MaxIntensity")
  df$NM_QC_Ratings_Failures=length(nmqccols)-rowSums( df[,nmqccols], na.rm=T )
  df$NM_QC_Ratings_Failures[ is.na( df$NM2DMT_NM_max ) ]=NA
  nmqccols=c(nmqccols,'NM_QC_Ratings_Failures')
  # Optionally, return tables of ratings and their proportions
  ratings_table <- table(df$NM_QC_Ratings_Failures)
  if ( verbose ) {
    print( ratings_table )
  }
  return( df[,nmqccols] )
  # list(df = df, ratings_table = ratings_table, ratings_proportion = ratings_proportion)
}


#' Eliminate Non-Unique Columns from a Matrix
#'
#' This function takes a matrix as input and returns a new matrix with all non-unique columns removed.
#' It first transposes the matrix to work with rows instead of columns, making it easier to apply the
#' `duplicated` function. Then, it finds and removes duplicated rows in the transposed matrix, which correspond
#' to non-unique columns in the original matrix. Finally, it transposes the matrix back to its original
#' orientation, now with only unique columns remaining.
#'
#' @param matrix A numeric or character matrix from which to remove non-unique columns.
#' @return A matrix with non-unique columns removed, maintaining the original order of the remaining columns.
#' @examples
#' exampleMatrix <- matrix(c(1, 2, 3, 1, 2, 3, 1, 2, 3, 4, 5, 6), nrow = 3, byrow = TRUE)
#' print(exampleMatrix)
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    2    3    1
#' # [2,]    2    3    1    2
#' # [3,]    3    1    2    4
#' 
#' uniqueMatrix <- eliminateNonUniqueColumns(exampleMatrix)
#' print(uniqueMatrix)
#' #      [,1] [,2]
#' # [1,]    1    1
#' # [2,]    2    2
#' # [3,]    3    4
#' @export
eliminateNonUniqueColumns <- function(matrix) {
  # Transpose the matrix to work with rows instead of columns for easier application of the 'duplicated' function
  t_matrix <- t(matrix)
  
  # Find unique rows in the transposed matrix (which correspond to unique columns in the original matrix)
  unique_t_matrix <- t_matrix[!duplicated(t_matrix), ]
  
  # Transpose back to the original orientation, now with non-unique columns removed
  result_matrix <- t(unique_t_matrix)
  
  return(result_matrix)
}


#' Interpret SiMLR Vector
#'
#' This function interprets a vector from SiMLR (similarity-driven multivariate linear reconstruction)
#' specifically focusing on a given variable (e.g., a specific principal component or cluster). It extracts and normalizes the vector associated 
#' with the specified SiMLR variable, sorts it to identify the top elements, and optionally filters out non-significant values. 
#' This function is useful for understanding the contribution of different features in the context of the SiMLR analysis.
#'
#' @param simlrResult A list containing SiMLR analysis results, which should include a matrix `v` 
#' representing vectors of interest (e.g., principal components).
#' @param simlrMats A list of matrices associated with SiMLR analysis, where each matrix corresponds 
#' to a different modality or data type analyzed by SiMLR.
#' @param simlrVariable A string specifying the variable within `simlrResult` to interpret. The variable 
#' name should include both an identifier (e.g., "PC" for principal component) and a numeric index.
#' @param n2show An integer specifying the number of top elements to show from the sorted, normalized vector. 
#' Defaults to 5. If `NULL` or greater than the length of the vector, all elements are shown.
#' @param shortnames boolean
#' @param return_dataframe boolean
#' @return A named vector of the top `n2show` elements (or all if `n2show` is `NULL` or too large), 
#' sorted in decreasing order of their absolute values. Elements are named according to their identifiers 
#' in `simlrMats` and filtered to exclude non-significant values (absolute value > 0).
#' @examples
#' # This example assumes you have SiMLR result `simlrResult`, matrices `simlrMats`, and you want to 
#' # interpret the first principal component "PC1".
#' # simlrResult <- list(v = list(PC = matrix(runif(20), ncol = 2)))
#' # simlrMats <- list(PC = matrix(runif(100), ncol = 10))
#' # simlrVariable <- "PC1"
#' # interpretedVector <- interpret_simlr_vector(simlrResult, simlrMats, simlrVariable)
#' # print(interpretedVector)
#' @importFrom stringr str_match str_extract
#' @export
interpret_simlr_vector <- function( simlrResult, simlrMats, simlrVariable, n2show = 5, shortnames=TRUE, return_dataframe=FALSE ) {

  split_string_correctly <- function(input_string) {
    # Extract the leading alphabetic characters (possibly including numbers within the alphabetic segment)
    alpha_part <- str_match(input_string, "([A-Za-z0-9]+(?=[A-Za-z]+[0-9]+$))[A-Za-z]*")[,1]
    
    # Extract the numeric part at the end
    numeric_part <- str_extract(input_string, "[0-9]+$")
    
    c( alpha_part, numeric_part)
  }
  varparts = split_string_correctly( simlrVariable )
  varparts[1]=gsub("PC","",varparts[1])
  nmslist=list()
  if ( shortnames ) {
    for ( k in 1:length(simlrMats) ) 
      nmslist[[names(simlrMats)[k]]]=shorten_pymm_names(colnames(simlrMats[[k]]))
  } else {
    for ( k in 1:length(simlrMats) ) 
      nmslist[[names(simlrMats)[k]]]=colnames(simlrMats[[k]])
  }

  # Extract the vector for the given modality and region, and normalize it
  t1vec <- abs(simlrResult$v[[varparts[1]]][, as.integer(varparts[2])])
  t1vec=t1vec/max(t1vec)
  
  # Assign names to the vector elements from the names list
  names(t1vec) <- nmslist[[varparts[1]]]
  
  # Sort the vector in decreasing order and select the top 'n2show' elements
  # If 'n2show' is NULL or greater than the length of t1vec, use the length of t1vec
  n_items_to_show <- if (is.null(n2show)) length(t1vec) else min(c(n2show, length(t1vec)))
  t1vec_sorted <- head(t1vec[order(t1vec, decreasing = TRUE)], n_items_to_show)
  
  # Filter out non-significant values (absolute value > 0)
  t1vec_filtered <- t1vec_sorted[abs(t1vec_sorted) > 0]
  if ( return_dataframe ) {
    t1vec_filtered=data.frame( anat=names(t1vec_filtered), values=t1vec_filtered)
  }
  return(t1vec_filtered)
}



#' Interpret SiMLR Vector
#'
#' This function interprets a vector from SiMLR (similarity-driven multivariate linear reconstruction)
#' specifically focusing on a given variable (e.g., a specific principal component or cluster). It extracts and normalizes the vector associated 
#' with the specified SiMLR variable, sorts it to identify the top elements, and optionally filters out non-significant values. 
#' This function is useful for understanding the contribution of different features in the context of the SiMLR analysis. Assumes this input is generated by antspymm_simlr.
#'
#' @param simlrResult a specific v matrix out of SiMLR
#' @param simlrVariable A string specifying the variable within `simlrResult` to interpret. The variable 
#' name should include both an identifier (e.g., "PC" for principal component) and a numeric index.
#' @param n2show An integer specifying the number of top elements to show from the sorted, normalized vector. 
#' Defaults to 5. If `NULL` or greater than the length of the vector, all elements are shown.
#' @param shortnames boolean
#' @param return_dataframe boolean
#' @return A named vector of the top `n2show` elements (or all if `n2show` is `NULL` or too large), 
#' sorted in decreasing order of their absolute values. Elements are named according to their identifiers 
#' in `simlrMats` and filtered to exclude non-significant values (absolute value > 0).
#' @examples
#' # This example assumes you have SiMLR result `simlrResult`, matrices `simlrMats`, and you want to 
#' # interpret the first principal component "PC1".
#' # simlrResult <- list(v = list(PC = matrix(runif(20), ncol = 2)))
#' # simlrMats <- list(PC = matrix(runif(100), ncol = 10))
#' # simlrVariable <- "PC1"
#' # interpretedVector <- interpret_simlr_vector2(simlrResult$v[[1]] )
#' # print(interpretedVector)
#' @importFrom stringr str_match str_extract
#' @export
interpret_simlr_vector2 <- function( simlrResult, simlrVariable, n2show = 5, shortnames=TRUE, return_dataframe=FALSE ) {

  split_string_correctly <- function(input_string) {
    # Extract the leading alphabetic characters (possibly including numbers within the alphabetic segment)
    alpha_part <- str_match(input_string, "([A-Za-z0-9]+(?=[A-Za-z]+[0-9]+$))[A-Za-z]*")[,1]
    
    # Extract the numeric part at the end
    numeric_part <- str_extract(input_string, "[0-9]+$")
    
    c( alpha_part, numeric_part)
  }
  varparts = split_string_correctly( simlrVariable )
  varparts[1]=gsub("PC","",varparts[1])
  if ( shortnames ) {
    nmslist=shorten_pymm_names( rownames( simlrResult ) )
  } else {
    nmslist = rownames( simlrResult )
  }

  # Extract the vector for the given modality and region, and normalize it
  t1vec <- abs(simlrResult[, simlrVariable ])
  t1vec=t1vec/max(t1vec)
  
  # Assign names to the vector elements from the names list
  names(t1vec) <- nmslist
  
  # Sort the vector in decreasing order and select the top 'n2show' elements
  # If 'n2show' is NULL or greater than the length of t1vec, use the length of t1vec
  n_items_to_show <- if (is.null(n2show)) length(t1vec) else min(c(n2show, length(t1vec)))
  t1vec_sorted <- head(t1vec[order(t1vec, decreasing = TRUE)], n_items_to_show)
  
  # Filter out non-significant values (absolute value > 0)
  t1vec_filtered <- t1vec_sorted[abs(t1vec_sorted) > 0]
  if ( return_dataframe ) {
    t1vec_filtered=data.frame( anat=names(t1vec_filtered), values=t1vec_filtered)
  }
  return(t1vec_filtered)
}


#' Linear Mixed Effects Model Analysis with ANOVA and Effect Size Calculation
#'
#' This function fits two linear mixed effects models: a base model without the main predictor of interest
#' and a full model including the main predictor. It then compares these models using ANOVA to test the
#' significance of adding the main predictor. Additionally, it calculates the effect sizes for the predictor
#' in the full model. This function is designed to facilitate the analysis of data where both fixed and
#' random effects are present, accommodating complex experimental designs.
#' NOTE: this function will scale variables internally to aid coefficient estimate visualization.
#'
#' @param data A data frame containing the variables referenced in the model formulas.
#' @param outcome The name of the dependent variable (outcome) as a string.
#' @param predictor The name of the main predictor variable as a string.
#' @param fixed_effects A string specifying the fixed effects to be included in the model, excluding the main predictor.
#' @param random_effects A string specifying the random effects to be included in the model. we assume that a subject ID is the first entry if this is a vector.
#' @param predictoroperator either a \code{+} or \code{*}
#' @param verbose boolean
#' @return A list containing the fitted full model object, ANOVA model comparison, calculated effect sizes for the
#' predictor, coefficients of the full model, and the count of unique levels in the random effects variable.
#' @examples
#' # Assuming 'data' is your dataset with columns 'outcome', 'predictor', 'fixed_var1', ...,
#' # and 'subject' as the random effect:
#' # results <- lmer_anv_p_and_d(data, "outcome", "predictor", "fixed_var1 + fixed_var2", "subject")
#' # summary(results$full_model) # Full model summary
#' # results$model_comparison # ANOVA comparison
#' # results$effect_sizes # Effect sizes
#' @importFrom lmerTest lmer
#' @importFrom stats anova
#' @importFrom stringr str_extract str_match
#' @importFrom effectsize t_to_d
#' @export
lmer_anv_p_and_d <- function(data, outcome, predictor, fixed_effects, random_effects, predictoroperator='*', verbose=FALSE ) {
  # Validate input to ensure variables exist in the data frame
  stopifnot(is.character(outcome), is.character(predictor), is.character(fixed_effects), is.character(random_effects))

  hasConverged <- function (mm) {
#    if ( !inherits(mm, "merMod")) stop("Error: must pass a lmerMod object")
    retval <- NULL    
    if(is.null(unlist(mm@optinfo$conv$lme4))) {
      retval = 1
    }
    else {
      if (isSingular(mm)) {
        retval = 0
      } else {
        retval = -1
      }
    }
    return(retval)
  }
  # Construct the model formulas directly
  base_model_formula <- paste( outcome, "~", convert_to_random_effects(random_effects), "+", fixed_effects )
  if ( verbose ) {
    print("base_model_formula")
    print( base_model_formula )
    }
  full_model_formula <- paste( base_model_formula, predictoroperator, predictor )
  base_model_formula = as.formula( base_model_formula )
  full_model_formula = as.formula( full_model_formula )
  all_vars = all.vars( full_model_formula )
  missing_vars <- setdiff(all_vars, colnames(data))
  if ( verbose ) {
    print( full_model_formula )
    }
  if (length(missing_vars) > 0) {
    stop("The following variables are missing in the data: ", paste(missing_vars, collapse = ", "), ".")
  }
  
  
  # Subset data to exclude rows with NA values for relevant variables
  datasub <- na.omit(data[, all_vars])
  datasub = scale_variables_in_equation( datasub, full_model_formula )

  # Fit the linear mixed models using lmer from the lme4 package
  base_model = tryCatch({
      lmer(base_model_formula, data = datasub, REML = FALSE)
    }, error = function(e) {
      NULL
    })
  if ( is.null(base_model) ) return(NULL)
  full_model = tryCatch({
    lmer(full_model_formula, data = datasub, REML = FALSE)
    }, error = function(e) {
      NULL
    })
  if ( is.null(full_model)  ) return(NULL)
  if ( hasConverged(full_model) != 1 ) return(NULL)
  
  # Perform ANOVA to compare the models
  model_comparison <- anova(base_model, full_model)
  
  # Calculate effect sizes for the full model
  coefs <- summary(full_model)$coefficients
  ndf <- length(unique(datasub[[random_effects[1]]])) # Now using datasub for N calculation
  effect_sizes <- effectsize::t_to_d(coefs[, "t value"], rep(ndf, nrow(coefs)))
  effect_sizes <- data.frame(effect_sizes)
  rownames(effect_sizes) <- rownames(coefs)
  effect_sizes <- effect_sizes[grep(predictor, rownames(effect_sizes)), ]

  # Prepare and return the results
  results <- list(
    full_model = full_model,
    model_comparison = model_comparison,
    effect_sizes = effect_sizes,
    coefficients = coefs,
    n = length(unique(datasub[[random_effects[1]]]))
  )
  
  return(results)
}



#' Set Seed Based on Current Time
#'
#' This function sets the random number generator seed based on the current time,
#' with a fine resolution of seconds. This ensures a different seed is used each time
#' the function is called, provided calls are at least one second apart.
#'
#' @return The numeric value used as the seed, derived from the current time in seconds.
#' @examples
#' seedValue <- setSeedBasedOnTime()
#' print(seedValue)
#' @export
setSeedBasedOnTime <- function() {
  op <- options(digits.secs = 8)
  # Get the current time
  currentTime <- Sys.time()
  
  # Convert the current time to a numeric value
  # numericTime <- as.numeric(currentTime, units = "secs")
  numericTime = as.integer(substr(as.character(Sys.time()),22,200))
  # Use the numeric time as the seed
  set.seed(numericTime)
  
  # Optionally, return the seed value used
  return(numericTime)
}



#' Plot Categorical Data
#'
#' This function creates a visualization for categorical data using ggplot2. It takes a dataset and a column name,
#' generates a frequency table for the categorical data, and plots the frequencies.
#'
#' @param datatoplot A data frame containing the dataset to be plotted.
#' @param columnName The name of the column in the dataset that contains categorical data.
#'
#' @return A ggplot object representing the frequency of each category in the specified column.
#' @export
#'
#' @examples
#' # plotCategoricalData(myData, "myCategoryColumn")
plotCategoricalData <- function(  datatoplot, columnName ) {
  # Convert the frequency table to a data frame for ggplot

  freqDataFrame <- as.data.frame(datatoplot[[columnName]]$FrequencyTable)
  names(freqDataFrame) <- c("Category", "Frequency")
  
  # Add a column to indicate the subject's category
  freqDataFrame$SubjectCategory <- freqDataFrame$Category == datatoplot[[columnName]]$SubjectCategory
  
  # Plot
  ggplot(freqDataFrame, aes(x = Category, y = Frequency, fill = SubjectCategory)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "lightblue")) +
    theme_minimal() +
    labs(title = paste("Dist.", columnName, "/ Subject (red)"),
         y = "Frequency",
         x = columnName) +
    coord_flip()  # Flip for horizontal bars
}


#' Convert Character Vector to Random Effects Terms for lme4
#'
#' This function takes a vector of strings representing grouping factors and converts 
#' each string into a random effect term formatted for use in lme4 model formulas. 
#' Specifically, it formats them as random intercepts for the specified grouping factors.
#'
#' @param variables A character vector where each element is the name of a grouping factor
#' for which a random intercept will be included in the model.
#'
#' @return A character vector of random effect terms, each formatted as `(1|group)`, 
#' where `group` is replaced with each string from the input vector.
#'
#' @examples
#' grouping_factors <- c("school", "class")
#' random_effects <- convert_to_random_effects(grouping_factors)
#' print(random_effects)
#'
#' @export
convert_to_random_effects <- function(variables) {
  # Check if the input is indeed a vector of strings
  if (!is.character(variables)) {
    stop("The input should be a character vector.")
  }
  
  # Create the random effect terms
  random_effects <- paste0("(1|", variables, ")")
  
  return( paste(random_effects,collapse="+"))
}


#' Normative Summary
#'
#' This function provides a normative summary for a given subject within a dataset. It can be configured
#' to focus on specific columns and can be zoomed into particular data points.
#'
#' @param data A data frame containing the dataset for analysis.
#' @param subjectRow The row number or identifier for the subject of interest within the dataset.
#' @param columns A vector of column names in the dataset that should be summarized.
#' @param zoom A parameter to specify the focus or zoom level of the summary.
#' @param idcolumn column name for the unique subject ID
#' @param verbose A logical flag indicating whether to print detailed output.
#'
#' @return A summary object containing normative statistics for the specified columns of the subject.
#' @export
#'
#' @examples
#' # normativeSummary(myData, 1, c("Column1", "Column2"), zoom = 1, verbose = TRUE)
normativeSummary <- function(data, subjectRow, columns, zoom, idcolumn='commonID', verbose=TRUE) {
  if(!is.data.frame(data)) stop("The 'data' input must be a data frame.")
  if(!all(columns %in% names(data))) stop("All specified columns must exist in the data frame.")
  if(subjectRow > nrow(data) || subjectRow < 1) stop("Subject row is out of bounds.")
  
  if ( ! missing( zoom ) ) {
    dataz=find_closest_subjects( data[subjectRow,], data, k=zoom, 'commonSex', 'commonAge')
    data = do.call(rbind, dataz)
    subjectRow=1
  }
  succcolor='deepskyblue4'#  'dodgerblue1'
  summaryList <- list()
  histList = list()
  
  for (col in columns) {
    columnData <- data[[col]]
    if ( subtyper::fs(antspymm_vartype(col) %in% c("T1","T2Flair","DTI")) & 'brainVolume' %in% colnames(data)) {
      columnData=columnData/data$brainVolume
      if ( verbose ) {
        print(paste("normalize",col,'by BV'))
      }
    }
    isNumeric <- is.numeric(columnData)
    if (isNumeric) {
      # Process numeric data
      meanVal <- mean(columnData, na.rm = TRUE)
      sdVal <- (stats::sd(columnData, na.rm = TRUE))
      subjectScore <- columnData[subjectRow]
      zScore <- (subjectScore - meanVal) / sdVal
      if ( verbose ) {
        print( paste( col, meanVal, sdVal, subjectScore, zScore ) )
      }
      if ( !is.na( subjectScore) ) {
        tTestResult <- t.test(columnData, alternative = "greater", mu = subjectScore)
      } else tTestResult=list(estimate=NA,p.value=NA,statistic=NA)

      ttl=paste(shorten_pymm_names(col),'sub. (blue) vs pop.')
      histdf=data.frame( vid=rep("pop",length(columnData)), value=columnData )
      if ( max( histdf$value,na.rm=T ) < 1 ) {
        scl=100/(range( histdf$value,na.rm=T )[2]-range( histdf$value,na.rm=T  )[1])
        histdf$value=histdf$value * scl
        ttl=paste(ttl,"\n : vals scaled by",insight::format_value(scl))
      }
      histdf[subjectRow,'vid']='subject'

      histList[[col]]=gghistogram(histdf, x = 'value',  
        add.params=list(size=1.25,linetype = "dashed"),
        add = "mean", add_density = FALSE, title=ttl, fill='vid', legend='none')

      summaryList[[col]] <- list(
        Mean = meanVal,
        SD = sdVal,
        SubjectScore = subjectScore,
        ZScore = zScore,
        TTestPValue = tTestResult$p.value
      )
    } else {
      # Process categorical data
      freqTable <- table(columnData)
      subjectCategory <- as.character(columnData[subjectRow])
      
      summaryList[[col]] <- list(
        FrequencyTable = freqTable,
        SubjectCategory = subjectCategory
      )
    }
  }
  
  # Assuming 'summaryList' contains the processed data
  for (col in columns) {
    if (!is.numeric(data[[col]])) {
      # Assume 'summaryList' is available from the earlier processing
      freqTable <- summaryList[[col]]$FrequencyTable
      subjectCategory <- summaryList[[col]]$SubjectCategory
      histList[[col]]=plotCategoricalData( summaryList, col)
    }
  }

  # Print summary

  # Visualization: For simplicity, focusing on numeric columns for z-score plot
  numericColumns <- sapply(summaryList, function(x) "ZScore" %in% names(x))
  if (any(numericColumns)) {
    zScores <- sapply(summaryList[numericColumns], function(x) x$ZScore)
    names(zScores) <- names(summaryList)[numericColumns]
    zScoreDataFrame <- data.frame(Column = names(zScores), ZScore = zScores)
    zScoreDataFrame$Column = shorten_pymm_names(zScoreDataFrame$Column )
    if ( verbose )
      print( zScoreDataFrame )
    histList[[paste0(col,".z")]]=ggplot(zScoreDataFrame, aes(x = Column, y = ZScore, fill = ZScore)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_gradient2(low = succcolor, high = "red", mid = "white", midpoint = 0) +
      theme_minimal() +
      labs(title = "Z-Scores vs. pop.",
           y = "Z-Score",
           x = "") +
      coord_flip() + theme(legend.key.width = unit(0.005, "npc"))
  } else {
    cat("No numeric data for z-score plot.\n")
  }

  grid.arrange(grobs=histList,ncol=round(sqrt(length(histList))),top=paste("Normative Results:",data[subjectRow,idcolumn]))
  return(summaryList)
}


#' Find Closest Subjects
#'
#' This function identifies the closest subjects in a target dataset to a reference dataset
#' based on specified criteria such as sex and age.
#'
#' @param reference_df A data frame containing the reference dataset.
#' @param target_df A data frame containing the target dataset to search within.
#' @param k The number of closest subjects to find.
#' @param sex_col_name The name of the column in the datasets that contains subjects' sex.
#' @param age_col_name The name of the column in the datasets that contains subjects' age.
#'
#' @return A data frame with the rows from `target_df` that are closest to `reference_df` based on the criteria.
#' @export
#'
#' @examples
#' # find_closest_subjects(referenceData, targetData, 5, "Gender", "Age")
find_closest_subjects <- function(reference_df, target_df, k, sex_col_name, age_col_name) {
  # Error checking for column names
  if(!(sex_col_name %in% names(reference_df)) | !(age_col_name %in% names(reference_df))) {
    stop("Reference data frame must contain the specified 'sex' and 'age' column names.")
  }
  if(!(sex_col_name %in% names(target_df)) | !(age_col_name %in% names(target_df))) {
    stop("Target data frame must contain the specified 'sex' and 'age' column names.")
  }
  if(!is.numeric(k) | k <= 0 | length(k) != 1) {
    stop("'k' must be a positive integer.")
  }
  
  # Initialize an empty list to store the results
  results <- list()
  
  # Loop through each row in the reference data frame
  for (i in 1:nrow(reference_df)) {
    # Extract the current row's sex and age
    current_sex <- reference_df[[sex_col_name]][i]
    current_age <- reference_df[[age_col_name]][i]
    
    # Filter the target data frame for subjects of the same sex
    same_sex_df <- target_df[target_df[[sex_col_name]] == current_sex,]
    
    # Calculate the absolute difference in age between the reference subject and all target subjects
    same_sex_df$age_diff <- abs(same_sex_df[[age_col_name]] - current_age)
    
    # Sort the data frame by age difference
    sorted_subjects <- same_sex_df[order(same_sex_df$age_diff), ]
    
    # Select the top 'k' closest subjects, handling cases where there are fewer than 'k' subjects
    num_rows_to_select <- min(nrow(sorted_subjects), k)
    closest_subjects <- sorted_subjects[1:num_rows_to_select, ]
    
    # Store the result
    results[[i]] <- closest_subjects
  }
  
  # Return the final result
  return(results)
}




#' Replace Values in a Vector
#'
#' This function replaces specified values within a vector with new values.
#'
#' @param vec The original vector in which values need to be replaced.
#' @param old_values A vector containing the old values to be replaced.
#' @param new_values A vector containing the new values to replace the old ones.
#'
#' @return A vector with the old values replaced by the new values.
#' @export
#'
#' @examples
#' replace_values(c(1, 2, 3, 4), c(2, 4), c(20, 40))
replace_values <- function(vec, old_values, new_values) {
  # Ensure old_values and new_values are vectors of the same length
  if (length(old_values) != length(new_values)) {
    stop("old_values and new_values must have the same length")
  }
  
  # Check for NA in old_values to handle it specifically
  for (i in seq_along(old_values)) {
    if (is.na(old_values[i])) {
      vec[is.na(vec)] <- new_values[i]
    } else {
      vec[vec == old_values[i]] <- new_values[i]
    }
  }
  
  return(vec)
}





#' Perform SiMLR Analysis on Multimodal ANTsPyMM Data
#'
#' This function processes multimodal data using SiMLR. It is designed
#' to be flexible, allowing for various preprocessing steps and analysis options. The analysis
#' can be adjusted through multiple parameters, offering control over the inclusion of certain
#' data types, permutation testing, and more.
#'
#' @param blaster A dataframe containing multimodal data for analysis.
#' @param select_training_boolean boolean vector to define which entries are in training data
#' @param connect_cog Vector of column names to be treated as a special target matrix;  often used for cognitive data and in a superivsed variant of simlr.  Exclude this argument if this is unclear.
#' @param energy The type of energy model to use for similarity analysis. Defaults to 'reg'.
#' @param nsimlr Number of similarity analyses to perform. Defaults to 5.
#' @param covariates any covariates to adjust training matrices. if covariates is set to 'mean' then the rowwise mean will be factored out of each matrix.
#' @param myseed Seed for random number generation to ensure reproducibility. Defaults to 3.
#' @param doAsym integer 0 for FALSE, 1 for TRUE and 2 for separate matrices for asymm variables.
#' @param returnidps Logical indicating whether to return the intermediate processing steps' results. Defaults to FALSE.
#' @param restrictDFN Logical indicating whether to restrict analysis to default network features. Defaults to FALSE.
#' @param resnetGradeThresh image quality threshold (higher better).
#' @param doperm Logical indicating whether to perform permutation tests. Defaults to FALSE.  Will randomize image features in the training data and thus leads to "randomized" but still regularized projections.
#' @param exclusions vector of strings to exclude from predictors
#' @param inclusions vector of strings to include in predictors
#' @param sparseness vector or scalar value to set sparseness
#' @param iterations int value to set max iterations
#' @param verbose boolean
#' @return A list containing the results of the similarity analysis and related data.
#' @export
#' @examples
#' # Example usage:
#' # result <- antspymm_simlr(dataframe)
antspymm_simlr = function( blaster, select_training_boolean, connect_cog,  energy=c('cca','reg','lrr'), nsimlr=5, covariates='1', myseed=3,  doAsym=TRUE, returnidps=FALSE, restrictDFN=FALSE, 
resnetGradeThresh=1.02, doperm=FALSE, 
exclusions=NULL, inclusions=NULL, sparseness=NULL, iterations=NULL, verbose=FALSE ) 
{
  safegrep <- function(pattern, x, ...) {
    result <- grep(pattern, x, ...)
    if (length(result) == 0) {
      return(1:length(x))
    }
    return(result)
  }
  safeclean = function( pattern, x,fixed=FALSE, exclude=TRUE) {
    mysub=grep(pattern,x,fixed=fixed)
    if ( length(mysub) == 0 ) return( x )
    if ( exclude ) return( x[-mysub] ) else return( x[mysub] )
  }
  idps=antspymm_predictors(blaster,TRUE,TRUE)
  rsfnames = idps[ grepl("rsfMRI",idps) ]
  if ( length(rsfnames) > 0 ) rsfnames = rsfnames[ safegrep("_2_",rsfnames)]
  if ( !all(grepl("rsfMRI", rsfnames )) ) rsfnames=c()
  if ( !is.null(exclusions)) {
    for ( x in exclusions ) {
      idps=safeclean(x,idps)
      rsfnames=safeclean(x,rsfnames)
    }
  }
  if ( !is.null(inclusions)) {
    for ( x in inclusions ) {
      idps=safeclean(x,idps,exclude=FALSE)
      rsfnames=safeclean(x,rsfnames,exclude=FALSE)
    }
  }
  idps=idps[ -multigrep(antspymm_nuisance_names()[-3],idps)]
  if ( doAsym == 0 ) {
    idps=safeclean("Asym",idps)
    } else {
    idps=safeclean("Asymcit168",idps)
  }
  idps=safeclean("cleanup",idps)
  idps=safeclean("snseg",idps)
  idps=safeclean("_deep_",idps)
  idps=safeclean("peraf",idps)
  idps=safeclean("alff",idps)
  idps=safeclean("LRAVGcit168",idps)
  idps=safeclean("_l_",idps,fixed=TRUE)
  idps=safeclean("_r_",idps,fixed=TRUE)
  if ( restrictDFN ) {
    rsfnames = rsfnames[ safegrep("Default",rsfnames)]
  } else {
#    rsfnames = rsfnames[ multigrep( c("imbic","TempPar"),rsfnames)]
  }
  perfnames = idps[ multigrep( c("perf_cbf_mean_"),idps,intersect=TRUE)]
  t1names = idps[ multigrep( c("T1Hier"),idps,intersect=TRUE)]
  dtnames = unique( c( 
    idps[ multigrep( c("mean_fa","DTI"),idps,intersect=TRUE)],
    idps[ multigrep( c("mean_md","DTI"),idps,intersect=TRUE)] ))

  t1asymnames=c()
  dtasymnames=c()
  pfasymnames=c()
  if ( doAsym == 2 ) {
    t1nms = t1names
    t1asymnames = t1nms[ grep("Asym",t1nms)]
    t1names = t1nms[ !( t1nms %in%  t1asymnames ) ]

    dtnms = dtnames
    dtasymnames = dtnms[ grep("Asym",dtnms)]
    dtnames = dtnames[ !( dtnames %in%  dtasymnames ) ]

    pfnms = perfnames
    pfasymnames = pfnms[ grep("Asym",pfnms)]
    perfnames = pfnms[ !( pfnms %in%  pfasymnames ) ]
    }

  idps=unique(t1names)
  idplist = list()
  idplist[["t1"]]=t1names
  if ( length(dtnames) > 0 ) {
    idps = c( idps, unique(dtnames) )
    idplist[["dt"]]=dtnames
  }

  if ( length(rsfnames) > 0 ) {
    idps = c( idps, unique(rsfnames) )
    idplist[["rsf"]]=rsfnames
  }
  if ( length(perfnames) > 0 ) {
    idps = c( idps, unique(perfnames) )
    idplist[["perf"]]=perfnames
  }
  if ( length(t1asymnames) > 0 ) {
    idps = c( idps, unique(t1asymnames) )
    idplist[["t1a"]]=t1asymnames
  }

  if ( length(dtasymnames) > 0 ) {
    idps = c( idps, unique(dtasymnames) )
    idplist[["dta"]]=dtasymnames
  }

  if ( length(pfasymnames) > 0 ) {
    idps = c( idps, unique(pfasymnames) )
    idplist[["pfa"]]=pfasymnames
  }
  if ( !missing( connect_cog ) ) { 
    idplist[["cg"]]=connect_cog
  }
  if ( verbose ) {
    print(names(idplist))
    print(sample(idps,10))
  }
  if ( returnidps ) return(idps)
  allnna=select_training_boolean[  blaster$T1Hier_resnetGrade >= resnetGradeThresh ]
  blaster2=blaster[  blaster$T1Hier_resnetGrade >= resnetGradeThresh, ]
  stopifnot( min(dim(blaster2)) > 3 )
  if ( verbose ) {
    print("dim( subsetdataframe)")
    print(dim(blaster2) )
  }
  #################################################
  nperms=0
  matsFull = list()
  mats = list()
  for ( kk in 1:length(idplist)) {
      matsFull[[ names(idplist)[kk] ]] = blaster[,idplist[[kk]]]
      mats[[ names(idplist)[kk] ]] = antsrimpute( blaster2[allnna,idplist[[kk]]] )
      }
  if ( verbose ) print("mats done")
  if ( doperm ) {
    nada=setSeedBasedOnTime()
    sss=sample( 1:nrow( matsFull[[1]]  ))
    for ( jj in 1:length( mats ) ) {
        ss=sample( 1:nrow( mats[[jj]]  ))
        mats[[jj]]=mats[[jj]][sample( 1:nrow( mats[[jj]]  )),]
    }
  }
  nms = names(mats)
  regs0 = list()

  update_residuals <- function(mats, x, covariate, blaster2, allnna) {
    if ( is.null(covariate) ) return(mats[[x]])
    if ( covariate == 'robust' ) {
      return( robustMatrixTransform( data.matrix( mats[[x]] ) ) )
    }
    if ( covariate == 'mean' ) {
      mymean=rowMeans(  data.matrix( mats[[x]] ) )
      covariate2='mymean'
    } else covariate2=covariate
    formula <- as.formula(paste("data.matrix(mats[[", x, "]]) ~ ", covariate2))
    fit <- lm(formula, data = blaster2[allnna, ])
    residuals(fit)
    }

  if ( verbose) print("setting up regularization")
  for ( x in 1:length(mats)) {
      if ( verbose ) {
        if ( x == 1 ) print(paste("training n= ",nrow(mats[[x]])))
        cat(paste0(names(mats)[x],"..."))
      }
      mats[[x]]=update_residuals( mats, x, covariates, blaster2, allnna )
      mats[[x]]=data.matrix(mats[[x]])
      mycor = cor( mats[[x]] )
      mycor[mycor < 0.8]=0
      regs0[[x]]=data.matrix(mycor)
      }
  names(regs0)=names(mats)
  regs = regs0 # regularizeSimlr( mats, fraction=0.05, sigma=rep(2.0,length(mats)) )
  if ( verbose ) print("regularizeSimlr done")
  names(regs0)=names(mats)
  names(regs)=names(mats)
  if ( !missing( connect_cog ) ) {
    regs[["cg"]] = Matrix(regs0[["cg"]], sparse = TRUE) 
    print("regularize cg")
  }

#  if ( !doperm )
#    for ( pp in 1:length(regs)) plot(image(regs[[pp]]))

  if ( verbose ) print("loop mat")
  for ( k in 1:length(mats)) {
    if ( ncol(mats[[k]]) != ncol(regs[[k]]) ) {
      regs[[k]]=Matrix(regs0[[k]], sparse = TRUE) 
      msg=paste("regularization cols not equal",k,ncol(mats[[k]]),ncol(regs[[k]]),names(mats)[k])
      message(msg)
      # stop( )
      }
    }
  if ( verbose ) print("loopmatdone")
  ########### zzz ############
  myjr = T
  prescaling = c( 'center', 'np' )
  optimus = 'lineSearch'
  maxits = 1000
  if ( ! is.null( iterations ) ) maxits = iterations
  if ( verbose ) print( paste( "maxits",maxits) )
  ebber = 0.99
  pizzer = rep( "positive", length(mats) )
  objectiver='cca';mixer = 'pca'
  if ( energy == 'reg') {
    objectiver='regression';mixer = 'ica'
  }
  if ( energy == 'lrr') {
    objectiver='lowRankRegression';mixer = 'pca'
  }
  if ( verbose ) print("sparseness begin")
  sparval = rep( 0.8, length( mats ))
  if ( ! is.null( sparseness ) ) {
    if ( length( sparseness ) == length(mats) ) {
      sparval = sparseness
    } else sparval = rep( sparseness[1], length( mats ))
    if ( verbose ) {
      print('sparseness')
      print(sparseness)
    }
  }
  
  if ( nsimlr < 1 ) {
    ctit=0
    for ( jj in 1:length(mats) ) {
      ctit=ctit+ncol(mats[[jj]])
      sparval[jj] = 1.0 - 20/ncol(mats[[jj]])
      if ( sparval[jj] < 0 ) sparval[jj] = 0.5
    }
    nsimlr = round( ctit * nsimlr )
    message(paste("nsimlr",nsimlr))
#    print(paste("nsimlr",nsimlr))
#    print(sparval)
  }

  if ( verbose ) {
    print("initu begin")
  }
  initu = initializeSimlr(
      mats,
      nsimlr,
      jointReduction = myjr,
      zeroUpper = FALSE,
      uAlgorithm = "pca",
      addNoise = 0 )
  if ( verbose ) print("initu done")

  # initu = initu[,(ncol(initu)-nsimlr):ncol(initu)]
  # initu = initu[,1:nsimlr]

  if ( ! missing( connect_cog ) ) {
    clist = list()
    inflammNums=which(names(mats)=='cg')
    for ( j in 1:length( mats ) ) clist[[j]] = inflammNums
    for ( j in inflammNums )
      clist[[j]] = (1:length(mats))[ -inflammNums ]
    } else clist=NULL
    

  simlrX = simlr( mats, regs, 
    iterations=maxits, 
    verbose= !doperm,
    randomSeed = myseed,
    mixAlg=mixer,
    energyType=objectiver,
    scale = prescaling,
    sparsenessQuantiles=sparval,
    expBeta = ebber,
    positivities = pizzer, 
    connectors=clist,
    optimizationStyle=optimus,
    initialUMatrix=initu )
  if ( verbose ) print('simlr done')
  for ( kk in 1:length(mats) ) {
    rownames(simlrX$v[[kk]])=idplist[[kk]]
    temp = simlrX$v[[kk]]
    if ( pizzer[kk] == 'positive' ) {
#      for ( n in 1:ncol(temp)) temp[,n]=abs(temp[,n])/max(abs(temp[,n]))
#      simlrX$v[[kk]]=eliminateNonUniqueColumns(temp)
      }
    }

  #################
  nsimx=nsimlr
  nms = names( simlrX$v ) = names(mats)
  simmat = data.matrix(matsFull[[1]] )%*% abs( simlrX$v[[1]] )
  colnames( simmat ) = paste0(nms[1],colnames( simmat ))
  for ( j in 2:length(mats)) {
      if (names(mats)[j]=='cg' & pizzer[j] != 'positive' ) {
        temp = data.matrix(matsFull[[j]] ) %*% ( simlrX$v[[j]])
      } else temp = data.matrix(matsFull[[j]] ) %*% abs( simlrX$v[[j]] )
      colnames( temp ) = paste0(nms[j],colnames( temp ))
      simmat = cbind( simmat, temp )
  }
  blaster2sim = cbind( blaster, simmat )
  nsim = ncol( simlrX$v[[1]] )
  simnames = colnames(simmat)
  kk=1
  nmats=1:length(matsFull)
  matsB=mats
  for ( kk in 1:length(mats)) matsB[[kk]]=data.matrix(matsB[[kk]])
  kk=length(mats)
  temp = predictSimlr( matsB, simlrX, targetMatrix=kk, 
        sourceMatrices=nmats[nmats!=kk] )

  return( list( demog=blaster2sim, mats=matsFull, simnames=simnames, simlrX=simlrX, energy=energy, temp=temp ) )
  ################
  }


#' Match Two Data Frames Based on Nearest Neighbor Matching
#'
#' This function matches two data frames on specified variables using the nearest neighbor
#' matching method and checks if the key variable is statistically significantly different 
#' between the matched data frames.
#'
#' @param df1 A data frame.
#' @param df2 A data frame.
#' @param match_vars A character vector of variable names to match on.
#' @return A list containing the matched data frames and the result of the t-test.
#' @import proxy
#' @export
#' @examples
#' set.seed(123)
#' df1 <- data.frame(
#'   id = 1:100,
#'   age = rnorm(100, mean = 30, sd = 5),
#'   gender = sample(c("Male", "Female"), 100, replace = TRUE),
#'   score = rnorm(100, mean = 75, sd = 10)
#' )
#'
#' df2 <- data.frame(
#'   id = 101:200,
#'   age = rnorm(100, mean = 30, sd = 5),
#'   gender = sample(c("Male", "Female"), 100, replace = TRUE),
#'   score = rnorm(100, mean = 70, sd = 15)
#' )
#'
#' result <- match_data_frames(df1, df2, match_vars = c("age", "gender") )
match_data_frames <- function(df1, df2, match_vars) {
  ocolnames = intersect( colnames(df1), colnames(df2))
  # Convert categorical variables to numeric
  for (var in match_vars) {
    if (is.factor(df1[[var]]) || is.character(df1[[var]])) {
      levels <- unique(c(df1[[var]], df2[[var]]))
      df1[[paste0(var, "_numeric")]] <- as.numeric(factor(df1[[var]], levels = levels))
      df2[[paste0(var, "_numeric")]] <- as.numeric(factor(df2[[var]], levels = levels))
    } else {
      df1[[paste0(var, "_numeric")]] <- df1[[var]]
      df2[[paste0(var, "_numeric")]] <- df2[[var]]
    }
  }

  # Normalize the numeric variables
  normalize <- function(x) {
    return((x - min(x,na.rm=T)) / (max(x,na.rm=T) - min(x,na.rm=T)))
  }
  normalize <- function(x) {
    eps=0.001
    myq = quantile(x,c(1.0-eps,eps), na.rm=T)
    return((x - min(x,na.rm=T)) / (myq[1]-myq[2]))
  }
  
  for (var in match_vars) {
    norm_var <- paste0(var, "_numeric")
    df1[[paste0(norm_var, "_normalized")]] <- normalize(df1[[norm_var]])
    df2[[paste0(norm_var, "_normalized")]] <- normalize(df2[[norm_var]])
  }
  
  
  # Calculate distances and find the nearest neighbor
  distance_vars <- paste0(match_vars, "_numeric_normalized")
  distances <- proxy::dist(df1[, distance_vars], df2[, distance_vars], method = "Euclidean")
  nearest_neighbors <- apply(as.matrix(distances), 1, which.min)
  
  # Create matched data frames based on nearest neighbors
  matched_df2 <- df2[nearest_neighbors, ]
  return(matched_df2[ ,ocolnames] )
  }


#' Write a list of data frames to disk with specific naming convention
#'
#' This function writes each data frame in a list to a separate CSV file on disk,
#' using the names of each data frame to create unique filenames.
#'
#' @param data_list A list of data frames to write to disk.
#' @param file_prefix A character string to use as the prefix for the filenames.
#'
#' @return No return value, called for side effects.
#' @examples
#' mysim <- list(simlrX = list(v = list(
#'   data1 = data.frame(matrix(rnorm(147 * 171), nrow = 147, ncol = 171)),
#'   data2 = data.frame(matrix(rnorm(156 * 171), nrow = 156, ncol = 171))
#' )))
#' write_simlr_data_frames(mysim$simlrX$v, "output")
#' @export
write_simlr_data_frames <- function(data_list, file_prefix) {
  for (i in seq_along(data_list)) {
    # Generate a filename using the index
    file_name <- paste0(file_prefix, "_", names(data_list)[i], "_simlr.csv")
    
    # Write the data frame to disk
    write.csv(data_list[[i]], file_name, row.names = TRUE)
  }
}

#' Read a list of data frames from disk with specific naming convention
#'
#' This function reads a list of data frames from disk into a list,
#' assuming the files are named with a common prefix and the names of the data frames.
#' It converts the column named `X` to the row names of the read data frame.
#'
#' @param file_prefix A character string used as the prefix for the filenames.
#' @param data_names A character vector of names for the data frames.
#'
#' @return A list of data frames read from disk with the column named `X` set as row names.
#' @examples
#' # data_names <- c("data1", "data2")
#' # data_list <- read_simlr_data_frames(file_prefix = "output", data_names = data_names)
#' # dim(data_list[[1]])
#' # dim(data_list[[2]])
#' @export
read_simlr_data_frames <- function(file_prefix, data_names) {
  data_list <- list()
  
  for (name in data_names) {
    # Generate the filename using the prefix and data names
    file_name <- paste0(file_prefix, "_", name, "_simlr.csv")
    
    # Read the data frame from disk
    if ( file.exists( file_name ) ) {
      df <- read.csv(file_name, row.names = 1)
      
      # Convert the column named `X` to row names, if it exists
      if ("X" %in% colnames(df)) {
        rownames(df) <- df$X
        df <- df[ , !colnames(df) %in% "X"]
      }
      
      # Store the data frame in the list
      data_list[[name]] <- df
    }
  }
  
  return(data_list)
}



#' Apply simlr matrices to an existing data frame and combine the results
#'
#' This function takes a list of matrices, applies each matrix via matrix multiplication
#' to an existing data frame, and combines the resulting projections with the original data frame.
#'
#' @param existing_df An existing data frame to which the matrices will be applied.
#' @param matrices_list A list of matrices read from CSV files.
#' @param n_limit NULL or integer that can limit the number of projections
#' @param robust boolean
#' @param verbose boolean
#'
#' @return A list including (entry one) data frame with the original data frame combined with the projections (entry two) the new column names
#' @export
#' @examples
#' matrices_list <- list(
#'   matrix1 = matrix(rnorm(147 * 171), nrow = 147, ncol = 171),
#'   matrix2 = matrix(rnorm(147 * 156), nrow = 147, ncol = 156)
#' )
#' existing_df <- data.frame(matrix(rnorm(147 * 5), nrow = 147, ncol = 5))
#' # combined_df <- apply_simlr_matrices(existing_df, matrices_list)
apply_simlr_matrices <- function(existing_df, matrices_list, n_limit=NULL, robust=FALSE, verbose=FALSE ) {
  newnames=c()
  for (name in names(matrices_list)) {
    if ( verbose ) print(name)
    # Ensure the matrix multiplication is valid
    locnames = rownames( matrices_list[[name]] )
    edfnames = colnames(existing_df) 
    inames = intersect( locnames, edfnames )
    if ( length(inames) > 0 ) {
      # Perform matrix multiplication
      imat = data.matrix(existing_df[,inames])
      if ( robust ) imat = robustMatrixTransform( imat )
      projection <- as.data.frame(
        imat %*% data.matrix(matrices_list[[name]][inames,]))
      ##################################################
      # Update column names to reflect the matrix name #
      colnames(projection) = paste0( name, colnames( matrices_list[[name]] ) )
      # Combine the projections with the existing data frame
      if ( !is.null(n_limit )  ) {
        projection=projection[,1:n_limit]
      }
      newnames=c(newnames,colnames(projection))
      existing_df <- cbind(existing_df, projection)
      if ( verbose ) {
        print( inames )
        print( colnames(projection) )
        print(tail(colnames(existing_df)))
      }
    } else {
      warning(paste("Number of columns in existing data frame does not match number of rows in matrix", name))
    }
  }
  
  return( list(extendeddf=existing_df, newcolnames=newnames))
}



#' Select Important Variables Using Partial Correlation Matrix
#'
#' This function computes the partial correlation matrix for the specified columns in a dataset and selects the most important variables based on the sum of absolute values of partial correlations.
#'
#' @param data A data frame containing the dataset.
#' @param cols A vector of column names or indices to be included in the analysis.
#' @param threshold A numeric value between 0 and 1 indicating the proportion of variables to retain. Defaults to 0.2 (20%).
#'
#' @return A vector of column names corresponding to the most important variables.
#' @export
#'
#' @examples
#' set.seed(123)
#' qqq <- data.frame(matrix(rnorm(100 * 10), ncol = 10))
#' colnames(qqq) <- paste0("Var", 1:10)
#' tempcols <- colnames(qqq)
#' important_vars <- select_important_variables(qqq, tempcols, threshold = 0.3)
#' print(important_vars)
select_important_variables <- function(data, cols, threshold = 0.2) {
  # Compute the correlation matrix
  mycor <- cor(na.omit(data[, cols]))
  
  # Compute the inverse of the correlation matrix (precision matrix)
  precision_matrix <- solve(mycor)
  
  # Compute partial correlations from the precision matrix
  partial_cor_matrix <- -cov2cor(precision_matrix)
  diag(partial_cor_matrix) <- 1  # Set the diagonal to 1 for partial correlations
  
  # Sum the absolute values of partial correlations for each variable
  variable_importance <- apply(abs(partial_cor_matrix), 1, sum)
  
  # Determine the threshold for selection
  num_vars_to_select <- round(length(variable_importance) * threshold)
  
  # Select the most important variables
  important_vars <- names(sort(variable_importance, decreasing = TRUE)[1:num_vars_to_select])
  
  return(important_vars)
}