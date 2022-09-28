####################################################
####   R code for GENERAL ISA FUNCTIONS    #########
####################################################
####   Andrea, April 2017                  #########
####################################################

# options(java.parameters = "-Xmx128000m")

# library(pheatmap)
# library(RColorBrewer)
# library(gplots)
# library(xlsx)
# library(ggplot2)
# library(scales) 
# library(splines)
# library(MASS)
# library(gridExtra)
# 
# library(vegan)
# library(mclust)
# library(psych)
# library(GPArotation)
# library(nFactors)
# library(FactoMineR)
# 
# library(XLConnect)
# library(ISAnalytics)
# library(VennDiagram)
# library(reshape2)
# library(tagcloud)

###############################################################
## functions
###############################################################

###############################################################
#' @title Rename columns of a matrix (df) by metadata
#' 
#' @author Andrea Calabria
#' @details version 1.0 
#'
#' @rdname renameDfColumnsByMetadataField
#' @docType methods
#' @aliases renameDfColumnsByMetadataField
#'
#' @param df an input dataframe of IS matrix.
#' @param metadata_df a dataframe of metadata, row names must be an overset of the columns of the input df.
#' @param key_field a string/character for the field name (a colname of the metadata_df)
#' @param starting_data_col_index an integer value >0 with the index of the data columns (excluding annotations, etc.)
#'
#' @return a dataframe with the same size of the input df but with different colnames, taken from metadata_df
#' @usage todo
#'
###############################################################
###############################################################
renameDfColumnsByMetadataField <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # goal: given the input df and the metadata df, for each colname, look for that value in the metadata index and return the key field to aggregate. then apply the aggregation function of plyr
  require(plyr)
  # colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:length(colnames(df))]) # get column names of df into a new df
  # rownames(colID_df) <- colnames(df)[starting_data_col_index:length(colnames(df))]
  # colmetadata_df <- merge(colID_df, label_metadata, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # rownames(colmetadata_df) <- colmetadata_df$Row.names
  # # get names from the key field
  # colmetadata_df_outnames <- colmetadata_df[colnames(df)[starting_data_col_index:length(colnames(df))],c(key_field)]
  # # create the new out df
  # tmp_df <- df
  # colnames(tmp_df) <- c(colnames(df)[1:starting_data_col_index-1], as.character(colmetadata_df_outnames))
  # # return
  # return (tmp_df)
  
  # check that all cols of df are in metadata
  cols_to_rename <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  if (length(intersect(rownames(metadata_df), cols_to_rename)) == length(cols_to_rename)) {  # perfect case, go ahead
    names(df) <- c(colnames(df[1:starting_data_col_index-1]), 
                   as.character(metadata_df[cols_to_rename, key_field])
                   )
    return (df)
  } else {
    message(paste("[AP]\tERROR: not all cols of the file are present in metadata."))
  }
}

aggregateDfColumnsByName <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and aggregate by name
  # so far,  this is a by-hand procedure... waiting for better ideas
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  if (length(grep(key_field, colnames(colmetadata_df)))>1) {
    message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  }
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

###############################################################
#' @title Get a new df of IS where columnd are collapsed by metadata, collapsing by using a custom
#' 
#' @author Andrea Calabria
#' @details version 1.0 
#'
#' @rdname aggregateDfColumnsByName_customFun
#' @docType methods
#' @aliases aggregateDfColumnsByName_customFun
#'
#' @param df an input dataframe patient_iss.
#' @param metadata_df a df of metadata such as the association file.
#' @param key_field a string for the column name (aka field) to search for in the metadata df.
#' @param starting_data_col_index an integer value, the index of the first column to analyze (1..n)
#' @param number_of_last_cols_to_remove an integer value (0..n) of number of last columns to remove from the analysis
#' @param myfunction is the function name for the aggregation (default: "sum")
#' 
#' @return df
#' @usage TODO
#' @description TODO
#'
###############################################################
aggregateDfColumnsByName_customFun <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0, myfunction = sum, ...) {
  # given the input df, look up the column  names and aggregate by name
  # so far,  this is a by-hand procedure... waiting for better ideas
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  if (length(grep(key_field, colnames(colmetadata_df)))>1) {
    message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  }
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        # out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {myfunction(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        # this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {myfunction(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

###############################################################
#' @title Get a new df of IS where columnd are collapsed by metadata, collapsing by using a custom (mean)
#' 
#' @author Andrea Calabria
#' @details version 1.0 
#'
#' @rdname aggregateDfColumnsByName_mean
#' @docType methods
#' @aliases aggregateDfColumnsByName_customFun
#'
#' @param df an input dataframe patient_iss.
#' @param metadata_df a df of metadata such as the association file.
#' @param key_field a string for the column name (aka field) to search for in the metadata df.
#' @param starting_data_col_index an integer value, the index of the first column to analyze (1..n)
#' @param number_of_last_cols_to_remove an integer value (0..n) of number of last columns to remove from the analysis
#' @param myfunction is the function name for the aggregation (default: "sum")
#' 
#' @return df
#' @usage TODO
#' @description TODO
#'
###############################################################
aggregateDfColumnsByName_mean <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and aggregate by name
  # so far,  this is a by-hand procedure... waiting for better ideas
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  if (length(grep(key_field, colnames(colmetadata_df)))>1) {
    message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  }
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        # out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {mean(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        # this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {mean(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

computeAbundancePercentage_fromMetadata <- function(df = df, metadata_df, key_field = "QuantificationSum", starting_data_col_index = 1, last_data_col_index = NA){
  # given a matrix of data only, compute abundance
  message(paste("[AP]\tCalculating abundance from input matrix"))
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
  }
  # first check if you have any emopty column (no rows with values >0)
  overall_sum <- apply(df[starting_data_col_index:last_data_col_index], 2, function(x) {sum(x, na.rm = TRUE)})
  if (min(overall_sum) == 0) {
    message(paste("[AP]\t\tWARNING: You have empty cols (cols with no values)"))
  }
  # do abundance
  abundance_df <- data.frame(sapply(colnames(df[starting_data_col_index:last_data_col_index]), function(x) {
      # standard relative %:
      .quantification_sum <- metadata_df[x, key_field]
      if (.quantification_sum > 0) {
        ( df[,x]*100 / (metadata_df[x, key_field]) )
      } else {
        ( df[,x]*0 )
      }
    } )
  ) # do the abundance
  rownames(abundance_df) <- rownames(df) # change names
  # do warnings
  if (min(abundance_df)<0) {
    message(paste("[AP]\t\tThe abundance is producing values <0 (!!!)"))
  }
  # return out df
  return (abundance_df)
}

computeAbundancePercentage <- function(df = df, metadata_df = NULL, key_field = "QuantificationSum", starting_data_col_index = 1, last_data_col_index = NA){
  # given a matrix of data only, compute abundance
  message(paste("[AP]\tCalculating abundance from input matrix"))
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
  }
  source_df <- df[starting_data_col_index:last_data_col_index]
  # compute internally the sum df, and call the column of sum as key_field
  # fill NA with 0
  source_df[is.na(source_df)] <- 0
  metadata_df <- data.frame( 
    "QuantificationSum" = apply(source_df, 2, function(x) {sum(x, na.rm=T)} ),
    "NumIS" = apply(source_df, 2, function(x) {length(x[x>0])} )
  )
  names(metadata_df) <- c(key_field, "NumIS")
  # first check if you have any emopty column (no rows with values >0)
  overall_sum <- apply(source_df, 2, function(x) {sum(x, na.rm = TRUE)})
  if (min(overall_sum) == 0) {
    message(paste("[AP]\t\tWARNING: You have empty cols (cols with no values)"))
  }
  # do abundance
  abundance_df <- data.frame(sapply(colnames(source_df), function(x) {
    # standard relative %:
    .quantification_sum <- metadata_df[x, key_field]
    if (.quantification_sum > 0) {
      ( source_df[,x]*100 / (metadata_df[x, key_field]) )
    } else {
      ( source_df[,x]*0 )
    }
  } )
  ) # do the abundance
  rownames(abundance_df) <- rownames(df) # change names
  colnames(abundance_df) <- colnames(source_df)
  # do warnings
  if (min(abundance_df)<0) {
    message(paste("[AP]\t\tThe abundance is producing values <0 (!!!)"))
  }
  # return out df
  return (abundance_df)
}

sortRowsByMaxValue <- function(df = df, starting_data_col_index = 1, last_data_col_index = NA, annotation_columns = c()){
  # TODO: fix a problem with single colum data and single row data
  # given a df, compute the max value by row and return the df sorted by max value
  message(paste("[AP]\tCalculating max from input matrix"))
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
  }
  out_df <- df[starting_data_col_index:last_data_col_index]
  # compute internally the sum df, and call the column of sum as key_field
  # fill NA with 0
  out_df[is.na(out_df)] <- 0
  # single column df?
  if (ncol(out_df) == 1) {
    message(paste("[AP]\t-- Single data col input matrix"))
    single_col_name <- colnames(df[starting_data_col_index:last_data_col_index])
    # add a foo col to manage df
    out_df <- cbind(out_df, data.frame("removeThisCol" = rep(999, nrow(out_df))))
    # out_df <- data.frame(order(-out_df[single_col_name]), )
    # single row df?
    if (nrow(out_df) > 1) {
      out_df <- out_df[order(-out_df[single_col_name]), ]
      names(out_df) <- single_col_name
    } else {
      message(paste("[AP]\t-- Single data col input matrix and single row (!!!)"))
    }
    
  } else {
    out_df <- cbind(out_df, data.frame("df_row_max_tmp" = apply(out_df, 1, max)) )
    out_df <- out_df[order(-(out_df["df_row_max_tmp"])), setdiff(colnames(out_df), c("df_row_max_tmp")) ]
  } # if (ncol(out_df) == 1)
  # add annotation
  rows_to_keep <- rownames(out_df)
  if (length(annotation_columns) > 0) {
    out_df <- cbind(df[rows_to_keep, c(annotation_columns)], out_df)
  }
  # return df
  return (out_df[setdiff(colnames(out_df), c("removeThisCol"))])
}

# sortRowsByColumnValues <- function(df, value_orientation = "descending", compact = TRUE, scan_columns_from_left = TRUE, replace_0_with_NA = FALSE) {
#   # given the input df, sort row df by column, orientation user defined. Column order will remain the same
#   # message("[AP]\tConverting input df by adding 0 to NA")
#   # df[is.na(df)] <- 0 # avoid here inside NA
#   # if (compact) {
#   #   df <- compactDfByColumns(df)
#   # }
#   list_container <- list() # init the list of resulting object
#   message(paste("[AP]\tNow looping over columns to comput sharing results"))
#   c_index <- 1
#   for (c in colnames(df)) {
#     message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
#     slice_df <- df[which(df[c]>0),] # slice df
#     list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])})) ) # the relativ eprecentage of contaminations
#     c_index <- c_index + 1
#   }
#   names(list_container) <- colnames(df) # rename list objects
#   r <- as.data.frame(list_container)
#   rownames(r) <- colnames(slice_df)
#   ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
#   if (left_to_rigth_reading_output) {
#     message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
#     # if need to prune_selected_rows
#     if (prune_selected_rows & length(pruning_rows_labels)>0) {
#       message(paste("[AP]\tRemoving selected rows"))
#       r <- r[,!colnames(r) %in% pruning_rows_labels]
#       return (t(r))
#     } else {
#       return (t(r))
#     }
#   } 
#   else {
#     message(paste("[AP]\tOutput matrix orientation: Top-Down"))
#     if (prune_selected_rows & length(pruning_rows_labels)>0) {
#       message(paste("[AP]\tRemoving selected rows"))
#       r <- r[!rownames(r) %in% pruning_rows_labels,]
#       return (r)
#     } else {
#       return (r)
#     }
#   }
# }

getTopIS <- function(df, starting_data_col_index = 1, last_data_col_index = NA, output_IS_annotation_col = c("GeneName"), output_IS_annotation_colname = c("Gene Name"), output_IS_value_colname = c("Quantification"), number_of_output_elements = 10) {
  ### given a source df (from any quantification or abundance), this function returns the top N IS as dataframe with the associated gene name and its value
  # do prelimiar operations and checks
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
    message(paste("[AP]\tConsidering df columns from", starting_data_col_index, colnames(df[starting_data_col_index]), "to", last_data_col_index, colnames(df[last_data_col_index])))
  }
  if (output_IS_annotation_col %in% colnames(df)) {
    message(paste("[AP]\tOutput column for annotation, found", output_IS_annotation_col))
  } else {
    message(paste("[AP]\tERROR: Output column for annotation NOT found", output_IS_annotation_col))
  }
  
  # in case of NA elements, fill with 0
  df[is.na(df)] <- 0 # avoid here inside NA
  
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  
  # get the top N genes (by value)
  top_genes_byabundance_withvalue <- NULL 
  for (c in colnames(df[starting_data_col_index:last_data_col_index])) {
    message(paste("[AP]\t-> processing the colum", c, "(", as.character(c_index), "of", as.character(last_data_col_index-starting_data_col_index), as.character(round(c_index*100/(last_data_col_index-starting_data_col_index),2)), "%]"))
    slice_df <- df[which(df[c]>0), c(output_IS_annotation_col, c)] # slice df
    slice_df_sorted <- slice_df[order(-slice_df[c]), ]
    slice_df_topgenes <- slice_df_sorted[1:number_of_output_elements, ]
    
    # get output with value
    if (length(top_genes_byabundance_withvalue) == 0) {
      # init out df
      top_genes_byabundance_withvalue <- slice_df_topgenes
      names(top_genes_byabundance_withvalue) <- c(paste(c, output_IS_annotation_colname), paste(c, output_IS_value_colname))
    }
    else {
      actual_colnames <- colnames(top_genes_byabundance_withvalue)
      top_genes_byabundance_withvalue <- cbind(top_genes_byabundance_withvalue, slice_df_topgenes )
      names(top_genes_byabundance_withvalue) <- c(actual_colnames, paste(c, output_IS_annotation_colname), paste(c, output_IS_value_colname))
    }
    
    c_index <- c_index + 1
  }
  
  # out
  return (top_genes_byabundance_withvalue)
  
  # # get the top 50 genes (by abundance)
  # top_genes_byabundance <- NULL 
  # top_genes_byabundance_withperc <- NULL 
  # first_n_elements_to_return <- 30000
  # for (k in levels(factor(abundance_df_melt_dataonly$UniqueSample)) ) {
  #   message(paste("[AP]\tProcessing sample", k))
  #   # slice df by value
  #   slice <- abundance_df_melt_dataonly[which(abundance_df_melt_dataonly$UniqueSample == k),]
  #   # keep the genes
  #   sorted_slice <- slice[order(-slice$Abundance), ]
  #   sorted_slice_topgenes <- sorted_slice[1:first_n_elements_to_return, c("closest_gene")]
  #   sorted_slice_topgenes_withperc <- sorted_slice[1:first_n_elements_to_return, c("closest_gene", "Abundance")]
  #   # get output
  #   if (length(top_genes_byabundance) == 0) {
  #     # init out df
  #     top_genes_byabundance <- as.data.frame(sorted_slice_topgenes)
  #     names(top_genes_byabundance) <- k
  #   }
  #   else {
  #     actual_colnames <- colnames(top_genes_byabundance)
  #     top_genes_byabundance <- cbind(top_genes_byabundance, as.data.frame(sorted_slice_topgenes) )
  #     names(top_genes_byabundance) <- c(actual_colnames, k)
  #   }
  #   # get output with percentage
  #   if (length(top_genes_byabundance_withperc) == 0) {
  #     # init out df
  #     top_genes_byabundance_withperc <- as.data.frame(sorted_slice_topgenes_withperc)
  #     names(top_genes_byabundance_withperc) <- c(paste(k, "Gene symbol"), paste(k, "Abundance"))
  #   }
  #   else {
  #     actual_colnames <- colnames(top_genes_byabundance_withperc)
  #     top_genes_byabundance_withperc <- cbind(top_genes_byabundance_withperc, as.data.frame(sorted_slice_topgenes_withperc) )
  #     names(top_genes_byabundance_withperc) <- c(actual_colnames, paste(k, "Gene symbol"), paste(k, "Abundance"))
  #   }
  # }
  # write.xlsx2(top_genes_byabundance, file = paste("analyses/01.abundance/MLD.preclinical.aggregatedTriplicates_source_postAllColllisions.top", first_n_elements_to_return, "IS_onSamplePerc.xlsx", sep = ""), sheetName = paste("Top", first_n_elements_to_return, "IS abundance"), append = TRUE, row.names = FALSE)
  # write.xlsx2(top_genes_byabundance_withperc, file = paste("analyses/01.abundance/MLD.preclinical.aggregatedTriplicates_source_postAllColllisions.top", first_n_elements_to_return, "IS_onSamplePerc.xlsx", sep = ""), sheetName = paste("Top", first_n_elements_to_return, "IS abundance with val"), append = TRUE, row.names = FALSE)
  # 
}


getTopIS_stdAnnotation_returnID <- function(df, starting_data_col_index = 1, 
                                            last_data_col_index = NA, 
                                            annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                                            threshold_of_first_top = 0, 
                                            number_of_output_elements = 10) {
  ### given a source df (from any quantification or abundance), this function returns the top N IS as dataframe with the associated gene name and its value
  # do prelimiar operations and checks
  # evaluate column length
  if (is.na(last_data_col_index)) {
    last_data_col_index <- length(colnames(df))
    message(paste("[AP]\tConsidering df columns from", starting_data_col_index, colnames(df[starting_data_col_index]), "to", last_data_col_index, colnames(df[last_data_col_index])))
  }
  if (length(annotation_cols) == 0) {
    message(paste("[AP]\tERROR: Annotation NOT found, paramenter 'annotation_cols'"))
    stopifnot(length(annotation_cols) > 0)
  }
  
  # in case of NA elements, fill with 0
  df[is.na(df)] <- 0 # avoid here inside NA
  
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  
  # get the top N genes (by value)
  top_genes_byabundance_withvalue <- NULL 
  for (c in colnames(df[starting_data_col_index:last_data_col_index])) {
    message(paste("[AP]\t-> processing the colum", c, "(", as.character(c_index), "of", as.character(last_data_col_index-starting_data_col_index), as.character(round(c_index*100/(last_data_col_index-starting_data_col_index),2)), "%]"))
    slice_df <- df[which(df[c]>0), c(annotation_cols, c)] # slice df
    slice_df_sorted <- slice_df[order(-slice_df[c]), ]
    slice_df_topgenes <- slice_df_sorted[1:number_of_output_elements, ]
    max_abundance_value <- max(slice_df_topgenes[c], na.rm = T)
    #annotation_cols_tomelt <- c(annotation_cols, c("OncoGene", "TumorSuppressor", "Onco1_TS2", "ClinicalRelevance", "DOIReference", "KnownClonalExpension"))
    slice_df_topgenes_melt <- melt(data = slice_df_topgenes, id.vars = annotation_cols, variable.name = "Sample", na.rm = T, value.name = "Abundance")
    
    c_index <- c_index + 1
    
    # get output with value
    if (length(top_genes_byabundance_withvalue) == 0) {
      # init out df
      if (max_abundance_value >= threshold_of_first_top) {
        top_genes_byabundance_withvalue <- slice_df_topgenes_melt
      } # if (max_abundance_value >= threshold_of_first_top) 
    }
    else {
      if (max_abundance_value >= threshold_of_first_top) {
        top_genes_byabundance_withvalue <- rbind(top_genes_byabundance_withvalue, slice_df_topgenes_melt )
      } # if (max_abundance_value >= threshold_of_first_top) 
    } # if (length(top_genes_byabundance_withvalue) == 0)
    
    
  } # for (c in colnames(df[starting_data_col_index:last_data_col_index]))
  
  # out
  return (top_genes_byabundance_withvalue)
}


reduceDfbyRowID <- function(df, id_df) {
  # given the input dfs, merge by rows
  message(paste("[AP]\tReducing rows of the source df (df) by the guide df (id_df)"))
  out_df <- merge(df, id_df, by=0, all.y = TRUE) # get the output merge df of annotations with all x metadata
  rownames(out_df) <- c(out_df$Row.names)
  # return df
  return (out_df)
}

getTopHitGenes <- function(df, data_column_to_use, group_by_colums, maxGenesToReturn = 10, column_suffix = c("Gene Symbol", "Annotated IS")) {
  # input: df, data columns to extract results
  # logics: get the top HIT genes (not by abundance)
  # todo: control on max number of returning elements:: if dat aelements are < than max -> check it
  # output: df with values
  result_df <- NULL 
  for (k in data_column_to_use) {
    message(paste("[AP]\tPocessing column(s)", k))
    aggregated_slice <- aggregate(df[which(df[k] > 0), c(k)], by=list(df[which(df[k] > 0), group_by_colums]), FUN=function(x) {length(x[x>0])})
    aggregated_slice_sorted <- aggregated_slice[order(-c(aggregated_slice$x)),]
    # get output with values
    if (length(result_df) == 0) {
      # init out df
      result_df <- as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ])
      names(result_df) <- c(paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
    else {
      actual_colnames <- colnames(result_df)
      result_df <- cbind(result_df, as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ]) )
      names(result_df) <- c(actual_colnames, paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
  }
  return (result_df)
}

getTopCIS_withIScount <- function(df, data_column_to_use, group_by_colum, maxGenesToReturn = 10, column_suffix = c("Gene Symbol", "Annotated IS")) {
  # input: df of melted data (!!!!), data values to extract results
  # logics: get first the top CIS and their MAX score, then combine their number of annotated landed IS (IS count)
  # todo: control on max number of returning elements:: if dat aelements are < than max -> check it
  # output: df with values
  result_df <- NULL 
  for (k in data_column_to_use) {
    message(paste("[AP]\tPocessing column(s)", k))
    aggregated_slice <- aggregate(df[which(df[k] > 0), c(k)], by=list(df[which(df[k] > 0), group_by_colums]), FUN=function(x) {length(x[x>0])})
    aggregated_slice_sorted <- aggregated_slice[order(-c(aggregated_slice$x)),]
    # get output with values
    if (length(result_df) == 0) {
      # init out df
      result_df <- as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ])
      names(result_df) <- c(paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
    else {
      actual_colnames <- colnames(result_df)
      result_df <- cbind(result_df, as.data.frame(aggregated_slice_sorted[1:maxGenesToReturn, ]) )
      names(result_df) <- c(actual_colnames, paste(k, column_suffix[1]), paste(k, column_suffix[2]))
    }
  }
  return (result_df)
}


computeBackgroundRow_byThreshold <- function( source_df = df, threshold, background_row_label = "Background", no_negative_values = TRUE ) {
  # from a df, given the input threshold, 
  message(paste("[AP]\tCompute background row given a threshold", threshold))
  out_df <- NULL # init our df
  # do checks and perform selection
  if (ncol(source_df) == 1) {
    #abovet <- data.frame('AboveThreshold' = source_df[source_df >= threshold])
    single_data_col <- colnames(source_df) # the data col
    source_df <- cbind(source_df, data.frame("AboveThreshold" = rep(1, nrow(source_df))) )
    abovet <- source_df[which(source_df[single_data_col] >= threshold), ]
    # keep only IS/region/ID to track
    df_totrack <- abovet[single_data_col]
    # compute background row
    bg_row <- as.data.frame( round(100 - sum(df_totrack), digits = 10) )
    names(bg_row) <- colnames(df_totrack)
    rownames(bg_row) <- background_row_label
    # add background row to the sliced df (transposing the bg df)
    out_df <- rbind(df_totrack, bg_row)
    
  } else {
    # compute what is above threshold
    abovet <- data.frame('AboveThreshold' = apply(source_df, 1, function(x) {if (!is.na(x)) { ifelse(max(x, na.rm = TRUE) >= threshold, 1, 0)}}) )
    # bind dataframes
    df_comp <- cbind(source_df, abovet)
    # keep only IS/region/ID to track
    df_totrack <- df_comp[which(df_comp$AboveThreshold == 1), colnames(source_df)]
    # compute background row
    bg_row <- as.data.frame( round(100 - apply(df_totrack, 2, sum), digits = 10) )
    names(bg_row) <- background_row_label
    # add background row to the sliced df (transposing the bg df)
    out_df <- rbind(df_totrack, t(bg_row)) 
  }
  # return df if it does not contain errors
  if (min(out_df) < 0 & no_negative_values) {
    message(paste("[AP]\t[ERROR]\tThe matrix contains values <0 (!!!)"))
    return (NA)
  } else {
    return (out_df)
  }
}

addSharingColumn <- function ( source_df = df, column_label = "NumberSharedIS"  ) {
  # given the input df, ad dth ecolumn of shared elements with a custom label
  # from a df, given the input threshold, 
  message(paste("[AP]\tCompute sharing from the input df"))
  # fill 0s 
  source_df[is.na(source_df)] <- 0
  # do sharing
  out_df <- as.data.frame(c(apply(source_df, 1, function(x) { ifelse( length(x[x>0])==0, NA, length(x[x>0]) ) } ) ))
  # rename the column name
  names(out_df) <- column_label
  # replace NAs
  out_df[out_df==0] <- NA
  # get output
  return (out_df)
}

computeBackgroundRow_bySharedIS <- function( source_df = df, threshold = 1, background_row_label = "Background", shared_column_label = "NumberSharedIS", remove_shared_column = TRUE, no_negative_values = TRUE) {
  # sharing must be already available
  if (shared_column_label %in% colnames(source_df)) {
    # ok proceed
    message(paste("[AP]\tCompute background row given the number of shared IS above the min threshold (>)", threshold))
    
    # keep only IS/region/ID to track
    df_totrack <- source_df[which(source_df[shared_column_label] > threshold), ]
    # compute background row
    bg_row <- as.data.frame( round(100 - apply(df_totrack, 2, sum), digits = 10) )
    names(bg_row) <- background_row_label
    
    # add background row to the sliced df (transposing the bg df)
    out_df <- rbind(df_totrack, t(bg_row)) 
    
    # return df if it does not contain errors
    if (min(out_df) < 0 & no_negative_values) {
      message(paste("[AP]\t[ERROR]\tThe matrix contains values <0 (!!!)"))
      return (NA)
    } else {
      
      # remove shared IS column?
      if (remove_shared_column) {
        return (out_df[, setdiff(colnames(out_df), c(shared_column_label))])
      } else {
        return (out_df)
      }
      
    }
  } else {
    message(paste("[AP]\tERROR: the required column of the shared elements is missing. Check the input label name (here provided", shared_column_label, ") or just run the appropriate function to retrieve the shared column."))
  } # if else
}

computeBackgroundRow_byThresholdORSharedIS <- function( source_df = df, sharing_threshold = 1, percentage_threshold = 0.5, background_row_label = "Background", shared_column_label = "NumberSharedIS", remove_additional_column = TRUE, no_negative_values = TRUE) {
  # sharing must be already available
  if (shared_column_label %in% colnames(source_df)) {
    # ok proceed
    message(paste("[AP]\tCompute background row given the number of shared IS above the min threshold (>)", sharing_threshold, " OR percentage threshold", percentage_threshold))
    
    # compute what is above threshold
    abovet <- data.frame('AboveThreshold' = apply(source_df[, setdiff(colnames(source_df), c(shared_column_label))], 1, function(x) {if (!is.na(x)) { ifelse(max(x, na.rm = TRUE) >= percentage_threshold, 1, 0)}}) )
    # bind dataframes
    df_comp <- cbind(source_df, abovet)
    
    # keep only IS/region/ID to track
    df_totrack <- df_comp[which(df_comp[shared_column_label] > sharing_threshold | df_comp$AboveThreshold == 1), ]
    # compute background row
    bg_row <- as.data.frame( round(100 - apply(df_totrack, 2, sum), digits = 10) )
    names(bg_row) <- background_row_label
    
    # add background row to the sliced df (transposing the bg df)
    out_df <- rbind(df_totrack, t(bg_row)) 
    
    # remove shared IS column?
    if (remove_additional_column) {
      # return df if it does not contain errors
      out_df <- out_df[, setdiff(colnames(out_df), c(shared_column_label, "AboveThreshold"))]
      if (min(out_df) < 0 & no_negative_values) {
        message(paste("[AP]\t[ERROR]\tThe output matrix contains values <0 (!!!)"))
        return (out_df)
      } else {  
        return (out_df[, setdiff(colnames(out_df), c(shared_column_label, "AboveThreshold"))])
      } 
    } else {
      return (out_df)
    }
    
    # return df if it does not contain errors
    if (min(out_df) < 0 & no_negative_values) {
      message(paste("[AP]\t[ERROR]\tThe output matrix contains values <0 (!!!)"))
      return (out_df)
    } else {
      
    }
  } else {
    message(paste("[AP]\tERROR: the required column of the shared elements is missing. Check the input label name (here provided", shared_column_label, ") or just run the appropriate function to retrieve the shared column."))
  } # if else
}

# addBackgroundRow_byNumberOfIS <- function( source_df = df, is_to_visualize = 10 ) {
#   # from a df, given the input number of IS to visualize, 
#   message(paste("[AP]\tCompute background row given the numnber of IS to visualize independently from each sample", is_to_visualize))
#   
#   return (out_df)
# }

getSharedIS <- function(df, compact = TRUE, left_to_rigth_reading_output = TRUE, prune_selected_rows = FALSE, pruning_rows_labels = c()) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  for (c in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[c]>0),] # slice df
    list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])})) ) # the relativ eprecentage of contaminations
    c_index <- c_index + 1
  }
  names(list_container) <- colnames(df) # rename list objects
  r <- as.data.frame(list_container)
  rownames(r) <- colnames(slice_df)
  ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
  if (left_to_rigth_reading_output) {
    message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
    # if need to prune_selected_rows
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[,!colnames(r) %in% pruning_rows_labels]
      return (t(r))
    } else {
      return (t(r))
    }
  } 
  else {
    message(paste("[AP]\tOutput matrix orientation: Top-Down"))
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[!rownames(r) %in% pruning_rows_labels,]
      return (r)
    } else {
      return (r)
    }
  }
}

fullDistinctColTrackingIS_amongColumnsByName_0mat <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and track IS by col name
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  # Same of trackIS_amongColumnsByName BUT:
  # - the INPUT matrix must have 0s -> no biniary version of the matrix and no SUM but LENGTH (count)
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])) # get column names of df into a new df
  rownames(colID_df) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[out_df>=1] <- 1 # binary data
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {length(x[x>0])}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[slice_cols]
        this_slice[this_slice>=1] <- 1
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {length(x[x>0])})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and track IS by col name
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  
  # binarize df, do not remove here the initial cols, this is not required
  message(paste("[AP]\tBinarizing df"))
  df_binarycount <- data.frame(t(apply(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)], 1, function(x) {ifelse(x>0, 1, 0)})))
  names(df_binarycount) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df_binarycount)) # get column names of df into a new df
  rownames(colID_df) <- colnames(df_binarycount)
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df_binarycount[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df_binarycount[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName_0mat <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and track IS by col name
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  # Same of trackIS_amongColumnsByName BUT:
  # - the INPUT matrix must have 0s -> no biniary version of the matrix and no SUM but LENGTH (count)
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])) # get column names of df into a new df
  rownames(colID_df) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[out_df>=1] <- 1 # binary data
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {length(x[x>0])}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[slice_cols]
        this_slice[this_slice>=1] <- 1
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {length(x[x>0])})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName_getRatio <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0, replace_NA_with_0 = TRUE) {
  # given the input df, look up the column  names and track IS by col name -> return the ratio of tracked IS among available samples
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  
  # binarize df, do not remove here the initial cols, this is not required
  if (replace_NA_with_0) {
    message(paste("[AP]\tConvert df NA as 0"))
    df[is.na(df)] <- 0
  } 
  df_binarycount <- df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)] # the name binary is not meaningful here
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df_binarycount)) # get column names of df into a new df
  rownames(colID_df) <- colnames(df_binarycount)
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df_binarycount[slice_cols]
        out_df[is.na(out_df)] <- 0
        out_df[out_df>0] <- 1
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {length(x[x>0])/length(slice_cols)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df_binarycount[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        this_slice[this_slice>0] <- 1
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {length(x[x>0])/length(slice_cols)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName_getRatio_doBinMat <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column  names and track IS by col name -> return the ratio of tracked IS among available samples
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  
  # binarize df, do not remove here the initial cols, this is not required
  message(paste("[AP]\tBinarizing df"))
  df_binarycount <- data.frame(t(apply(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)], 1, function(x) {ifelse(x>0, 1, 0)})))
  names(df_binarycount) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df_binarycount)) # get column names of df into a new df
  rownames(colID_df) <- colnames(df_binarycount)
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df_binarycount[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)/length(slice_cols)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df_binarycount[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df_binarycount[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)/length(slice_cols)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName_getMax <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column names and track IS by col name -> return the MAX value
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])) # get column names of df into a new df
  rownames(colID_df) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {max(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {max(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

trackIS_amongColumnsByName_getFoldMaxonSecond <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0) {
  # given the input df, look up the column names and track IS by col name -> return the MAX value
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])) # get column names of df into a new df
  rownames(colID_df) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        out_df[out_df>0] <- 1 # all values >1: place 1 that is the fold
        names(out_df) <- k
      } 
      else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {
          .max = max(x, na.rm = TRUE)
          tmp <- sort(x, decreasing = T)
          ifelse(length(x[x>0])>1, .max/tmp[2], 1)
        }))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } # if (length(out_df) == 0) { # if this is the first loop
    else {
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        this_slice[this_slice>0] <- 1
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } 
      else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {
          .max = max(x, na.rm = TRUE)
          tmp <- sort(x, decreasing = T)
          ifelse(length(x[x>0])>1, .max/tmp[2], 1)
        })) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

disaggregateDfColumnsByName <- function(df, metadata_df, key_field, groupedby_value, row_id_list, annotation_columns_to_include) {
  # given the input df, look up the column  names and trace back the columns to disaggregate by name (1 to N)
  metadata_output_cols <- rownames(metadata_df[which(metadata_df[key_field] == groupedby_value), ])
  output_cols <- c( intersect(annotation_columns_to_include, colnames(df)) , intersect(metadata_output_cols, colnames(df)) )
  message(paste("[AP]\tThese will be the output columns:"))
  message(paste("\t\n", output_cols))
  # return dataframe
  if (length(output_cols) > 1) {
    out_df <- df[row_id_list, output_cols]
  } else {
    out_df <- data.frame(df[row_id_list, output_cols])
    rownames(out_df) <- row_id_list
    names(out_df) <- output_cols
  }
  return (out_df)  
}

compactDfByColumns <- function(df) {
  # compct df by pruning columns without any value (all rows of these columns are NA or 0)
  # 1. identify null columns
  # 2. return the same df witout these columns
  colnames(df) 
  notnull_cols <- apply(df, 2, function(x) {nrow(df) - length(x[ x == 0 | is.na(x) ])})
  null_col_names <- as.vector(names(notnull_cols[notnull_cols==0]))
  null_col_index <- colnames(df) %in% null_col_names
  return (df[!null_col_index])
}

compactDfByRows <- function(df, data_columns, annotation_columns = NULL) {
  # compact df by pruning rows without any value (meaning that no value for any columns, all NA or 0)
  # 1. identify null rows
  # 2. return the same df witout these rows
  # Example usage: t <- compactDfByRows(gdf, data_columns = grep("BThal", colnames(gdf)))
  message("[AP]\tCompacting by rows (pruning rows without values)")
  if (length(data_columns) == 0) {
    stop("[AP]\tERROR: Arguments data_columns must be specified. Try with grep, example: compactDfByRows(gdf, data_columns = grep('BThal', colnames(gdf)))")
  } 
  col_sum <- data.frame("SUM" = apply(df[data_columns], 1, function(x) {abs(sum(x, na.rm = TRUE))}), 
                        "ID" = rownames(df[data_columns])
                        ) # use absolute values so that only positive values are the rows to keep
  rows_to_keep <- rownames(col_sum[which(col_sum$SUM > 0),])
  if (length(annotation_columns) > 0) {
    return (df[rows_to_keep, c(annotation_columns, data_columns)])
  } else {
    return (df[rows_to_keep, data_columns]) 
  }
}


# valore di riferimento: numero di IS condivise tra campioni -> fissato un paziente, per ogni tempo|campione fai rapporto di quante sono le IS condivise prma e dopo (diventa un rapporto N_condivise/N_osservate)
getSharedISratio <- function(df, compact = TRUE, left_to_rigth_reading_output = TRUE, prune_selected_rows = FALSE, pruning_rows_labels = c()) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  for (c in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[c]>0),] # slice df
    list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])/nrow(slice_df)})) ) # the relativ eprecentage of contaminations
    c_index <- c_index + 1
  }
  names(list_container) <- colnames(df) # rename list objects
  r <- as.data.frame(list_container)
  rownames(r) <- colnames(slice_df)
  ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
  if (left_to_rigth_reading_output) {
    message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
    # if need to prune_selected_rows
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[,!colnames(r) %in% pruning_rows_labels]
      return (t(r))
    } else {
      return (t(r))
    }
  } 
  else {
    message(paste("[AP]\tOutput matrix orientation: Top-Down"))
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[!rownames(r) %in% pruning_rows_labels,]
      return (r)
    } else {
      return (r)
    }
  }
}

# come la gemella getSharedISratio, cambia il tipo di calcolo
getSharedISnumber <- function(df, compact = TRUE, left_to_rigth_reading_output = TRUE, prune_selected_rows = FALSE, pruning_rows_labels = c()) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  list_container <- list() # init the list of resulting object
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  for (c in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", c, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[c]>0),] # slice df
    list_container <- c( list_container, data.frame(c = apply(slice_df, 2, function(x) {length(x[x>0])})) ) # the relativ eprecentage of contaminations
    c_index <- c_index + 1
  }
  names(list_container) <- colnames(df) # rename list objects
  r <- as.data.frame(list_container)
  rownames(r) <- colnames(slice_df)
  ### NB: in questo momento la matrice (df r) ha lettura alto basso (!!!!), e non sinistra destra. in base all'opzione invertila nel return
  if (left_to_rigth_reading_output) {
    message(paste("[AP]\tTranspose output matrix: Left->Right orientation"))
    # if need to prune_selected_rows
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[,!colnames(r) %in% pruning_rows_labels]
      return (t(r))
    } else {
      return (t(r))
    }
  } 
  else {
    message(paste("[AP]\tOutput matrix orientation: Top-Down"))
    if (prune_selected_rows & length(pruning_rows_labels)>0) {
      message(paste("[AP]\tRemoving selected rows"))
      r <- r[!rownames(r) %in% pruning_rows_labels,]
      return (r)
    } else {
      return (r)
    }
  }
}

getQuantileElements_byProbs <- function(df, columnNames_to_exclude, annotation_count_column = "Onco1_TS2", gene_column = "closest_gene", probs = c(0, 0.25, 0.5, 0.75, 1) ) {
  # given an input df (only data columns), compute the quartiles and extract elements (row IDs) by classes of quartiles
  out_df <- NULL
  onco_byqclass <- NULL
  df_dataonly <- df[setdiff(colnames(df), columnNames_to_exclude)]
  # 1. compute quantiles
  message(paste("[AP]\tCompute quantile by data column(s):", colnames(columnNames_to_exclude)))
  quantile_df <- apply(df_dataonly, 2, function(x) {quantile(x[x>0], probs = probs)})
  # 2. get sorted
  for (k in colnames(df_dataonly)) {
    message(paste("[AP]\tPocessing column(s)", k))
    slice_df <- df[which(df[k]>0), c(columnNames_to_exclude, k)]
    # do quantiles 
    for (index_qclass in seq(1, length(as.numeric(gsub("%", "", rownames(quantile_df))))-1) ) { # if you need to convert in number: as.numeric(gsub("%", "", rownames(quantile_df))) 
      qclass <- as.numeric(gsub("%", "", rownames(quantile_df)[index_qclass]))/100
      next_qclass <- as.numeric(gsub("%", "", rownames(quantile_df)[index_qclass + 1]))/100
      slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile(slice_df[,k], probs = probs)[index_qclass] & slice_df[k] < quantile(slice_df[,k], probs = probs)[index_qclass + 1])
      if (index_qclass == length(as.numeric(gsub("%", "", rownames(quantile_df))))-1 ) {
        slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile(slice_df[,k], probs = probs)[index_qclass])
      }
      # check this slice first
      if( nrow(slice_df_byqclass) == 0 ) {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = qclass,
                                       "ThresholdClassValue_lower" = quantile(slice_df[,k], probs = probs)[index_qclass],
                                       "ThresholdClassValue_upper" = quantile(slice_df[,k], probs = probs)[index_qclass + 1],
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = NA,
                                       "Ratio_OncoTSOnElements" = NA
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = qclass,
                                                            "ThresholdClassValue_lower" = quantile(slice_df[,k], probs = probs)[index_qclass],
                                                            "ThresholdClassValue_upper" = quantile(slice_df[,k], probs = probs)[index_qclass + 1],
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = NA,
                                                            "Ratio_OncoTSOnElements" = NA
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      else {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = qclass,
                                       "ThresholdClassValue_lower" = quantile(slice_df[,k], probs = probs)[index_qclass],
                                       "ThresholdClassValue_upper" = quantile(slice_df[,k], probs = probs)[index_qclass + 1],
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                       "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]) / nrow(slice_df_byqclass)
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = qclass,
                                                            "ThresholdClassValue_lower" = quantile(slice_df[,k], probs = probs)[index_qclass],
                                                            "ThresholdClassValue_upper" = quantile(slice_df[,k], probs = probs)[index_qclass + 1],
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                                            "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]) / nrow(slice_df_byqclass)
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      
    } # for (index_qclass in seq(1, length(as.numeric(gsub("%", "", rownames(quantile_df))))-1) )
  } # for (k in colnames(df_dataonly))
  # get results
  return (onco_byqclass)
}

getQuantileElements <- function(df, columnNames_to_exclude, annotation_count_column = "Onco1_TS2", gene_column = "closest_gene") {
  # given an input df (only data columns), compute the quartiles and extract elements (row IDs) by classes of quartiles
  out_df <- NULL
  onco_byqclass <- NULL
  df_dataonly <- df[setdiff(colnames(df), columnNames_to_exclude)]
  # 1. compute quantiles
  message(paste("[AP]\tCompute quantile by data column(s):", colnames(columnNames_to_exclude)))
  quantile_df <- apply(df_dataonly, 2, function(x) {quantile(x[x>0])})
  # 2. get sorted
  for (k in colnames(df_dataonly)) {
    message(paste("[AP]\tPocessing column(s)", k))
    slice_df <- df[which(df[k]>0), c(columnNames_to_exclude, k)]
    # do quantiles 
    for (index_qclass in seq(1, length(as.numeric(gsub("%", "", rownames(quantile_df))))-1) ) { # if you need to convert in number: as.numeric(gsub("%", "", rownames(quantile_df))) 
      qclass <- as.numeric(gsub("%", "", rownames(quantile_df)[index_qclass]))/100
      next_qclass <- as.numeric(gsub("%", "", rownames(quantile_df)[index_qclass + 1]))/100
      slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile(slice_df[,k], qclass) & slice_df[k] < quantile(slice_df[,k], next_qclass))
      if (index_qclass == length(as.numeric(gsub("%", "", rownames(quantile_df))))-1 ) {
        slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile(slice_df[,k], qclass))
      }
      # check this slice first
      if( nrow(slice_df_byqclass) == 0 ) {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = qclass,
                                       "ThresholdClassValue_lower" = quantile(slice_df[,k], qclass),
                                       "ThresholdClassValue_upper" = quantile(slice_df[,k], next_qclass),
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = NA,
                                       "Ratio_OncoTSOnElements" = NA
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = qclass,
                                                            "ThresholdClassValue_lower" = quantile(slice_df[,k], qclass),
                                                            "ThresholdClassValue_upper" = quantile(slice_df[,k], next_qclass),
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = NA,
                                                            "Ratio_OncoTSOnElements" = NA
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      else {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = qclass,
                                       "ThresholdClassValue_lower" = quantile(slice_df[,k], qclass),
                                       "ThresholdClassValue_upper" = quantile(slice_df[,k], next_qclass),
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                       "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]) / nrow(slice_df_byqclass)
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = qclass,
                                                            "ThresholdClassValue_lower" = quantile(slice_df[,k], qclass),
                                                            "ThresholdClassValue_upper" = quantile(slice_df[,k], next_qclass),
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                                            "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]) / nrow(slice_df_byqclass)
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      
    } # for (index_qclass in seq(1, length(as.numeric(gsub("%", "", rownames(quantile_df))))-1) )
  } # for (k in colnames(df_dataonly))
  # get results
  return (onco_byqclass)
}

getQuantileElements_byInterval <- function(df, columnNames_to_exclude, annotation_count_column = "Onco1_TS2", gene_column = "closest_gene", number_of_intervals = 4 ) {
  # given an input df (only data columns), compute the quartiles and extract elements (row IDs) by classes of quartiles. in this case, find the MIN and MAX that set up the interval, then divide the interval in segments and compute the annotations
  out_df <- NULL
  onco_byqclass <- NULL
  df_dataonly <- df[setdiff(colnames(df), columnNames_to_exclude)]
  # 1. compute quantiles
  message(paste("[AP]\tCompute quantile by data column(s):", colnames(columnNames_to_exclude)))
  # get the intervals for each data column. This matrix will have rownames with the integer interval and the column with the sample, cells contains the separation value between two adjacent intervals
  quantile_df <- apply(df_dataonly, 2, function(x) {
    step <- (max(x[x>0]) - min(x[x>0])) / number_of_intervals;
    seq(min(x[x>0]), max(x[x>0]), step);
  })
  # 2. get sorted
  for (k in colnames(df_dataonly)) {
    message(paste("[AP]\tPocessing column(s)", k))
    slice_df <- df[which(df[k]>0), c(columnNames_to_exclude, k)]
    # do quantiles 
    for (index_qclass in seq(1, nrow(quantile_df)-1) ) { # if you need to convert in number: as.numeric(gsub("%", "", rownames(quantile_df))) 
      qclass <- quantile_df[index_qclass, k]
      next_qclass <- quantile_df[index_qclass + 1, k]
      # do slice, but note that if this is the last loop cycle, include also the last element!
      slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile_df[index_qclass, k] & slice_df[k] < quantile_df[index_qclass + 1, k])
      if (index_qclass == nrow(quantile_df)-1 ) {
        slice_df_byqclass <- subset(slice_df, slice_df[k] >= quantile_df[index_qclass, k])
      }
      # check this slice first
      if( nrow(slice_df_byqclass) == 0 ) {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = index_qclass,
                                       "ThresholdClassValue_lower" = quantile_df[index_qclass, k],
                                       "ThresholdClassValue_upper" = quantile_df[index_qclass + 1, k],
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = NA,
                                       "Ratio_OncoTSOnElements" = NA
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = index_qclass,
                                                            "ThresholdClassValue_lower" = quantile_df[index_qclass, k],
                                                            "ThresholdClassValue_upper" = quantile_df[index_qclass + 1, k],
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = NA,
                                                            "Ratio_OncoTSOnElements" = NA
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      else {
        # add results in the dataframe
        if (length(onco_byqclass) == 0 ) {
          onco_byqclass <- data.frame( "Sample" = k,
                                       "QuantileClass" = index_qclass,
                                       "ThresholdClassValue_lower" = quantile_df[index_qclass, k],
                                       "ThresholdClassValue_upper" = quantile_df[index_qclass + 1, k],
                                       "TotalElements" = nrow(slice_df_byqclass),
                                       "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                       "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ])/nrow(slice_df_byqclass)
                                       # "ElementList" = as.vector( slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), gene_column] )
          )
        }
        else { # if (length(onco_byqclass) == 0 )
          actual_colnames <- rownames(onco_byqclass)
          onco_byqclass <- rbind(onco_byqclass, data.frame( "Sample" = k,
                                                            "QuantileClass" = index_qclass,
                                                            "ThresholdClassValue_lower" = quantile_df[index_qclass, k],
                                                            "ThresholdClassValue_upper" = quantile_df[index_qclass + 1, k],
                                                            "TotalElements" = nrow(slice_df_byqclass),
                                                            "AnnotatedElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ]),
                                                            "Ratio_OncoTSOnElements" = nrow(slice_df_byqclass[which(slice_df_byqclass[annotation_count_column]>0), ])/nrow(slice_df_byqclass)
          )
          )
        } # if (length(onco_byqclass) == 0 )   
      } # if( nrow(slice_df_byqclass) == 0 )
      
    } # for (index_qclass in seq(1, length(as.numeric(gsub("%", "", rownames(quantile_df))))-1) )
  } # for (k in colnames(df_dataonly))
  # get results
  return (onco_byqclass)
}

getMinDate_amongColumnsInDateByName_binmat <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0, col_prefix = "TrackingMinDate_") {
  # given the input df, look up the column  names and track IS by col name
  # just binarize the matrix and then leave the function operating in SUM as aggregateDfColumnsByName
  # for each col name, slice df and apply function
  # Same of trackIS_amongColumnsByName BUT:
  # - the INPUT matrix must have 0s -> no biniary version of the matrix and no SUM but LENGTH (count)
  
  # acquire metadata
  colID_df <- data.frame("colID" = colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])) # get column names of df into a new df
  rownames(colID_df) <- colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])
  # merge with metadata
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_searchkey <- colmetadata_df[,c(key_field)]
  
  # now get the keys to seach
  search_keys <- as.character(levels(factor(colmetadata_df_searchkey)))
  
  #message(paste("[AP]\tThese will be the output/aggregated columns:", search_keys))
  out_df <- NULL
  for (k in search_keys) {
    # get col IDs from metadata
    #label_metadata[which(label_metadata[grep(key_field, colnames(label_metadata))] == k), c("colID")]
    message(paste("[AP]\tPocessing search key", k))
    # search the columns
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[key_field] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      # if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
      #     out_df <- df[slice_cols]
      #     names(out_df) <- paste(col_prefix, k, sep = '')
      # }
      # else {
      #     out_df <- data.frame(
      #     "slicemindate" = c(apply(df[slice_cols], 1, function(x) { ifelse(min(as.Date(x), na.rm = T)==Inf, NA, min(as.Date(x), na.rm = T)) } ))
      #     )
      #     names(out_df) <- paste(col_prefix, k, sep = '')
      # } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
      out_df <- data.frame(
        "slicemindate" = c(apply(df[slice_cols], 1, function(x) { ifelse(min(as.Date(x), na.rm = T)==Inf, NA, min(as.Date(x), na.rm = T)) } ))
      )
      names(out_df) <- paste(col_prefix, k, sep = '')
      
    } else {
      # if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
      #     actual_colnames <- colnames(out_df)
      #     this_slice <- df[slice_cols]
      #     names(this_slice) <- paste(col_prefix, k, sep = '')
      #     out_df <- cbind(out_df, this_slice)
      # }
      # else {
      #     actual_colnames <- colnames(out_df)
      #     this_slice <- data.frame(
      #     "slicemindate" = c(apply(df[, slice_cols], 1, function(x) { ifelse(min(as.Date(x), na.rm = T)==Inf, NA, min(as.Date(x), na.rm = T)) } ))
      #     )
      #     names(this_slice) <- paste(col_prefix, k, sep = '')
      #     out_df <- cbind(out_df, this_slice )
      # } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
      actual_colnames <- colnames(out_df)
      this_slice <- data.frame(
        "slicemindate" = c(apply(df[slice_cols], 1, function(x) { ifelse(min(as.Date(x), na.rm = T)==Inf, NA, min(as.Date(x), na.rm = T)) } ))
      )
      names(this_slice) <- paste(col_prefix, k, sep = '')
      out_df <- cbind(out_df, this_slice )
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}

rewriteMatrix <- function(full_df, project_id, patient_file_xlsx, id_cols, extra_cols = c(), return_df = FALSE ) {
  message(paste("[AP]\tCreating Matrix structure and file for Project: ", project_id))
  search_path <- paste("^", project_id, "_", sep = "") # decide the search pattern
  sourcecol <- names(full_df) %in% c( id_cols, grep(search_path, names(full_df), value=TRUE), grep("^Collision against", names(full_df), value=TRUE), extra_cols )
  #sourcecol <- names(full_df) %in% c( id_cols, grep(search_path, names(full_df), value=TRUE), grep("^Collision against", names(full_df), value=TRUE) )
  patient_column <- paste("SequenceCountSum_", project_id, sep = "") # create the variable with the name of the predefined column
  col_number <- grep(patient_column, colnames(full_df)) # find corresponding column numbers within the data structure
  patient_formatrix <- full_df[ which(full_df[,col_number]>0), sourcecol ] # slice source data frame (full df) both by rows (only this patient IS) and columns
  names(patient_formatrix)[names(patient_formatrix) %in% c( grep(paste("^Collision against ", project_id, sep = ""), names(patient_formatrix), value=TRUE) )] <- "all" # now rename the patient specific column with the name "all"
  patient_formatrix[patient_formatrix==0] <- NA
  write.xlsx2(patient_formatrix, file=patient_file_xlsx, row.names=FALSE, sheetName = "ISmatrix_gauss_&collisions")
  
  # once finished to export dataset you may export/return the df
  if (return_df) {
    return (patient_formatrix)
  }
}

rewriteMatrixForISA <- function(full_df, project_id, patient_file_xlsx, id_cols ) {
  message(paste("[AP]\tCreating Matrix structure and file for Project: ", project_id))
  search_path <- paste("^", project_id, "_", sep = "") # decide the search pattern
  sourcecol <- names(full_df) %in% c( id_cols, 'GeneName', 'GeneStrand', grep(search_path, names(full_df), value=TRUE), grep("^Collision against", names(full_df), value=TRUE) )
  #sourcecol <- names(full_df) %in% c( id_cols, grep(search_path, names(full_df), value=TRUE), grep("^Collision against", names(full_df), value=TRUE) )
  patient_column <- paste("SequenceCountSum_", project_id, sep = "") # create the variable with the name of the predefined column
  col_number <- grep(patient_column, colnames(full_df)) # find corresponding column numbers within the data structure
  patient_formatrix <- full_df[ which(full_df[,col_number]>0), sourcecol ] # slice source data frame (full df) both by rows (only this patient IS) and columns
  names(patient_formatrix)[names(patient_formatrix) %in% c( grep(paste("^Collision against ", project_id, sep = ""), names(patient_formatrix), value=TRUE) )] <- "all" # now rename the patient specific column with the name "all"
  patient_formatrix[patient_formatrix==0] <- NA
  write.xlsx2(patient_formatrix, file=patient_file_xlsx, row.names=FALSE, sheetName = "ISmatrix_gauss_&collisions")
}

# given an input vector of labels (matrix column labels), get all metadata info in the same input order
queryMetadata <- function(labels, label_metadata_df, metadata_to_export = c("patient_id", "tissue", "sample", "timepoint", "cell_type", "relative_blood_perc", "separation_purity", "VCN","cell_marker", "pool_id", "n_reactions", "dna_ng_used", "dna_extraction_date", "linear_pcr_date", "firstpcr_pcr_date", "fusionprimer_pcr_date", "fusionprimer_pcr_date", "pool_date", "ltrvalid_reads_noplasmid", "ISS_IDENTIFIED", "Patient_SelectedDate"), transpose_output = FALSE) {
  return_df <- label_metadata_df[labels, metadata_to_export]
  if (transpose_output) {
    return (t(return_df))
  } else {
    return (return_df)
  }
  
}

# slice df to extract data for DB import
getTSVForDB <- function(full_df, project_id, assignment_label = "AssignedByDateTo_", selected_fields = c("Chr", "Pos", "Strand", "X_repeat", "X_cons_seq"), add_projec_id = FALSE ) {
  message(paste("[AP]\tCreating sequence dataset for IS only derived from patient", project_id))
  assignment_col <- grep(paste(assignment_label, project_id, sep = ""), colnames(full_df))
  df_patient_fordb <- full_df[ which(full_df[assignment_col] == 1 ), selected_fields ] # slice source data frame (full df) both by rows (only this patient IS) and columns
  df_repeats_groups <- data.frame(do.call('rbind', strsplit(as.character(df_patient_fordb[, "X_repeat"]), ";") )) # Alu_SINE:Alu
  df_repeats_class <- data.frame(do.call('rbind', strsplit(as.character(df_repeats_groups[,1]), ":") )) # Alu_SINE       Alu
  df_repeats_family <- data.frame(do.call('rbind', strsplit(as.character(df_repeats_class[,1]), "_") )) # Alu
  # return df with project_id or without
  if (add_projec_id) {
    df_project_id <- data.frame("Project_ID" = rep(project_id, nrow(df_patient_fordb)))
    df_result <- cbind(df_patient_fordb, df_repeats_family[,1], df_repeats_class[,1], df_repeats_groups[,1], df_project_id)
    names(df_result) <- c("Chr", "Pos", "Strand", "Repeat_Annotation", "Consensus_Sequence", "Repeat_Family", "Repeat_Classes", "Repeat_Group", "Project_ID")
    # v_repeat_annotation <- gsub(";", "", df_patient_fordb[, "X_repeat"])
    return (df_result)
  }
  else {
    df_result <- cbind(df_patient_fordb, df_repeats_family[,1], df_repeats_class[,1], df_repeats_groups[,1])
    names(df_result) <- c("Chr", "Pos", "Strand", "Repeat_Annotation", "Consensus_Sequence", "Repeat_Family", "Repeat_Classes", "Repeat_Group")
    # v_repeat_annotation <- gsub(";", "", df_patient_fordb[, "X_repeat"])
    return (df_result)
  }
}


### do intersection 
intersectMarkers_absoluteCounts <- function(df, with_prefix = "CD34", against_prefix = "CD15_PB", empty_element_zero = TRUE) {
  if (empty_element_zero) { # if the matrix contains 0 instead of null values (for unobseverd values)
    df_slice <- df[c(grep(with_prefix, colnames(df)), grep(against_prefix, colnames(df)))] # slice source df
    nIS_byCol <- data.frame("nIS" = apply(df_slice, 2, function(x) {length(x[x>0])})) # get a df of the number of IS per column
    for (i in seq(1:length(df_slice))) { # loop over sliced columns
      df_index_slice <- df_slice[which(df_slice[i]>0),] # select IS derived from this column only
      sharing_index_slice <- data.frame(apply(df_index_slice, 2, function(x) {length(x[x>0])})) # result of sharing counts among all the other columns
      names(sharing_index_slice) <- colnames(df_slice[i]) # rename cols
      # build the output df of counts
      if (i == 1) {
        return_df <- sharing_index_slice
      } else {
        return_df <- cbind(return_df, sharing_index_slice)
      }
    }
  } else { # else is implied to use NA
    df_slice <- df[c(grep(with_prefix, colnames(df)), grep(against_prefix, colnames(df)))] # slice source df
    nIS_byCol <- data.frame("nIS" = apply(df_slice, 2, function(x) {length(x[!is.na(x)])})) # get a df of the number of IS per column
    for (i in seq(1:length(df_slice))) { # loop over sliced columns
      df_index_slice <- df_slice[which(df_slice[i]>0),] # select IS derived from this column only
      sharing_index_slice <- data.frame(apply(df_index_slice, 2, function(x) {length(x[!is.na(x)])})) # result of sharing counts among all the other columns
      names(sharing_index_slice) <- colnames(df_slice[i]) # rename cols
      # build the output df of counts
      if (i == 1) {
        return_df <- sharing_index_slice
      } else {
        return_df <- cbind(return_df, sharing_index_slice)
      }
    } 
  }
  return (list("sharing_counts" = return_df, "nIS" = nIS_byCol))
}

### do intersection 
intersectByName_absoluteCounts <- function(df, col_names, empty_element_zero = TRUE) {
  if (empty_element_zero) { # if the matrix contains 0 instead of null values (unobseverd)
    df_slice <- df[col_names %in% colnames(df)]
    nIS_byCol <- data.frame("nIS" = apply(df_slice, 2, function(x) {length(x[x>0])})) # get a df of the number of IS per column
    for (i in seq(1:length(df_slice))) { # loop over sliced columns
      df_index_slice <- df_slice[which(df_slice[i]>0),] # select IS derived from this column only
      sharing_index_slice <- data.frame(apply(df_index_slice, 2, function(x) {length(x[x>0])})) # result of sharing counts among all the other columns
      names(sharing_index_slice) <- colnames(df_slice[i])
      # build the output df of counts
      if (i == 1) {
        return_df <- sharing_index_slice
      } else {
        return_df <- cbind(return_df, sharing_index_slice)
      }
    }
  } else { # else is implied to use NA
    df_slice <- df[col_names %in% colnames(df)]
    nIS_byCol <- data.frame("nIS" = apply(df_slice, 2, function(x) {length(x[!is.na(x)])})) # get a df of the number of IS per column
    for (i in seq(1:length(df_slice))) { # loop over sliced columns
      df_index_slice <- df_slice[which(df_slice[i]>0),] # select IS derived from this column only
      sharing_index_slice <- data.frame(apply(df_index_slice, 2, function(x) {length(x[!is.na(x)])})) # result of sharing counts among all the other columns
      names(sharing_index_slice) <- colnames(df_slice[i])
      # build the output df of counts
      if (i == 1) {
        return_df <- sharing_index_slice
      } else {
        return_df <- cbind(return_df, sharing_index_slice)
      }
    } 
  }
  return (list("sharing_counts" = return_df, "nIS" = nIS_byCol))
}

# do ratio of the intersection on the union using the intersection results
# prerequisite: results of intersectMarkers_absoluteCounts or intersectByName_absoluteCounts
computeRatioOnUnion <- function(sharing_object) {
  for (col in colnames(sharing_object$sharing_counts)) {
    rows <- c()
    for (row in rownames(sharing_object$sharing_counts)) {
      rows <- c(rows, sharing_object$sharing_counts[row, col]*100/(sharing_object$nIS[col,] + sharing_object$nIS[row,] - sharing_object$sharing_counts[row, col]))
    }
    if (col == colnames(sharing_object$sharing_counts)[1]) {
      col_ratio <- as.data.frame(rows)
      rownames(col_ratio) <- rownames(sharing_object$sharing_counts)
    } else {
      col_ratio <- cbind(col_ratio, as.data.frame(rows)) 
    }
  }
  names(col_ratio) <- colnames(sharing_object$sharing_counts)
  return (col_ratio)
}

# do ratio of the intersection on one of the sets, the columns
# prerequisite: results of intersectMarkers_absoluteCounts or intersectByName_absoluteCounts
# logics:: given set X and set Y, the intersection produces the sets A, B, C as A = X&!Y, B = X&Y, C = !X&Y. here the denominator is Y (=> B+C)
computeRatioOnColumn <- function(sharing_object) {
  col_ratio <- data.frame(sapply(colnames(sharing_object$sharing_counts), function(x) { (sharing_object$sharing_counts[,x]/sharing_object$nIS[x,])*100 } ))
  rownames(col_ratio) <- colnames(sharing_object$sharing_counts)
  return (col_ratio)
}

# do cumulative counts (given a direction that depends on the input given order of the columns)
getBinaryCumulative_ProgressiveUnion <- function(df) {
  # to overcome the problem of missing values (non 0 values) for a column (resulting in apply errors), I coded a dummy column... it works but it is very rude
  df_fake <- cbind(df, data.frame("thisTMPsum" = apply(df, 1, function(x) {length(x)})))
  cumulated_progressive <- as.data.frame(t(apply(df_fake, 1, function(x) {.first_obs <- min(which(x>0)); c(rep(0, .first_obs-1), rep(1, length(df_fake)-.first_obs+1)) })))
  cumulated_progressive <- cumulated_progressive[-length(cumulated_progressive)]
  names(cumulated_progressive) <- colnames(df)
  return (cumulated_progressive)
}

getBinaryCumulative_ProgressiveUnion_byColumnsName <- function(df, metadata_df, key_field, starting_data_col_index = 2, number_of_last_cols_to_remove = 0, decreasing_sort = F) {
  # given the input df, look up the column names and extract cols of the same group, then sort them and for each group do union cumulative count
  # so far,  this is a by-hand procedure... waiting for better ideas
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  # get names from the key field
  colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  
  # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  if (length(grep(key_field, colnames(colmetadata_df)))>1) {
    message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  }
  
  output_cols <- as.character(levels(factor(colmetadata_df_outnames)))
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        out_df[out_df > 0] <- 1
        # names(out_df) <- k
      } else {
        # first sort columns by name
        slice_cols <- sort(x = slice_cols, decreasing = decreasing_sort)
        # out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)}))
        out_df <- getBinaryCumulative_ProgressiveUnion(df = df[slice_cols])
        # names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        # this_slice <- df[, slice_cols] # corrected on Nov 3
        this_slice <- df[slice_cols]
        this_slice[is.na(this_slice)] <- 0
        this_slice[this_slice > 0] <- 1
        cumulative_df <- getBinaryCumulative_ProgressiveUnion(df = this_slice)
        out_df <- cbind(out_df, cumulative_df )
        # names(out_df) <- c(actual_colnames, k)
      } else {
        # first sort columns by name
        slice_cols <- sort(x = slice_cols, decreasing = decreasing_sort)
        # actual_colnames <- colnames(out_df)
        cumulative_df <- getBinaryCumulative_ProgressiveUnion(df = df[, slice_cols])
        out_df <- cbind(out_df, cumulative_df )
        # names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
}


# do cumulative counts with ONLY new elements (given a direction that depends on the input given order of the columns)
getBinaryCumulative_ProgressiveNew <- function(df) {
  # to overcome the problem of missing values (non 0 values) for a column (resulting in apply errors), I coded a dummy column... it works but it is very rude
  df_fake <- cbind(df, data.frame("thisTMPsum" = apply(df, 1, function(x) {length(x)})))
  cumulated_newonly <- as.data.frame(t(apply(df_fake, 1, function(x) {.first_obs <- min(which(x>0)); c(rep(0, .first_obs-1), x[.first_obs], rep(0, length(df_fake)-.first_obs)) })))
  cumulated_newonly <- cumulated_newonly[-length(cumulated_newonly)]
  names(cumulated_newonly) <- colnames(df)
  return (cumulated_newonly)
}

# get union for columns by grep/text search
getColumnUnion <- function(df, column_search_key) {
  df_slice <- df[grep(column_search_key, colnames(df))]
  union_df <- data.frame(apply(df_slice, 1, function(x) {rsum <- sum(x, na.rm = TRUE); ifelse(rsum > 0, 1, rsum)}))
  names(union_df) <- paste(column_search_key, "_union", sep = "")
  return (union_df)
}
# get sum for columns by grep/text search
getColumnSum <- function(df, column_search_key) {
  df_slice <- df[grep(column_search_key, colnames(df))]
  sum_df <- data.frame(apply(df_slice, 1, function(x) {sum(x, na.rm = TRUE)}))
  names(sum_df) <- paste(column_search_key, "_sum", sep = "")
  return (sum_df)
}

# fill missing columns (as time point) with NA or 0 values for the actual rows
# output_suffix_tohave:: the last fidl of the label names (that is the time point string as text)
completeDFwithMissingColumns <- function(df, output_suffix_tohave, fill_value = 0) {
  present_col_suffix <- levels(factor(t(as.data.frame(strsplit(colnames(df), "_")))[,3])) # colnames(df)
  missing_col_suffix_to_include <- setdiff(output_suffix_tohave, present_col_suffix)
  # message(paste(missing_col_suffix_to_include, collapse = " "))
  # for each suffix create a dedicated column of 0 in the output as col_group * missing cols
  if (length(missing_col_suffix_to_include) == 0) {
    # if no missing values, return input
    return (df)
  } else {
    # message("I have to add columns...")
    # else add cols
    col_groups <- levels(factor(apply(t(as.data.frame(strsplit(colnames(df), "_"))), 1, function(x) {paste(x[1], "_", x[2], sep = "")})))
    # message(paste(col_groups, collapse = " "))
    list_container <- list()
    list_container_names <- c()
    for (col in col_groups) {
      for (missing_timepoint in missing_col_suffix_to_include) {
        new_col <- rep(fill_value, nrow(df))
        list_container <- c(list_container, data.frame(new_col))
        list_container_names <- c(list_container_names, paste(col, missing_timepoint, sep = "_"))
      }
    }
    added_cols_to_df <- as.data.frame(list_container)
    names(added_cols_to_df) <- list_container_names
    # now bind df and sort it
    complete_df <- cbind(df, added_cols_to_df)
    sorted_complete_df <- complete_df[,order(names(complete_df))]
    return (sorted_complete_df)
  }
}

# group columns by summming values, using search key
groupColumns_byPrefix <- function(df, search_key, new_label = "NEW_LABEL") {
  # use a prefix to search for df columns
}
# group columns by summming values, listing all input columns
groupColumns_byName <- function(df, col_names, new_label = "NEW_LABEL") {
  # use a prefix to search for df columns
}


getElementCategories <- function(df, starting_data_col_index = 2, number_of_last_cols_to_remove = 0, prune_known_categories = c(0)) {
  # questa funzione si potrebbe rimpiazzare alla grande con la semplice riga "k <- as.data.frame(apply(df, 2, function(x) {table(x)} ))" se l'input avesse tutte le categorie per ogni colonna
  out_df <- NULL
  # from df, find categories, that are the unique elements
  unique_elements <- unique(unlist(df))
  message(paste(c("[AP]\tFound these categories:", unique_elements), collapse=" ") )
  message(paste(c("[AP]\tThe following categories are user-defined un-wanted:", prune_known_categories), collapse=" ") )
  # if required by user, remove un-wanted categories
  unique_elements <- setdiff(unique_elements, prune_known_categories)
  message(paste(c("[AP]\tThe following categories will be analyzed (these will be the output columns):", unique_elements), collapse=" ") )
  for (step in unique_elements) {
    message(paste("[AP]\t\tRetrieving elements in category", step) )
    this_step_df <- as.data.frame(apply(df, 2, function(x) {length(x[x==step])} ))
    names(this_step_df) <- as.character(step)
    if (length(out_df) == 0) { # if this is the first loop
      out_df <- this_step_df
    } else {
      out_df <- cbind(out_df, this_step_df)
    } # if else if (length(out_df) == 0)
    this_step_df <- NULL
  }
  return (out_df)
}

###############################################################
#' @title Recalibrate IS using a BED map
#' 
#' @author Andrea Calabria
#' @details version 1.1 (used in clonal tracking), 2018-04-19
#'
#' @rdname recalibrateIS
#' @docType methods
#' @aliases recalibrateIS
#'
#' @param df_is an input dataframe even with annotation
#' @param recalibration_map a df importing the BED file annotation from reference, thus having the following structure of columns: "from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance"
#' @param recalibrate_IS_distance integer of the distance in bp to distinguish tow ajacent IS. Default = 7, that is <=7.
#' @param id_cols a vector of column names for annotation. Default: c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
#'
#' @return a df with the same number of columns of the input df (df_is) but at most the same number of rows.
#' @usage TODO
#' @note : 
#' 
#' The following code is for testing the function, assuming that in the folder "data" you have the appropriate files.
#' t.expected_recalibratedIS <- read.csv("source/test/test.matrix.foo.expected_after_recalibration.csv", check.names = F, header = T, sep = "\t", row.names = 1) # load the recalibrated expected matrix
#' id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
#' t.in <- read.csv("source/test/test.matrix.foo.new_format.csv", check.names = F, header = T, sep = "\t")
#' rownames(t.in) <- apply(t.in[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
#' t.rmap <- read.csv("source/test/test.matrix.foo.new_format.withRef_WindowMatrix.bed", check.names = F, header = F, sep = "\t")
#' rownames(t.rmap) <- t.rmap$V4
#' names(t.rmap) <- c("from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance")
#' t.recalibratedIS <- recalibrateIS(df_is = t.in, recalibration_map = t.rmap)
#' 
#'
###############################################################
recalibrateIS <- function(df_is, recalibration_map, recalibrate_IS_distance = 7, id_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand") ) {
  library(gplots)
  ##### -------- RECALIBRATION ------- #####
  message(paste("[AP]\tRecalibration"))
  # strategy:
  # - slice df with only ID to remap / recalibrate (the "from" field)
  # - correct IS (apply the ID of the "to" field)
  # - extract the "to" fils IS ID
  # - fuse IS with (potential) existing IDs by summing elements (cells)
  # - remove IS ID from-to from dfallp
  # - add the new aggregated/fused matrix to dfallp
  
  # load recalibration matrix
  message(paste("[AP]\t\tLoading BED file and recalibration data, from file"))
  # define a col order
  col_order <- colnames(df_is)
  # recalibration_map <- read.csv(file_recalibration_map_bed, header=F, fill=T, sep='\t', check.names = FALSE)
  # rownames(recalibration_map) <- recalibration_map$V4
  names(recalibration_map) <- c("from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance")
  to_recalibrate_df <- recalibration_map[which(recalibration_map$from_start != recalibration_map$to_start & abs(recalibration_map$distance) <= recalibrate_IS_distance ),]
  # do the annotation from ID
  to_recalibrate_df_toID_annotation <- as.data.frame(t(as.data.frame(lapply(strsplit(as.character(to_recalibrate_df$to_id), "_", fixed = TRUE), function(x) {c(x)}))))
  names(to_recalibrate_df_toID_annotation) <- id_cols
  to_recalibrate_df <- cbind(to_recalibrate_df, to_recalibrate_df_toID_annotation)
  
  slice_to_recalibrate <- df_is[intersect(rownames(df_is), levels(factor(to_recalibrate_df$from_id))), setdiff(colnames(df_is), id_cols) ]
  
  # annotation from -> to
  message(paste("[AP]\t\tMerging data recalibration and source IS"))
  slice_to_recalibrate_annotation <- merge( 
    y = slice_to_recalibrate,
    x = to_recalibrate_df[rownames(slice_to_recalibrate), ],
    by = 0, all.y = TRUE 
  )
  rownames(slice_to_recalibrate_annotation) <- slice_to_recalibrate_annotation$Row.names
  # id_cols_remap <- c("chr", "to_start", "strand", "GeneName", "GeneStrand") # col names of id_cols are derived from the TO_ID -> thus the correct final annotation
  slice_recalibrated <- slice_to_recalibrate_annotation[c(id_cols, setdiff(colnames(df_is), id_cols))]
  # names(slice_recalibrated) <- c(id_cols, colnames(slice_recalibrated[grep(project_id, colnames(slice_recalibrated))]))
  # names(slice_recalibrated) <- gsub("to_start", "integration_locus", colnames(slice_recalibrated))
  
  # get existing IS, that will be fused with the recalibrated IS
  # set A
  slice_recalibrated_to_fuse <- df_is[intersect(rownames(df_is), levels(factor(to_recalibrate_df$to_id))), ]
  # this resulted empty ;) -> so you can change here the rownames in the final row names (recalibrate) and you will not have duplicated values
  # rownames(slice_recalibrated) <- apply(slice_recalibrated[id_cols], 1, function(x) { paste("chr",x[1],":",as.numeric(x[2]),"(",x[3],")_",x[4],"_",x[5], sep="") })
  
  message(paste("[AP]\t\tAggragate results of recalibration to fuse rows"))
  # fuse data alla cazzo, then with aggregate
  slice_recalibrated_mixed <- rbind(slice_recalibrated, slice_recalibrated_to_fuse)
  slice_recalibrated_mixed$chr <- gsub("chr", "", slice_recalibrated_mixed$chr)
  
  # aggregate data (fuse rows)
  slice_recalibrated_mixed_aggregated <- aggregate(slice_recalibrated_mixed[, setdiff(colnames(slice_recalibrated_mixed), id_cols)], by=list(slice_recalibrated_mixed$chr, slice_recalibrated_mixed$integration_locus, slice_recalibrated_mixed$strand, slice_recalibrated_mixed$GeneName, slice_recalibrated_mixed$GeneStrand), FUN=function(x) {sum(x, na.rm = T)})
  # rownames(slice_recalibrated_mixed_aggregated) <- apply(slice_recalibrated_mixed_aggregated, 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
  
  # slice_recalibrated_mixed_aggregated <- aggregate(slice_recalibrated_mixed[, setdiff(colnames(slice_recalibrated_mixed), id_cols)], by=list(slice_recalibrated_mixed$chr, slice_recalibrated_mixed$integration_locus, slice_recalibrated_mixed$strand), FUN=function(x) {sum(x, na.rm = T)})
  # # change name -> according to only _ based
  # rownames(slice_recalibrated_mixed_aggregated) <- apply(slice_recalibrated_mixed_aggregated, 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })
  # # to add the correct annotations, do the join with the source df, and first change rownames
  # rownames(slice_recalibrated_to_fuse) <- apply(slice_recalibrated_to_fuse[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })
  # a <- merge(x = slice_recalibrated_mixed_aggregated, y = slice_recalibrated_to_fuse[id_cols], by = 0, all = TRUE)
  
  # slice_recalibrated_topaste <- a[,c(id_cols, setdiff(colnames(a), c(id_cols, c("Row.names", "Group.1", "Group.2", "Group.3") )))]
  colnames(slice_recalibrated_mixed_aggregated)[1:5] <- id_cols
  slice_recalibrated_topaste <- slice_recalibrated_mixed_aggregated[col_order]
  rownames(slice_recalibrated_topaste) <- apply(slice_recalibrated_topaste[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
  
  # # do some checks
  # venn(list("source df_is" = rownames(df_is), "slice_recalibrated_topaste" = rownames(slice_recalibrated_topaste), "slice_recalibrated_mixed" = rownames(slice_recalibrated_mixed)))
  
  # update slice matrix
  # now merge data for the final df recalibrated
  out_df <- rbind(df_is[!(rownames(df_is) %in% rownames(slice_recalibrated_mixed)), col_order], slice_recalibrated_topaste)
  out_df$integration_locus <- as.numeric(out_df$integration_locus)
  return (out_df)
}

###############################################################
#' @title Merge matrixes
#' 
#' @author Andrea Calabria
#' @details version 1.1 (used in clonal tracking), 2020-09-03
#'
#' @rdname mergeAnnotatedMatrixes
#' @docType methods
#' @aliases recalibrateIS
#'
#' @param df_tomerge vector or input dataframe even with annotation
#' @param id_cols annotation cols
#'
#' @return a df with the same number of columns of the input df (df_is) but at most the same number of rows.
#' @usage TODO
#' @note : 
#'
###############################################################
mergeAnnotatedMatrixes <- function(df_tomerge, id_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand") ) {
  library(reshape2)
  message(paste("[AP]\tMerge matrixes"))
  # strategy: melt and cast
  molten_df <- NULL
  for (mat in df_tomerge) {
    mat[mat==0] <- NA
    this_molten_df <- melt(data = mat, id.vars = id_cols, na.rm = T, variable.name = "Sample")
    if (length(molten_df)>0) {
      molten_df <- rbind(molten_df, this_molten_df)
    } else {
      molten_df <- this_molten_df
    }
  }
  cast_df <- dcast(data = molten_df, chr + integration_locus + strand + GeneName + GeneStrand ~ Sample, value.var = "value", fun.aggregate = mean)
  cast_df[is.na(cast_df)] <- NA
  return (cast_df)
}


#' ###############################################################
#' #' @title Recalibrate IS using a BED map and Merge df
#' #' 
#' #' @author Andrea Calabria
#' #' @details version 1.0 (used in clonal tracking), 2018-10-23
#' #'
#' #' @rdname recalibrateIS_and_addNewIS
#' #' @docType methods
#' #' @aliases recalibrateIS_and_addNewIS
#' #'
#' #' @param df_is an input dataframe even with annotation
#' #' @param recalibration_map a df importing the BED file annotation from reference, thus having the following structure of columns: "from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance"
#' #' @param recalibrate_IS_distance integer of the distance in bp to distinguish tow ajacent IS. Default = 7, that is <=7.
#' #' @param id_cols a vector of column names for annotation. Default: c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
#' #'
#' #' @return a df with the same number of columns of the input df (df_is) but at most the same number of rows.
#' #' @usage TODO
#' #' @note This function works to recalibrate IS with distance <threshold and adds the IS with distance >threshold
#' #' 
#' #' Working test: 
#' #' 
#' #' the following code is for testing the function, assuming that in the folder "data" you have the appropriate files.
#' #' t.expected_recalibratedIS <- read.csv("source/test/test.matrix.foo.expected_after_recalibration.csv", check.names = F, header = T, sep = "\t", row.names = 1) # load the recalibrated expected matrix
#' #' id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")
#' #' t.in <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/data/test.matrix.foo.v2.csv", check.names = F, header = T, sep = "\t", row.names = 1)
#' #' rownames(t.in) <- apply(t.in[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
#' #' t.rmap <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/data/test.matrix.foo.v2.withRef.bed", check.names = F, header = F, sep = "\t")
#' #' rownames(t.rmap) <- t.rmap$V4
#' #' names(t.rmap) <- c("from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance")
#' #' t.recalibratedIS <- recalibrateIS(df_is = t.in, recalibration_map = t.rmap)
#' #' 
#' #'
#' ###############################################################
#' recalibrate_and_addNewIS <- function(df_is, recalibration_map, recalibrate_IS_distance = 7, id_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand") ) {
#'   library(gplots)
#'   ##### -------- RECALIBRATION ------- #####
#'   message(paste("[AP]\tRecalibration"))
#'   # strategy:
#'   # - slice df with only ID to remap / recalibrate (the "from" field)
#'   # - correct IS (apply the ID of the "to" field)
#'   # - extract the "to" fils IS ID
#'   # - fuse IS with (potential) existing IDs by summing elements (cells)
#'   # - remove IS ID from-to from dfallp
#'   # - add the new aggregated/fused matrix to dfallp
#'   
#'   # load recalibration matrix
#'   message(paste("[AP]\t\tLoading BED file and recalibration data, from file"))
#'   # define a col order
#'   col_order <- colnames(df_is)
#'   # recalibration_map <- read.csv(file_recalibration_map_bed, header=F, fill=T, sep='\t', check.names = FALSE)
#'   # rownames(recalibration_map) <- recalibration_map$V4
#'   names(recalibration_map) <- c("from_chr", "from_start", "from_end", "from_id", "from_score", "from_strand", "to_chr", "to_start", "to_end", "to_id", "to_score", "to_strand", "distance")
#'   to_recalibrate_df <- recalibration_map[which(recalibration_map$from_start != recalibration_map$to_start & abs(recalibration_map$distance) <= recalibrate_IS_distance ),]
#'   # do the annotation from ID
#'   to_recalibrate_df_toID_annotation <- as.data.frame(t(as.data.frame(lapply(strsplit(as.character(to_recalibrate_df$to_id), "_", fixed = TRUE), function(x) {c(x)}))))
#'   names(to_recalibrate_df_toID_annotation) <- id_cols
#'   to_recalibrate_df <- cbind(to_recalibrate_df, to_recalibrate_df_toID_annotation)
#'   
#'   slice_to_recalibrate <- df_is[intersect(rownames(df_is), levels(factor(to_recalibrate_df$from_id))), setdiff(colnames(df_is), id_cols) ]
#'   
#'   # annotation from -> to
#'   message(paste("[AP]\t\tMerging data recalibration and source IS"))
#'   slice_to_recalibrate_annotation <- merge( 
#'     y = slice_to_recalibrate,
#'     x = to_recalibrate_df[rownames(slice_to_recalibrate), ],
#'     by = 0, all.y = TRUE 
#'   )
#'   rownames(slice_to_recalibrate_annotation) <- slice_to_recalibrate_annotation$Row.names
#'   # id_cols_remap <- c("chr", "to_start", "strand", "GeneName", "GeneStrand") # col names of id_cols are derived from the TO_ID -> thus the correct final annotation
#'   slice_recalibrated <- slice_to_recalibrate_annotation[c(id_cols, setdiff(colnames(df_is), id_cols))]
#'   # names(slice_recalibrated) <- c(id_cols, colnames(slice_recalibrated[grep(project_id, colnames(slice_recalibrated))]))
#'   # names(slice_recalibrated) <- gsub("to_start", "integration_locus", colnames(slice_recalibrated))
#'   
#'   # get existing IS, that will be fused with the recalibrated IS
#'   # set A
#'   slice_recalibrated_to_fuse <- df_is[intersect(rownames(df_is), levels(factor(to_recalibrate_df$to_id))), ]
#'   # this resulted empty ;) -> so you can change here the rownames in the final row names (recalibrate) and you will not have duplicated values
#'   # rownames(slice_recalibrated) <- apply(slice_recalibrated[id_cols], 1, function(x) { paste("chr",x[1],":",as.numeric(x[2]),"(",x[3],")_",x[4],"_",x[5], sep="") })
#'   
#'   message(paste("[AP]\t\tAggragate results of recalibration to fuse rows"))
#'   # fuse data alla cazzo, then with aggregate
#'   slice_recalibrated_mixed <- rbind(slice_recalibrated, slice_recalibrated_to_fuse)
#'   slice_recalibrated_mixed$chr <- gsub("chr", "", slice_recalibrated_mixed$chr)
#'   
#'   # aggregate data (fuse rows)
#'   slice_recalibrated_mixed_aggregated <- aggregate(slice_recalibrated_mixed[, setdiff(colnames(slice_recalibrated_mixed), id_cols)], by=list(slice_recalibrated_mixed$chr, slice_recalibrated_mixed$integration_locus, slice_recalibrated_mixed$strand, slice_recalibrated_mixed$GeneName, slice_recalibrated_mixed$GeneStrand), FUN=function(x) {sum(x, na.rm = T)})
#'   # rownames(slice_recalibrated_mixed_aggregated) <- apply(slice_recalibrated_mixed_aggregated, 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
#'   
#'   # slice_recalibrated_mixed_aggregated <- aggregate(slice_recalibrated_mixed[, setdiff(colnames(slice_recalibrated_mixed), id_cols)], by=list(slice_recalibrated_mixed$chr, slice_recalibrated_mixed$integration_locus, slice_recalibrated_mixed$strand), FUN=function(x) {sum(x, na.rm = T)})
#'   # # change name -> according to only _ based
#'   # rownames(slice_recalibrated_mixed_aggregated) <- apply(slice_recalibrated_mixed_aggregated, 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })
#'   # # to add the correct annotations, do the join with the source df, and first change rownames
#'   # rownames(slice_recalibrated_to_fuse) <- apply(slice_recalibrated_to_fuse[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3], sep="") })
#'   # a <- merge(x = slice_recalibrated_mixed_aggregated, y = slice_recalibrated_to_fuse[id_cols], by = 0, all = TRUE)
#'   
#'   # slice_recalibrated_topaste <- a[,c(id_cols, setdiff(colnames(a), c(id_cols, c("Row.names", "Group.1", "Group.2", "Group.3") )))]
#'   colnames(slice_recalibrated_mixed_aggregated)[1:5] <- id_cols
#'   slice_recalibrated_topaste <- slice_recalibrated_mixed_aggregated[col_order]
#'   rownames(slice_recalibrated_topaste) <- apply(slice_recalibrated_topaste[id_cols], 1, function(x) { paste("chr",x[1],"_",as.numeric(x[2]),"_",x[3],"_",x[4],"_",x[5], sep="") })
#'   
#'   # # do some checks
#'   # venn(list("source df_is" = rownames(df_is), "slice_recalibrated_topaste" = rownames(slice_recalibrated_topaste), "slice_recalibrated_mixed" = rownames(slice_recalibrated_mixed)))
#'   
#'   # update slice matrix
#'   # now merge data for the final df recalibrated
#'   out_df <- rbind(df_is[!(rownames(df_is) %in% rownames(slice_recalibrated_mixed)), col_order], slice_recalibrated_topaste)
#'   out_df$integration_locus <- as.numeric(out_df$integration_locus)
#'   
#'   # now add the IS outside the boundaries
#'   to_add_df <- recalibration_map[which( (recalibration_map$from_start != recalibration_map$to_start & abs(recalibration_map$distance) > recalibrate_IS_distance) ),]
#'   out_df <- rbind(out_df, df_is[to_add_df$from_id,])
#'   
#'   return (out_df)
#' }


hsc_schnabel <- function(c1, c2, s) {
  message(paste("\tHSC Estimate, Petersen-Schnabel", "\n\t\tC1:", (length(c1) + 1), "\n\t\tC2:", (length(c2) + 1), "\n\t\tS:", (length(s) + 1)))
  n <- ( ((length(c1) + 1) * (length(c2) + 1)) / (length(s) + 1) ) - 1
  message(paste("\t\tEstimate:", round(n)))
  return (n)
}

lm_eqn <- function(df){
  # prototype example: p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}


createBEDdfFromMatrix <- function(df, matrix_format, id_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand")){
  # input: new format matrix as df with rownames as genomic ID: chrX_12345_+
  # output: df of a bed file format. NB:
  #   - end = start +1
  #   - name = rowname(df)
  # the ID will be always in the form of the NEW matrix: chrX_12345_+
  
  # define here some useful vars
  error <- NULL
  bed_fields <- c("chr", "start", "end", "id", "score", "strand") # for the output file format (sorted in this way)
  # id_cols <- c("chr", "integration_locus", "strand", "GeneName", "GeneStrand") # for the old matrix
  
  message(paste("[AP]\tCreat BED file from matrix"))
  if (matrix_format == "new") {
    message(paste("[AP]\t\tNew file format"))
    # decompose the genomic ID in its components and create the BED df
    row_ids <- rownames(df)
    if (length(grep("^chr", row_ids[1])) == 0 ) { # your IDs do not start with "chr"
      row_ids <- gsub("^", "chr", row_ids)
      }
    genomic_df <- data.frame(t(as.data.frame(strsplit(row_ids, "_"))))
    rownames(genomic_df) <- row_ids
    names(genomic_df) <- c("chr", "start", "strand")
    # now create a real bed df
    genomic_df <- cbind(genomic_df, 
                        data.frame("end" = apply(genomic_df, 1, function(x) {as.numeric(x[2])+1}),
                                   "id" = rownames(df),
                                   "score" = rep(100, nrow(df))
                                   )
                        )
    return (genomic_df[bed_fields])
  } else if (matrix_format == "old") {
    message(paste("[AP]\t\tOld file format"))
    # decompose the genomic ID in its components and create the BED df
    genomic_df <- df[id_cols[1:3]]
    if ( length(grep("^chr", genomic_df[1,"chr"])) == 0 ) { # your IDs do not start with "chr"
      genomic_df$chr <- gsub("^", "chr", genomic_df$chr)
    }
    names(genomic_df) <- c("chr", "start", "strand")
    # now create a real bed df
    genomic_df <- cbind(genomic_df, 
                        data.frame("end" = apply(genomic_df, 1, function(x) {as.numeric(x[2])+1}),
                                   "id" = apply(genomic_df, 1, function(x) { paste(x[1], as.numeric(x[2]), x[3], sep="_") }),
                                   "score" = rep(100, nrow(df))
                        )
    )
    return (genomic_df[bed_fields])
  } else {
    message(paste("[AP]\t\tERROR: matrix_format should be {new, old}"))
    error <- 1
    return (NULL)
  }
}

###############################################################
#' @title Calcolate a graph structure G(nodes, edges) from df of IS
#' 
#' @author Andrea Calabria
#' @details version 1.0 (used in hematopoietic tree graphs)
#'
#' @rdname getGraphFromMatrixes
#' @docType methods
#' @aliases getGraphFromMatrixes
#'
#' @param df an input dataframe patient_iss[setdiff(colnames(patient_iss), id_cols)]
#' @param stats_col_selection_chr a vector of columns: c("id", "patient_id", "cell_marker", "source", "cell_type", "pcr_method", "groupby_tissue_fumonths")
#' @param onRatio a bolean value (default TRUE), the output will be computed as ratio of sharing (default is abslute numbers: a triangular matrix) that is not a triangular matrix
#' @param colorBy a string of a choice among {cell_type, cell_marker} annotation
#' @param removeDiagonal a bolean value (default TRUE) to remove diagonal from input df
#' @param doTriangular a bolean value (default TRUE) to do the triangular version of the input df
#'
#' @return a list of 3 df an object: nodes, links, adjacency mattrix
#' @usage patient_graph <- getGraphFromMatrixes(df_iss = patient_iss[setdiff(colnames(patient_iss), id_cols)], df_stats = patient_stats, onRatio = TRUE, removeDiagonal = TRUE, colorBy = "cell_marker")
#'
###############################################################
# df: patient_iss[setdiff(colnames(patient_iss), id_cols)]
# stats_col_selection_chr: c("id", "patient_id", "cell_marker", "source", "cell_type", "pcr_method", "groupby_tissue_fumonths")
# stats_col_selection_num: c("followup", "followup_months", "nIS")
# onRatio: the output will be computed as ratio of sharing (default is abslute numbers: a triangular matrix) that is not a triangular matrix
# colorBy: select the coloring option by one of the available annotation: cell_type, cell_marker
###############################################################
getGraphFromMatrixes <- function(df_iss = patient_iss, df_stats = patient_stats, stats_col_selection_chr = c("id", "patient_id", "cell_marker", "source", "cell_type", "pcr_method", "groupby_tissue_fumonths"), stats_col_selection_num = c("followup", "followup_months", "nIS"), onRatio = TRUE, colorBy = "cell_type", removeDiagonal = TRUE, doTriangular = FALSE ) {
  
  # load libraries
  # library(visNetwork)
  library(reshape2)
  library(RColorBrewer)
  library(stringr)
  
  message(paste("[AP]\tCreating Graph from input dataframe"))
  
  ### -------- NODES ------- ###
  # select nodes
  message(paste("[AP]\t-> Nodes"))
  df_iss.nodes <- cbind(data.frame(lapply(df_stats[stats_col_selection_chr], as.character), stringsAsFactors=FALSE), df_stats[stats_col_selection_num])
  # add some general info for visNetwork
  # color info
  colorBy_numeric <- paste(colorBy, "_numeric", sep = '')
  if (colorBy == "cell_marker") {
    df_iss.nodes <- merge(x = df_iss.nodes, y = data.frame("cell_marker" = unique(df_iss.nodes$cell_marker), "cell_marker_numeric" = seq(1, length(unique(df_iss.nodes$cell_marker)))), by = c("cell_marker"), all.x = TRUE)
  } else {
    df_iss.nodes <- merge(x = df_iss.nodes, y = data.frame("cell_type" = unique(df_iss.nodes$cell_type), "cell_type_numeric" = seq(1, length(unique(df_iss.nodes$cell_type)))), by = c("cell_type"), all.x = TRUE)
  }
  df_iss.nodes$shape <- "dot"  
  df_iss.nodes$shadow <- TRUE # Nodes will drop shadow
  df_iss.nodes$title <- df_iss.nodes$id # Text on click
  df_iss.nodes$label <- df_iss.nodes$cell_marker # Node label
  df_iss.nodes$size <- round((df_iss.nodes$nIS)/40) # Node size
  df_iss.nodes$borderWidth <- 2 # Node border width
  # color_palette_spectral <- brewer.pal(nrow(unique(df_iss.nodes[colorBy_numeric])), "Spectral")
  color_palette_spectral <- brewer.pal(nrow(unique(df_iss.nodes[colorBy_numeric])), "Paired")
  color_palette_spectral_names <- colors()[ apply(col2rgb(color_palette_spectral), 2, function(z) {which.min( sapply(colors(), function(x) col.dist(inp=z, comp=x) ) )} ) ]
  df_iss.nodes$color.background <- color_palette_spectral_names[df_iss.nodes[, colorBy_numeric]]
  df_iss.nodes$color.border <- "black"
  df_iss.nodes$color.highlight.background <- "orange"
  df_iss.nodes$color.highlight.border <- "darkred"
  rownames(df_iss.nodes) <- apply(df_iss.nodes[c("cell_marker", "source", "followup_months")], 1, function(x){paste(x[1], x[2], str_pad(as.numeric(x[3]), 2, pad = 0), sep = "_")})
  
  ### -------- EDGES ------- ###
  message(paste("[AP]\t-> Edges"))
  # on Ratio?
  df_iss.links <- NULL
  df_iss.adjmat <- NULL
  if (onRatio) {
    message(paste("[AP]\t-> (weights on sharing ratio)"))
    # do all sharing among all IS and create the ratio and create a general graph structure
    df_iss.adjmat <- getSharedISratio(df = df_iss[setdiff(colnames(df_iss), id_cols)], compact = T, left_to_rigth_reading_output = FALSE)
  } else { # default on IS absolute shared numbers
    # do all sharing among all IS and create a general graph structure
    df_iss.adjmat <- getSharedISnumber(df = df_iss[setdiff(colnames(df_iss), id_cols)], compact = T, left_to_rigth_reading_output = FALSE)
  } # if (onRatio)
  # convert in Triangular?
  if (doTriangular) {
    df_iss.adjmat[lower.tri(x = df_iss.adjmat, diag = FALSE)] <- NA # do triangular
  } # if (doTriangular)
  # remove diagonal?
  if (removeDiagonal) {
    diag(df_iss.adjmat) <- NA
  } # if (removeDiagonal)
  # now melt data
  df_iss.adjmat_melt <- melt(as.matrix(df_iss.adjmat))
  names(df_iss.adjmat_melt) <- c("from", "to", "weight")
  df_iss.adjmat_melt <- df_iss.adjmat_melt[which(!is.na(df_iss.adjmat_melt$weight) & df_iss.adjmat_melt$weight > 0),]
  df_iss.adjmat_melt$type <- rep("hyperlink", nrow(df_iss.adjmat_melt)) # add the type column
  df_iss.links <- df_iss.adjmat_melt
  # do aesthetic
  df_iss.links$width <- 0.5 + df_iss.links$weight*10 # line width
  df_iss.links$color <- "gray"    # line color  
  df_iss.links$arrows <- "to" # arrows: 'from', 'to', or 'middle'
  df_iss.links$smooth <- TRUE    # should the edges be curved?
  df_iss.links$shadow <- FALSE    # edge shadow
  df_iss.links$label <- as.character(format(round(df_iss.links$weight, 2), nsmall = 3)) # edge label
  
  # what to return, a graph N, L
  return( list(df_iss.nodes, df_iss.links, df_iss.adjmat) )
}
###############################################################



###############################################################
#' @title Get a binary matrix of IS by metadata (group by, fixed position)
#' 
#' @author Andrea Calabria
#' @details version 1.0 (used in hematopoietic tree graphs binary vectors)
#'
#' @rdname binarizeColumnsByMetadata
#' @docType methods
#' @aliases binarizeColumnsByMetadata
#'
#' @param df an input dataframe patient_iss.
#' @param metadata_df a df of metadata such as the association file.
#' @param select_distinct_field a string for the column name (aka field) to search for all distinct values in the metadata df. The list of found values will create the vector of columns to search for. This means that each value found in this list will be statically ordered and searched as part of the column names. No multiple rows in output values are admitted -> otherwise you are required to aggregate first! Typical example is selecting cell marker.
#' @param groupby_field a string for the column name (aka field) to group by the rows in the metadata df. The idea is the replica of the SQL format, such that: SELECT distinct_field FROM metadata_df WHERE 1 GROUP BY groupby_field. Typical example is grouping by time point.
#' @param select_field a vector of values cotaining the selection to apply in the distinct_field clause. This is useful to filter data and not view all results. Example: c("CD15", "CD14")
#' @param starting_data_col_index an integer value, the index of the first column to analyze (1..n)
#' @param number_of_last_cols_to_remove an integer value (0..n) of number of last columns to remove from the analysis
#'
#' @return ??? 
#' @usage ???
#' @description This function is aimed at (1) acquiring metadata and completing the input matrix with missing observations (given N values obtained from distinct_field, and M values obtained from groupby_field, the returned matrix will be composed by NxM columns, <= number of real input columns, excluded the cols to omit), (2) converting the df in binary values, (3) acquiring the binary representation of the block of the N values (for each of the M groups) in a string (binary representation, or flag) keeping inaltered the order of the N values among the M groups, (4) producing for each row the integer conversion of the binary string and returning both the df of the M groups as columns with all input rows and in the cells the binary representation and its integer conversion. Moreover the ourput will provide the map df of the N values, their order and their integer number associated (flag).
#'
###############################################################
binarizeColumnsByMetadata <-function(df, 
                                     metadata_df, 
                                     groupby_field, 
                                     select_distinct_field = c("cell_marker"), 
                                     sep = "_", 
                                     select_field = "*", 
                                     starting_data_col_index = 2, 
                                     number_of_last_cols_to_remove = 0) {
  # idea is the replica of the SQL format, such that: SELECT distinct_field FROM metadata_df WHERE 1 GROUP BY groupby_field
  require(sqldf)
  require(dplyr)
  
  # (1) acquiring metadata and completing the input matrix with missing observations (given N values obtained from distinct_field, and M values obtained from groupby_field, the returned matrix will be composed by NxM columns, <= number of real input columns, excluded the cols to omit), 
  # (2) converting the df in binary values, 
  # (3) acquiring the binary representation of the block of the N values (for each of the M groups) in a string (binary representation, or flag) keeping inaltered the order of the N values among the M groups, 
  # (4) producing for each row the integer conversion of the binary string and returning both the df of the M groups as columns with all input rows and in the cells the binary representation and its integer conversion. Moreover the ourput will provide the map df of the N values, their order and their integer number associated (flag)
 
  
  # given the input df, look up the column  names and aggregate by name
  # for each col name, slice df and apply function
  colID_df <- data.frame("colID" = colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]) # get column names of df into a new df
  rownames(colID_df) <- colnames(df)[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)]
  colmetadata_df <- merge(colID_df, metadata_df, by=0, all.x = TRUE) # get the output merge df of annotations with all x metadata
  rownames(colmetadata_df) <- colmetadata_df$Row.names
  
  # # get names from the key field
  # colmetadata_df_outnames <- colmetadata_df[,c(key_field)]
  # 
  # # produce a warning IF two columns share the same pattern matching (prefix or suffix)
  # if (length(grep(key_field, colnames(colmetadata_df)))>1) {
  #   message(paste("[AP]\tWARNING: in your metadata file you have two columns with a very similar name, this is not allowed. I.e.: Group and TestGroup."))
  # }
  
  # get the groups
  groupby_field_allvalues <- levels(factor(colmetadata_df[,groupby_field]))
  distinct_field_allvalues <- levels(factor(colmetadata_df[,distinct_field]))
  
  # output_cols <- NULL
  # for (distinct_field_value in distinct_field_allvalues) { output_cols <- c(output_cols, paste(distinct_field_value, groupby_field_allvalues, sep = key_field_colseparator))}
  
  #message(paste("[AP]\tThese will be the output/aggregated columns:", output_cols))
  out_df <- NULL
  for (groupby_field_value in groupby_field_allvalues) {
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(groupby_field, colnames(colmetadata_df))] == groupby_field_value), c("colID")])  
  }
  
  
  for (k in output_cols) {
    # get col IDs from metadata
    #metadata_df[which(metadata_df[grep(key_field, colnames(metadata_df))] == k), c("colID")]
    message(paste("[AP]\tProcessing column(s)", k))
    slice_cols <- as.character(colmetadata_df[which(colmetadata_df[grep(key_field, colnames(colmetadata_df))] == k), c("colID")])
    if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        out_df <- df[slice_cols]
        out_df[is.na(out_df)] <- 0
        names(out_df) <- k
      } else {
        out_df <- as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)}))
        names(out_df) <- k
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    } else { # if (length(out_df) == 0) { # if this is the first loop
      if (length(slice_cols) == 1) { # if slice_cols contains only 1 value
        actual_colnames <- colnames(out_df)
        this_slice <- df[, slice_cols]
        this_slice[is.na(this_slice)] <- 0
        out_df <- cbind(out_df, this_slice )
        names(out_df) <- c(actual_colnames, k)
      } else {
        actual_colnames <- colnames(out_df)
        out_df <- cbind(out_df, as.data.frame(apply(df[, slice_cols], 1, function(x) {sum(x, na.rm = TRUE)})) )
        names(out_df) <- c(actual_colnames, k)
      } # if (len(slice_cols) == 1) { # if slice_cols contains only 1 value
    }
  } # else of if (length(out_df) == 0) { # if this is the first loop
  # return df
  return (out_df)
  
}

filterOutliersByPoolFragments <- function(metadata_df, 
                                          key_field_quantification = "QuantificationSum", 
                                          key_field_poolID = "pool_id", 
                                          per_pool_test = TRUE, 
                                          # select_key_field = "pcr_method", 
                                          # select_value_to_exclude = c("SLiM"),
                                          min_samples_per_pool = 30, 
                                          p_value_threshold = 0.05) {
  require(nortest)
  metadata_df_slice <- metadata_df[which(metadata_df[key_field_quantification]>0),]
  outlist_samplestofilter <- NULL
  df_outlist_samplestofilter <- NULL
  # metadata_df_slice$QuantificationSum_log2 <- log2(metadata_df_slice[key_field_quantification])
  if (per_pool_test) { # if you want to perform the analysis by pool:
    available_pools <- table(metadata_df_slice[key_field_poolID]) # levels(factor(metadata_df_slice[key_field_poolID]))
    selected_available_pools <- available_pools[available_pools >= min_samples_per_pool]
    for (pool_id in labels(selected_available_pools)[[1]]) {
      message(paste("[AP]\tProcessing pool", pool_id))
      retrieved_data_df <- metadata_df_slice[which(metadata_df_slice[key_field_poolID] == pool_id),] # subset(metadata_df_slice, PoolID == pool_id)
      if (nrow(retrieved_data_df) > 0) {
        message(paste("[AP]\t\tFound, go ahead with", pool_id))
        # # is the sample normal distributed
        # shapiro.test(retrieved_data_df$QuantificationSum_log2)
        scaled_data <- scale(retrieved_data_df[key_field_quantification])
        is_normal <- shapiro.test(scaled_data)
        lillie.test(scaled_data)
        shapiro.test(scaled_data)
        
        log2_data <- log2(retrieved_data_df[key_field_quantification])
        is_normal_log2 <- shapiro.test(log2_data[,key_field_quantification])
        lillie.test(log2_data[,key_field_quantification])
        # agostino.test(log2_data[,key_field_quantification])
        
        # if (is_normal$p.value >= p_value_threshold) {
        #   transformed_data <- scaled_data
        #   message(paste("[AP]\t\t!!Luky!!! The pool", pool_id, "follows a Normal distribution!"))
        #   avg <- mean(retrieved_data_df[,key_field_quantification])
        #   std <- sd(retrieved_data_df[,key_field_quantification])
        #   retrieved_data_df$zscore <- (retrieved_data_df[,key_field_quantification] - avg)/std
        #   count <- nrow(retrieved_data_df)
        #   retrieved_data_df$tstudent <- sqrt(((count*(count-2)*(retrieved_data_df$zscore)^2)/((count-1)^2-count*(retrieved_data_df$zscore)^2)))
        #   retrieved_data_df$tdist <- dt(retrieved_data_df$tstudent, df = count-2)
        # } else {
        #   message(paste("[AP]\t\tThe pool", pool_id, "does not follow a Normal distribution, skip it for the test of outliers."))
        # }
        if (is_normal_log2$p.value >= p_value_threshold) {
          transformed_data <- log2_data[,key_field_quantification]
          message(paste("[AP]\t\t!!Luky!!! The pool", pool_id, "follows a Normal distribution!"))
          avg <- mean(transformed_data)
          std <- sd(transformed_data)
          retrieved_data_df$zscore <- (transformed_data - avg)/std
          count <- nrow(retrieved_data_df)
          retrieved_data_df$tstudent <- sqrt(((count*(count-2)*(retrieved_data_df$zscore)^2)/((count-1)^2-count*(retrieved_data_df$zscore)^2)))
          retrieved_data_df$tdist <- dt(retrieved_data_df$tstudent, df = count-2)
          # select data to filter out
          outlist_samplestofilter <- c(outlist_samplestofilter, rownames(retrieved_data_df[which(retrieved_data_df$tdist < p_value_threshold & retrieved_data_df$zscore < 0),]))
          if (length(df_outlist_samplestofilter) == 0){ # first loop
            df_outlist_samplestofilter <- retrieved_data_df[which(retrieved_data_df$tdist < p_value_threshold & retrieved_data_df$zscore < 0),]
          } else {
            df_outlist_samplestofilter <- rbind(df_outlist_samplestofilter, retrieved_data_df[which(retrieved_data_df$tdist < p_value_threshold & retrieved_data_df$zscore < 0),])
          }
        } else {
          message(paste("[AP]\t\tThe pool", pool_id, "does not follow a Normal distribution, skip it for the test of outliers."))
        }
        
        # write.csv(retrieved_data_df, file = paste("workbench Sept 2017/", pool_id, ".csv", sep = ''), row.names = FALSE)
        
      } else
        message(paste("[AP]\t\tERROR: Not Found", pool_id))
    } 
  } else { # if (per_pool_test)
    message(paste("[AP]\t\tElse branch"))
  }
  
  return (df_outlist_samplestofilter)
  
}

# just flag
flagOutliersByPoolFragments <- function(metadata_df, 
                                        key_field_quantification = "QuantificationSum", 
                                        key_field_poolID = "pool_id", 
                                        transform_log2 = TRUE,
                                        per_pool_test = TRUE, 
                                        # select_key_field = "pcr_method", 
                                        # select_value_to_exclude = c("SLiM"),
                                        min_samples_per_pool = 5,
                                        run_normality_test = F,
                                        normality_p_value_threshold = 0.1,
                                        outlier_p_value_threshold = 0.05) {
  require(nortest)
  message(paste("[AP]\tRemoving missing values (NA) in the field", key_field_quantification))
  metadata_df_slice <- metadata_df[which(!(is.na(metadata_df[key_field_quantification])) ),]
  outlist_samplestofilter <- NULL
  df_outlist_samplestofilter <- NULL
  # metadata_df_slice$QuantificationSum_log2 <- log2(metadata_df_slice[key_field_quantification])
  if (per_pool_test) { # if you want to perform the analysis by pool:
    available_pools <- table(metadata_df_slice[key_field_poolID]) # levels(factor(metadata_df_slice[key_field_poolID]))
    selected_available_pools <- available_pools[available_pools >= min_samples_per_pool]
    for (pool_id in labels(selected_available_pools)[[1]]) {
      message(paste("[AP]\tProcessing pool", pool_id))
      retrieved_data_df <- metadata_df_slice[which(metadata_df_slice[key_field_poolID] == pool_id),] # subset(metadata_df_slice, PoolID == pool_id)
      # init a variable (transformed_data) that only slices useful data
      transformed_data <- retrieved_data_df[key_field_quantification]
      if (transform_log2) {
        message(paste("[AP]\t\tAsking for log2 transformation requires to remove values <= 0 in", pool_id))
        # first you need to remove 0 values x<=0
        retrieved_data_df <- metadata_df_slice[which(metadata_df_slice[key_field_poolID] == pool_id & 
                                                       metadata_df_slice[key_field_quantification] > 0 &
                                                       !(is.na(metadata_df_slice[key_field_quantification]))
        ),] # subset(metadata_df_slice, PoolID == pool_id)
        # transformed_data <- retrieved_data_df[key_field_quantification]
        transformed_data <- log2(retrieved_data_df[key_field_quantification])
      }
      # now evaluate the pool
      if (nrow(retrieved_data_df) > 0) {
        message(paste("[AP]\t\tFound, go ahead with", pool_id))
        if (run_normality_test) {
          # log2_data <- log2(retrieved_data_df[key_field_quantification])
          is_normal <- shapiro.test(transformed_data[,key_field_quantification])
          # lillie.test(log2_data[,key_field_quantification])
          # agostino.test(log2_data[,key_field_quantification])
          if (is_normal$p.value >= normality_p_value_threshold) {
            fields_to_return <- c("transformed_data", "normal_distribution", "zscore", "tstudent", "tdist", "to_remove") # init the fields that will be used along this process
            fields_to_process <- c("transformed_data", "zscore", "tstudent", "tdist") # init the fields that will be used along this process
            retrieved_data_df$transformed_data <- transformed_data[,key_field_quantification] # transformed_data <- log2_data[,key_field_quantification]
            retrieved_data_df$normal_distribution <- TRUE
            message(paste("[AP]\t\t!!Luky!!! The pool", pool_id, "follows a Normal distribution!"))
            # avg <- mean(log2_data[,key_field_quantification])
            # std <- sd(log2_data[,key_field_quantification])
            # retrieved_data_df$zscore <- (retrieved_data_df$log2data - avg)/std
            retrieved_data_df$zscore <- scale(retrieved_data_df$transformed_data)
            count <- nrow(retrieved_data_df)
            retrieved_data_df$tstudent <- sqrt(((count*(count-2)*(retrieved_data_df$zscore)^2)/((count-1)^2-count*(retrieved_data_df$zscore)^2)))
            retrieved_data_df$tdist <- dt(retrieved_data_df$tstudent, df = count-2)
            retrieved_data_df$to_remove <- apply(retrieved_data_df[fields_to_process], 1, function(x) {ifelse((x[4] < outlier_p_value_threshold & x[2] < 0), TRUE, FALSE)}) # [which(retrieved_data_df$tdist < outlier_p_value_threshold & retrieved_data_df$zscore < 0),] 
            # append results
            if (length(df_outlist_samplestofilter) == 0){ # first loop
              df_outlist_samplestofilter <- retrieved_data_df[fields_to_return]
            } else {
              df_outlist_samplestofilter <- rbind(df_outlist_samplestofilter, retrieved_data_df[fields_to_return])
            }
          } else {
            message(paste("[AP]\t\tThe pool", pool_id, "does not follow a Normal distribution, skip it for the test of outliers."))
            # retrieved_data_df$normal_distribution <- FALSE
            # fields_to_return <- c("transformed_data", "normal_distribution") # init the fields that will be used along this process
            # # append results
            # if (length(df_outlist_samplestofilter) == 0){ # first loop
            #   df_outlist_samplestofilter <- retrieved_data_df[fields_to_return]
            # } else {
            #   df_outlist_samplestofilter <- rbind(df_outlist_samplestofilter, retrieved_data_df[fields_to_return])
            # }
          } # if (is_normal$p.value >= normality_p_value_threshold)
        } else {
          
          fields_to_return <- c("transformed_data", "normal_distribution", "zscore", "tstudent", "tdist", "to_remove") # init the fields that will be used along this process
          fields_to_process <- c("transformed_data", "zscore", "tstudent", "tdist") # init the fields that will be used along this process
          retrieved_data_df$transformed_data <- transformed_data[,key_field_quantification] # transformed_data <- log2_data[,key_field_quantification]
          retrieved_data_df$normal_distribution <- TRUE
          message(paste("[AP]\t\t!!Luky!!! The pool", pool_id, "follows a Normal distribution!"))
          retrieved_data_df$zscore <- scale(retrieved_data_df$transformed_data)
          count <- nrow(retrieved_data_df)
          retrieved_data_df$tstudent <- sqrt(((count*(count-2)*(retrieved_data_df$zscore)^2)/((count-1)^2-count*(retrieved_data_df$zscore)^2)))
          retrieved_data_df$tdist <- dt(retrieved_data_df$tstudent, df = count-2)
          retrieved_data_df$to_remove <- apply(retrieved_data_df[fields_to_process], 1, function(x) {ifelse((x[4] < outlier_p_value_threshold & x[2] < 0), TRUE, FALSE)}) # [which(retrieved_data_df$tdist < outlier_p_value_threshold & retrieved_data_df$zscore < 0),] 
          # append results
          if (length(df_outlist_samplestofilter) == 0){ # first loop
            df_outlist_samplestofilter <- retrieved_data_df[fields_to_return]
          } else {
            df_outlist_samplestofilter <- rbind(df_outlist_samplestofilter, retrieved_data_df[fields_to_return])
          }
        } # if (run_normality_test) 
        # write.csv(retrieved_data_df, file = paste("workbench Sept 2017/", pool_id, ".csv", sep = ''), row.names = FALSE)
        
      } else {
        message(paste("[AP]\t\tERROR: Not rows remaining in pool", pool_id))
      } # if (nrow(retrieved_data_df) > 0)
    } # for (pool_id in labels(selected_available_pools)[[1]])
  } else { # if (per_pool_test)
    message(paste("[AP]\t\tElse branch"))
    retrieved_data_df <- metadata_df_slice[which(  metadata_df_slice[key_field_quantification] > 0 &
                                                     !(is.na(metadata_df_slice[key_field_quantification]))
                          ),] # subset(metadata_df_slice, PoolID == pool_id)
    # init a variable (transformed_data) that only slices useful data
    transformed_data <- metadata_df_slice[key_field_quantification]
    if (transform_log2) {
      message(paste("[AP]\t\tAsking for log2 transformation requires to remove values <= 0 in", pool_id))
      # first you need to remove 0 values x<=0
      retrieved_data_df <- metadata_df_slice[which(  metadata_df_slice[key_field_quantification] > 0 &
                                                     !(is.na(metadata_df_slice[key_field_quantification]))
                            ),] # subset(metadata_df_slice, PoolID == pool_id)
      # transformed_data <- retrieved_data_df[key_field_quantification]
      transformed_data <- log2(retrieved_data_df[key_field_quantification])
    }
    # now evaluate the pool
    if (nrow(retrieved_data_df) > 0) {
      # log2_data <- log2(retrieved_data_df[key_field_quantification])
      is_normal <- shapiro.test(transformed_data[,key_field_quantification])
      # lillie.test(log2_data[,key_field_quantification])
      # agostino.test(log2_data[,key_field_quantification])
      if (is_normal$p.value >= normality_p_value_threshold) {
        fields_to_return <- c("transformed_data", "normal_distribution", "zscore", "tstudent", "tdist", "to_remove") # init the fields that will be used along this process
        fields_to_process <- c("transformed_data", "zscore", "tstudent", "tdist") # init the fields that will be used along this process
        retrieved_data_df$transformed_data <- transformed_data[,key_field_quantification] # transformed_data <- log2_data[,key_field_quantification]
        retrieved_data_df$normal_distribution <- TRUE
        message(paste("[AP]\t\t!!Luky!!! Data follow Normal distribution!"))
        # avg <- mean(log2_data[,key_field_quantification])
        # std <- sd(log2_data[,key_field_quantification])
        # retrieved_data_df$zscore <- (retrieved_data_df$log2data - avg)/std
        retrieved_data_df$zscore <- scale(retrieved_data_df$transformed_data)
        count <- nrow(retrieved_data_df)
        retrieved_data_df$tstudent <- sqrt(((count*(count-2)*(retrieved_data_df$zscore)^2)/((count-1)^2-count*(retrieved_data_df$zscore)^2)))
        retrieved_data_df$tdist <- dt(retrieved_data_df$tstudent, df = count-2)
        retrieved_data_df$to_remove <- apply(retrieved_data_df[fields_to_process], 1, function(x) {ifelse((x[4] < outlier_p_value_threshold & x[2] < 0), TRUE, FALSE)}) # [which(retrieved_data_df$tdist < outlier_p_value_threshold & retrieved_data_df$zscore < 0),] 
        # append results
        if (length(df_outlist_samplestofilter) == 0){ # first loop
          df_outlist_samplestofilter <- retrieved_data_df[fields_to_return]
        } else {
          df_outlist_samplestofilter <- rbind(df_outlist_samplestofilter, retrieved_data_df[fields_to_return])
        }
      } else {
        message(paste("[AP]\t\tData do not follow Normal distribution, skip it for the test of outliers."))
      } # if (is_normal$p.value >= normality_p_value_threshold)
      
      # write.csv(retrieved_data_df, file = paste("workbench Sept 2017/", pool_id, ".csv", sep = ''), row.names = FALSE)
      
    } else {
      message(paste("[AP]\t\tERROR: Not rows remaining in data frame."))
    } # if (nrow(retrieved_data_df) > 0)
    
  } # else branch # if (per_pool_test)
  
  return (df_outlist_samplestofilter)
  
}



#####################
## stringCombination
#####################
##
##
## try the function:
# strings <- list("Pippo", "Pluto", "Topolino", "Paperino", "Merlino")
# stringCombination_2k(strings)

stringCombination_2k <- function(strings){
  if(length(strings) == 0){
    message("The provided list is empty")
    return()
  }
  if(class(strings) != "list"){
    message("The provided object is not a list")
    return()
  }
  
  N <- length(strings)
  N_comb <- factorial(N)/(factorial(N-2)*factorial(2))
  combList <- list()
  
  i <- 1
  for(j in 1:(N-1)){
    for(k in (j+1):N){
      combList[[i]] <- c(strings[[j]], strings[[k]])
      i <- i + 1
    }
  }
  
  combList
}

###############################################################
#' @title Mask matrix with a template matrix
#' 
#' @author Andrea Calabria
#' @details version 1.1 (used in dat afiltering), 2018-12-19
#'
#' @rdname maskMatrixFromOtherMatrix
#' @docType methods
#' @aliases maskMatrix
#'
#' @param df an input dataframe, even with annotation (5 columns). All columns must be numberic columns (except for annotation_cols).
#' @param df_mask a df that acts as mask. Columns and Rows must be of the same name (not necessary order) of the columns in df. It is assumed that rownames are defined as: <chr number>_<locus>_<strand>_<gene name>_<gene strand> (the canonical id_cols of annotation_cols) or the order and number of the annotation cols annotation_cols; if id_cols is c("chr", "integration_locus", "strand") then the rowname would be chr1_12344_-.
#' @param annotation_cols a vector of column names for annotation (if available). Default: c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"). If rows
#' @param replace_NA_with0 boolean, default TRUE.
#' @param colnames_to_exlude a vector containing the column names to exclude i.e. c("all"), the default.
#'
#' @return a filterd (masked) df.
#' @usage TODO
#' @note : 
#' The aim of the function is to acquire an input df (to be filtered) and a mask df (those 2 df MUST have the same number of rows and columns (checks to do)). Then, using a reshaping approach, the df_mask will be used to filter out data (anche cells) in the source df.
#' The following code is for testing the function, assuming that in the folder "data" you have the appropriate files. TODO
#' 
#'
###############################################################
maskMatrixFromOtherMatrix <- function(df, 
                                      df_mask, 
                                      annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                                      replace_NA_with0 = T,
                                      colnames_to_exlude = c("all"),
                                      apply_rownames_fromannotationcols = F,
                                      chunk_size = NULL,
                                      max_cols_limit_forcasting = 1600) {
  library(reshape2)
  library(dplyr)
  df_join <- NULL
  valid_cols <- NULL # init var of valid cols
  message(paste("[AP]\tMasking df with another df"))
  # check col size according to max_cols_limit_forcasting
  split_cast_by_col_number <- FALSE
  if (ncol(df)>max_cols_limit_forcasting | ncol(df_mask)>max_cols_limit_forcasting) {
    split_cast_by_col_number <- TRUE
  }
  if (length(annotation_cols) > 0) { # in case of having annotation cols
    if (apply_rownames_fromannotationcols) {
      message(paste("[AP]\t\tCheck df consistency, applying row ID based on annotations"))
      rownames(df) <- apply(df[annotation_cols], 1, paste, collapse="_") # {paste(x[1], x[2], x[3], x[4], sep = "_")}
      rownames(df_mask) <- apply(df_mask[annotation_cols], 1, paste, collapse="_") # {paste(x[1], x[2], x[3], x[4], sep = "_")}
    }
    if (length(intersect(rownames(df), rownames(df_mask))) == 0) {
      message(paste("[AP]\t\tERROR: Inconsistency on the ROW NAMES, depending on annotation_cols"))
    } else {
      valid_cols <- setdiff(intersect(colnames(df_mask), colnames(df)), colnames_to_exlude) # the columns in common between dfs
      message(paste("[AP]\t\tThe following columns are not present in df (but only in df_mask):\n", paste( setdiff(colnames(df_mask), c(colnames(df), colnames_to_exlude)), sep = '\n') ))
      
      message(paste("[AP]\t\tReshaping data"))
      df_melt <- NULL
      df_mask_melt <- NULL
      
      if (length(chunk_size) > 0) {
        # grant a min chunk size as the delta between row numbers
        delta_rows_df_and_dfmask <- ifelse(nrow(df) >= nrow(df_mask), nrow(df) - nrow(df_mask), 0)
        min_chunk_size <- delta_rows_df_and_dfmask + 1
        if (chunk_size < min_chunk_size) {
          message(paste("[AP]\t\tAdjusting chunk_size to", min_chunk_size, "since the number of rows of the input df is respetively (df, df_mask):", nrow(df), nrow(df_mask)))
          chunk_size <- min_chunk_size
        }
        steps <- seq(from = 1, to = nrow(df), by = chunk_size)
        for (i in seq(1, length(steps))) {
          message(paste("[AP]\t\t...processing chunk", i, "of", length(steps), ":", round(i*100/length(steps), digits = 2), "%"))
          if (i < length(steps)) { # all rows but NOT the last step
            df_melt_tmp <- melt(data = df[steps[i]:(steps[i+1]-1), valid_cols], id.vars = annotation_cols) # melt df
            df_melt_tmp <- df_melt_tmp[which(df_melt_tmp$value > 0 | !is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
            # df_melt_tmp <- df_melt_tmp[which(!is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
            df_melt <- rbind(df_melt, df_melt_tmp) # update data
            df_mask_melt_tmp <- melt(data = df_mask[steps[i]:(steps[i+1]-1), valid_cols], id.vars = annotation_cols) # melt mask df
            df_mask_melt_tmp <- df_mask_melt_tmp[which(df_mask_melt_tmp$value > 0 | !is.na(df_mask_melt_tmp$value)),] # remove missing data
            # df_mask_melt_tmp <- df_mask_melt_tmp[which(!is.na(df_mask_melt_tmp$value)),] # remove missing data
            df_mask_melt <- rbind(df_mask_melt, df_mask_melt_tmp) # update df
          } else { # last step
            df_melt_tmp <- melt(data = df[steps[i]:nrow(df), valid_cols], id.vars = annotation_cols)
            df_melt_tmp <- df_melt_tmp[which(df_melt_tmp$value > 0 | !is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
            # df_melt_tmp <- df_melt_tmp[which(!is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
            df_melt <- rbind(df_melt, df_melt_tmp)
            df_mask_melt_tmp <- melt(data = df_mask[steps[i]:nrow(df_mask), valid_cols], id.vars = annotation_cols)
            df_mask_melt_tmp <- df_mask_melt_tmp[which(df_mask_melt_tmp$value > 0 | !is.na(df_mask_melt_tmp$value)),] # remove missing data
            df_mask_melt_tmp <- df_mask_melt_tmp[which(!is.na(df_mask_melt_tmp$value)),] # remove missing data
            df_mask_melt <- rbind(df_mask_melt, df_mask_melt_tmp)
          }
        }
        # df_melt <- melt(data = df[valid_cols], id.vars = annotation_cols)
        # df_melt <- df_melt[which(df_melt$value > 0),]
        # df_mask_melt <- melt(data = df_mask[valid_cols], id.vars = annotation_cols)
        # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
      } else {
        if (nrow(df) > 1000000) { 
          message (paste("[AP\t\tWARNING: this step could fail due to the hign number of rows!"))
        } # if (nrow(df) > 1000000) { 
        df_melt <- melt(data = df[valid_cols], id.vars = annotation_cols)
        # df_melt <- df_melt[which(df_melt$value > 0),]
        df_mask_melt <- melt(data = df_mask[valid_cols], id.vars = annotation_cols)
        # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
      } # if (length(chunk_size) > 0)
      
      # df_melt <- df_melt[which(df_melt$value > 0),]
      # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
      message(paste("[AP]\t\tJoining data"))
      df_join <- merge(x = df_melt, y = df_mask_melt, by = c(annotation_cols, "variable"), all.x = T)
      names(df_join) <- c(annotation_cols, "Sample", "QuantificationDF", "QuantificationMASK")
      # df_join$QuantificationMASK <- round(df_join$QuantificationMASK)
      
      message(paste("[AP]\t\tFiltering results"))
      df_join <- df_join[which(!is.na(df_join$QuantificationMASK)),]
    } # if (length(intersect(rownames(df), rownames(df_mask))) == 0) 
   
  } else { # without annotation_cols
    message(paste("[AP]\t\tYou did not provide annotation cols, thus ROW names will be used as they are in data reshaping (consistency)"))
    valid_cols <- setdiff(intersect(colnames(df_mask), colnames(df)), colnames_to_exlude)
    message(paste("[AP]\t\tThe following columns are not present in df (but only in df_mask):\n", paste( setdiff(colnames(df_mask), c(colnames(df), colnames_to_exlude)), sep = '\n') ))
    
    message(paste("[AP]\t\tReshaping data"))
    # df_melt <- melt(data = as.matrix(df[valid_cols]))
    # df_melt <- df_melt[which(df_melt$value > 0),]
    # df_mask_melt <- melt(data = as.matrix(df_mask[valid_cols]))
    # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
    
    df_melt <- NULL
    df_mask_melt <- NULL
    if (length(chunk_size) > 0) {
      delta_rows_df_and_dfmask <- ifelse(nrow(df) >= nrow(df_mask), nrow(df) - nrow(df_mask), 0)
      min_chunk_size <- delta_rows_df_and_dfmask + 1
      if (chunk_size < min_chunk_size) {
        message(paste("[AP]\t\tAdjusting chunk_size to", min_chunk_size, "since the number of rows of the input df is respetively (df, df_mask):", nrow(df), nrow(df_mask)))
        chunk_size <- min_chunk_size
      }
      steps <- seq(from = 1, to = nrow(df), by = chunk_size)
      for (i in length(steps)) {
        message(paste("[AP]\t\t...processing chunk", i, "of", length(steps), ":", round(i*100/length(steps), digits = 2), "%"))
        if (i < length(steps)) { # all rows but NOT the last step
          # df_melt_tmp <- melt(data = df[steps[i]:(steps[i+1]-1), valid_cols], id.vars = annotation_cols) # melt df
          df_melt_tmp <- melt(data = as.matrix(df[steps[i]:(steps[i+1]-1), valid_cols]))
          df_melt_tmp <- df_melt_tmp[which(df_melt_tmp$value > 0 | !is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
          # df_melt_tmp <- df_melt_tmp[which(!is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
          df_melt <- rbind(df_melt, df_melt_tmp) # update data
          # df_mask_melt_tmp <- melt(data = df_mask[steps[i]:(steps[i+1]-1), valid_cols], id.vars = annotation_cols) # melt mask df
          df_mask_melt_tmp <- melt(data = as.matrix(df_mask[steps[i]:(steps[i+1]-1), valid_cols]))
          df_mask_melt_tmp <- df_mask_melt_tmp[which(df_mask_melt_tmp$value > 0 | !is.na(df_mask_melt_tmp$value)),] # remove missing data
          # df_mask_melt_tmp <- df_mask_melt_tmp[which(!is.na(df_mask_melt_tmp$value)),] # remove missing data
          df_mask_melt <- rbind(df_mask_melt, df_mask_melt_tmp) # update df
        } else { # last step
          df_melt_tmp <- melt(data = df[steps[i]:nrow(df), valid_cols], id.vars = annotation_cols)
          df_melt_tmp <- df_melt_tmp[which(df_melt_tmp$value > 0 | !is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
          # df_melt_tmp <- df_melt_tmp[which(!is.na(df_melt_tmp$value)),] # remove missing data (reduce data)
          df_melt <- rbind(df_melt, df_melt_tmp)
          df_mask_melt_tmp <- melt(data = df_mask[steps[i]:nrow(df_mask), valid_cols], id.vars = annotation_cols)
          df_mask_melt_tmp <- df_mask_melt_tmp[which(df_mask_melt_tmp$value > 0 | !is.na(df_mask_melt_tmp$value)),] # remove missing data
          # df_mask_melt_tmp <- df_mask_melt_tmp[which(!is.na(df_mask_melt_tmp$value)),] # remove missing data
          df_mask_melt <- rbind(df_mask_melt, df_mask_melt_tmp)
        }
      }
      # df_melt <- melt(data = df[valid_cols], id.vars = annotation_cols)
      # df_melt <- df_melt[which(df_melt$value > 0),]
      # df_mask_melt <- melt(data = df_mask[valid_cols], id.vars = annotation_cols)
      # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
    } else {
      if (nrow(df) > 1000000) { 
        message (paste("[AP\t\tWARNING: this step could fail due to the hign number of rows!"))
      } # if (nrow(df) > 1000000) { 
      df_melt <- melt(data = df[valid_cols], id.vars = annotation_cols)
      # df_melt <- df_melt[which(df_melt$value > 0),]
      df_mask_melt <- melt(data = df_mask[valid_cols], id.vars = annotation_cols)
      # df_mask_melt <- df_mask_melt[which(df_mask_melt$value > 0),]
    } # if (length(chunk_size) > 0)
    
    message(paste("[AP]\t\tJoining data"))
    df_join <- merge(x = df_melt, y = df_mask_melt, by = c(annotation_cols, "variable"), all.x = T)
    names(df_join) <- c(annotation_cols, "Sample", "QuantificationDF", "QuantificationMASK")
    # df_join$QuantificationMASK <- round(df_join$QuantificationMASK)
    
    message(paste("[AP]\t\tFiltering results"))
    df_join <- df_join[which(!is.na(df_join$QuantificationMASK)),]
  } # if (length(annotation_cols) > 0) 
  
  ### TODO ### actually, is working statically!!
  # rebuild the OUT matrix
  message(paste("[AP]\t\tCasting df (WARNING: SO FAR, this cast is STATIC (!!!) on the annotation cols 'chr + integration_locus + strand + GeneName + GeneStrand'. TODO generalization "))
  df_updated <- NULL
  if (split_cast_by_col_number) {
    # filter by column ID
    n_col_sets <- 2 # number of sets for splitting. now 2 since the join is easier with 2...!
    message(paste("[AP]\t\t\tSplitting columns in", n_col_sets, "df (by filtering data)"))
    col_list <- levels(factor(df_join$Sample)) # static list of columns
    this_col_list <- NULL
    if (length(df_updated) == 0) { # first loop
      this_col_list <- col_list[1:round(length(col_list)/n_col_sets)]
      df_join_slice <- df_join[which(df_join$Sample %in% this_col_list),]
      df_updated <- dcast(data = df_join_slice, chr + integration_locus + strand + GeneName + GeneStrand ~ Sample, value.var = "QuantificationDF", fun.aggregate = mean)
    } else {
      this_col_list <- setdiff(col_list, this_col_list)
      df_join_slice <- df_join[which(df_join$Sample %in% this_col_list),]
      df_updated_slice <- dcast(data = df_join_slice, chr + integration_locus + strand + GeneName + GeneStrand ~ Sample, value.var = "QuantificationDF", fun.aggregate = mean)
      # fai join
      df_updated <- merge(x = df_updated, y = df_updated_slice, by = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"), all = T)
      # df_updated <- full_join(x = df_updated, y = df_updated_slice, by = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"), copy = FALSE)
    } # if (length(df_updated) == 0)
  } else {
    df_updated <- dcast(data = df_join, chr + integration_locus + strand + GeneName + GeneStrand ~ Sample, value.var = "QuantificationDF", fun.aggregate = mean)
  } # if (split_cast_by_col_number) 
  
  if (replace_NA_with0) {
    message(paste("[AP]\t\tReplacing NA with 0"))
    df_updated[is.na(df_updated)] <- 0
  }
  
  return (df_updated)
}

###############################################################
#' @title Remove 0 from tails by IS in df
#' 
#' @author Andrea Calabria
#' @details version 0.1 (used in data filtering), 2019-02-28
#'
#' @rdname removeTail0byRow
#' @docType methods
#' @aliases remove0fromISRowTails
#'
#' @param df an input dataframe, even with annotation (5 columns). All columns must be numberic columns (except for annotation_cols).
#' @param annotation_cols a vector of column names for annotation (if available). Default: c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"). If rows
#' @param colnames_to_exlude a vector containing the column names to exclude i.e. c("all"), the default.
#'
#' @return a df with same input size but with 0 tails in rows removed
#' @usage TODO
#' @note : 
#' The aim of the function is to remove 0s from tails in each rows: [0, 0, 1, 0, 3, 0] -> [NA, NA, 1, 0, 3, NA].
#' The order of the columns is crucial! Input columns must be already in the desired order.
#'
###############################################################
removeTail0byRow <- function(df, 
                             annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                             colnames_to_exlude = c("all"),
                             test = FALSE) {
  data_cols <- setdiff(colnames(df), c(annotation_cols, colnames_to_exlude)) # the columns in common between dfs
  
  # if testing == TRUE
  df_sample <- NULL
  if (test) {
    # df_sample <- read.csv(file = paste("data/test.matrix.foo.for0TailRemoval.csv", sep = ""), header=TRUE, fill=T, sep='\t', check.names = F, row.names = 1)
    df_sample <- read.csv(file = paste("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/TIGET/Workbench/isatk/script/R/data/test.matrix.foo.for0TailRemoval.csv", sep = ""), header=TRUE, fill=T, sep='\t', check.names = F, row.names = 1)
    df <- df_sample
  }
  
  # df_notails <- data.frame(t(sapply(rownames(df[data_cols]), function(x) {
  #   ifelse(dfallp_projectcentric_tracking_date[x,]==dfallp_tracking_date_minDate[x,], 1, 0)
  # } )))
  df[is.na(df)] <- 0 # first grant that you have matrix with 0 (no NA)
  df_notails <- cbind(
      df[annotation_cols],
      data.frame(t(apply(df[data_cols], 1, function(x) {
            .val_vector <- x
            .index_firstVal <- head(which(x>0), 1)
            .index_lastVal <- tail(which(x>0), 1)
            for (k in seq(1, length(x))) {
              if (x[k] == 0 & (k < .index_firstVal | k > .index_lastVal)) {
                x[k] <- NA
              }
            }
            x # the output adjusted
          }
        ) # apply
      ) # t
    ) # data.frame
  ) # cbind
  
  # results
  return (df_notails)
} # removeTail0byRow




###############################################################
#' @title Reshape columns over time as theri flags given a matrix (df) and their metadata
#' 
#' @author Andrea Calabria
#' @details version 0.2 (30 July 2019) 
#'
#' @rdname flagColumns_byTime
#' @docType methods
#' @aliases flagColumns_byTime
#'
#' @param df an input dataframe of IS matrix.
#' @param metadata_df a dataframe of metadata, row names must be an overset of the columns of the input df.
#' @param time_field a string/character for the field of time name (a colname of the metadata_df)
#' @param starting_data_col_index an integer value >0 with the index of the data columns (excluding annotations, etc.)
#'
#' @return a dataframe with the same size of the input df but with different colnames, taken from metadata_df
#' @usage todo
#'
###############################################################
###############################################################
# flagColumnsByMetadata <- function(df, metadata_df, key_field, starting_data_col_index = 6, number_of_last_cols_to_remove = 0) {
flagColumns_byTime <- function(df, 
                               starting_data_col_index = 6, 
                               number_of_last_cols_to_remove = 1, 
                               timepoint_position_incolnames = 3, 
                               colname_sep = '_',
                               addMissingCols = T,
                               convertToDecimal = F) {
  # require(DescTools)
  message(paste("[AP]\tReshaping a df using Flags by time."))
  df_clean <- df[colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])]
  # be sure that this matrix is really full and sorted...
  if (addMissingCols) {
    df_clean <- addCombinatorialMissingColumns_byTime(df = df_clean, 
                                                   starting_data_col_index = 1, 
                                                   number_of_last_cols_to_remove = 0, 
                                                   timepoint_position_incolnames = timepoint_position_incolnames, 
                                                   colname_sep = colname_sep)
  }
  # binarize matrix
  df_clean[is.na(df_clean)] <-  0
  df_clean[df_clean > 0] <- 1
  # col details
  col_df <- as.data.frame(t(as.data.frame(lapply(strsplit(colnames(df_clean), "_", fixed = TRUE), function(x) {c(x)}))))
  # use the position of time to do loop for
  present_col_suffix <- levels(factor(col_df[, timepoint_position_incolnames]))
  present_col_prefix <- levels(factor( apply(col_df[, -(timepoint_position_incolnames)], 1, paste, collapse = colname_sep) ) )
  message(paste("[AP]\t-> This will be the ORDER of the binary (base 2) flag:\n\t\t", paste0(present_col_prefix, collapse = '\n\t\t') ))
  
  message(paste("[AP]\t-> Processing timepoint"))
  # outdf as empty df
  out_df <- NULL # the final numberof columns will be the same of length present_col_suffix
  unique_base2_strings <- NULL
  for (timepoint in present_col_suffix) { 
    message(paste("[AP]\t\t", timepoint))
    # slice matrix [same idea of apply(expand.grid(present_col_prefix, present_col_suffix), 1, paste, collapse = colname_sep) ]
    df_clean_timeslice <- df_clean[paste(present_col_prefix, timepoint, sep = colname_sep)]
    df_clean_timeslice_flag <- NULL
    if (convertToDecimal) {
      df_clean_timeslice_flag <- as.data.frame(apply(df_clean_timeslice, 1, function(x) {
                                                base2_string <- paste0(x, collapse = '')
                                                myBinToDec(base2_string)
                                                })
                                              )  
    } else {
      df_clean_timeslice_flag <- as.data.frame(apply(df_clean_timeslice, 1, paste, collapse = ''))
    } # if convertToDecimal
    
    unique_base2_strings <- c(unique_base2_strings, apply(df_clean_timeslice, 1, paste, collapse = '') )
    names(df_clean_timeslice_flag) <- timepoint
    if (length(out_df) == 0) {
      out_df <- df_clean_timeslice_flag
    } else {
      out_df <- cbind(out_df, df_clean_timeslice_flag)  
    } # if (length(out_df) == 0)
    
  } # for (timepoint in present_col_suffix) 
  # do the unique function on the vector
  unique_base2_strings <- unique(unique_base2_strings)
  
  # return a df of all the occurrences
  message(paste("[AP]\t-> Building the bag of cases (df)"))
  bagcases <- as.data.frame(t(as.data.frame(lapply(strsplit(unique_base2_strings, '', fixed = T), function(x) {c(x)}))))
  names(bagcases) <- present_col_prefix
  bagcases$NumberBase2 <- unique_base2_strings
  bagcases$NumberBase10 <- as.vector(sapply(unique_base2_strings, myBinToDec))
  rownames(bagcases) <- as.vector(sapply(unique_base2_strings, myBinToDec)) # gdf_hltfu_filtered_flags$base2_cases  BinToDec(gdf_hltfu_filtered_flags$base2_cases)
  
  return (list("flag_df" = out_df, 
               "timpoints" = present_col_suffix, 
               "bagcases_df" = bagcases,
               "sample_order" = present_col_prefix,
               "base2_cases" = unique_base2_strings
               ))
}


###############################################################
#' @title Get occurrences of elements by row in a given df.
#' 
#' @author Andrea Calabria
#' @details version 0.1 (30 July 2019) 
#'
#' @rdname searchOccurrencesInDf_byRow
#' @docType methods
#' @aliases searchOccurrencesInDf_byRow
#'
#' @param df an input dataframe of IS matrix.
#' @param return_as_df boolean, set T if you want return a dataframe instead of a matrix
#'
#' @return a dataframe or a matrix of occurrences
#' @description Similar to table(df) but works faster and by rows.
#' 
###############################################################
searchOccurrencesInDf_byRow <- function(df, return_as_df = T) {
  message(paste("[AP]\tCount occurrences of values contained in the input df, by row."))
  levels <- sort(unique(do.call(c, df))) #all unique values in df
  out <- sapply(levels, function(x) rowSums(df == x) ) #count occurrences of x in each row
  colnames(out) <- levels
  if (return_as_df) {
    return ( as.data.frame(out))
  } else {
    return (out)
  }
}

###############################################################
#' @title Binary stirng to Decimal value
#' 
#' @author Andrea Calabria
#' @details version 0.1 (30 July 2019) 
#'
#' @rdname myBinToDec
#' @docType methods
#' @aliases myBinToDec
#'
#' @param x a vector
#'
#' @return an integer
#' @description Dimilar to BinToDec of DescTools. Required to bypass this package, not available for older verions of R.
#' 
###############################################################
myBinToDec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}


###############################################################
#' @title Reshape columns over time as theri flags given a matrix (df) and their metadata
#' 
#' @author Andrea Calabria
#' @details version 0.1 (30 July 2019) 
#'
#' @rdname addCombinatorialMissingColumns_byTime
#' @docType methods
#' @aliases addCombinatorialMissingColumns_byTime
#'
#' @param df an input dataframe of IS matrix.
#' @param timepoint_position_incolnames a string/character for the field of time name (a colname of the metadata_df). defauult = 3.
#' @param starting_data_col_index an integer value >0 with the index of the data columns (excluding annotations, etc.). default = 6.
#' @param number_of_last_cols_to_remove integer: number of last columns to remove, default = 0.
#' @param colname_sep separator used in col names. default = '-'
#'
#' @return a dataframe with the same number of rows by with all combinatorial samples (by time) sort by sample and then by time. 
#' @usage todo
#' @description Given an input df (boundaries explicitly defined), this function will return a df expanded in columns by adding the missing timepoints (sorting columns by a fixed condition time based).
#'  Prerequisites: the df MUST have time points. Time points MUST have padding (or the sorting will NOT be correct). Separator is '_'.
#' 
###############################################################
addCombinatorialMissingColumns_byTime <- function(df, 
                                                  starting_data_col_index = 6, 
                                                  number_of_last_cols_to_remove = 0, 
                                                  timepoint_position_incolnames = 3, 
                                                  colname_sep = '_') {
  message(paste("[AP]\tAdding missing columns"))
  df_clean <- df[colnames(df[starting_data_col_index:(length(colnames(df))-number_of_last_cols_to_remove)])]
  # get df of details
  col_df <- as.data.frame(t(as.data.frame(lapply(strsplit(colnames(df_clean), colname_sep, fixed = TRUE), function(x) {c(x)}))))
  
  # fill missing columns (as time point) with NA or 0 values for the actual rows
  # output_suffix_tohave:: the last fidl of the label names (that is the time point string as text)
  present_col_suffix <- levels(factor(col_df[, timepoint_position_incolnames]))
  present_col_prefix <- levels(factor( apply(col_df[, -(timepoint_position_incolnames)], 1, paste, collapse = colname_sep) ) )
  
  # expand.grid(present_col_prefix, present_col_suffix)
  all_combinations <- apply(expand.grid(present_col_prefix, present_col_suffix), 1, paste, collapse = colname_sep) # this is THE OUTPUT ORDER
  
  # idea: add at the end all missing columns and then sort by name (as the previous order)
  missing_col_suffix_to_include <- setdiff(all_combinations, colnames(df_clean))
  df_missing_cols <- setNames(object = data.frame(matrix(ncol = length(missing_col_suffix_to_include), nrow = nrow(df_clean))), nm = missing_col_suffix_to_include)
  
  # bind all data, then sort
  out_df <- cbind(df_clean, df_missing_cols)
  out_df <- out_df[all_combinations]
  return (out_df)
}



### from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  require(ggplot2)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plot_cluster <- function(data, var_cluster, palette = "Spectral", xname = "V1", yname = "V2", plot_title)  
{
  ggplot(data, aes_string(x=xname, y=yname, color=var_cluster)) +
    geom_point(size=0.5) +
    guides(colour=guide_legend(override.aes=list(size=6))) +
    # xlab("") + ylab("") +
    ggtitle(plot_title) +
    theme_light(base_size=20) +
    theme(
      # axis.text.x=element_blank(),
      # axis.text.y=element_blank(),
      legend.direction = "horizontal", 
      legend.position = "bottom",
      legend.box = "horizontal") + 
    scale_colour_manual(values = colorRampPalette(brewer.pal(11, palette))(length(levels(factor(data[, var_cluster])))) )
  # scale_colour_brewer(palette = palette) 
}




###############################################################
#' @title Run tSNE on IS
#' 
#' @author Andrea Calabria
#' @details version 0.1 (03 August 2019) 
#'
#' @rdname tSNEonIS
#' @docType methods
#' @aliases tSNEonIS
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#' @param prefix_data_columns prefix to assign to output data columns
#' @param prefix_clustering_columns prefix to assign to output data clustering columns
#' @param factors_to_search number of factors to look for, default = 5.
#' @param do_plots boolean, do plots? default = T.
#' @param plot_filename full path or output plot file name.
#' @param plot_title plot title to use.
#' @param seed_number seed
#'
#' @return a dataframe with the same number of rows and new columns representing the clustering tSNE results.
#' @usage todo
#' @description ToDo
#' 
###############################################################
tSNEonIS <- function(df, 
                     data_columns,
                     prefix_data_columns = "Profiling", 
                     prefix_clustering_columns = "profiling", 
                     factors_to_search = 5,
                     do_plots = T,
                     plot_filename = "plot_file",
                     plot_title = "tSNE", 
                     seed_number = 123) {
  # do the same but removing not shared elements
  set.seed(seed_number)
  message(paste("[AP]\tCalculating tSNE on input data (only on selected columns:", paste(data_columns, collapse = " "), ")."))
  tsne_model_1 <- Rtsne(as.matrix(df[data_columns]), 
                        check_duplicates = FALSE, 
                        pca = TRUE, 
                        perplexity = 30, 
                        theta = 0.2, 
                        dims = 2, 
                        normalize = F,
                        num_threads = 0)
  
  ## getting the two dimension matrix
  d_tsne_1 <- as.data.frame(tsne_model_1$Y)
  rownames(d_tsne_1) <- rownames(df)
  datacols <- paste0(prefix_data_columns, colnames(d_tsne_1))
  names(d_tsne_1) <- datacols
  
  ## do clustering now
  message(paste("[AP]\tCalculating clusters on tSNE results:"))
  ## Creating k-means clustering model, and assigning the result to the data used to create the tsne
  message(paste("[AP]\t\t-> k-means with k =", factors_to_search))
  fit_cluster_kmeans <- kmeans(scale(d_tsne_1), factors_to_search)  # original
  # fit_cluster_kmeans=kmeans(scale(d_tsne_1_tracking_timeShared_no0), 5)  
  d_tsne_1$cluKmeans <- factor(fit_cluster_kmeans$cluster)
  names(d_tsne_1)[names(d_tsne_1) == 'cluKmeans'] <- paste(prefix_clustering_columns, "cluKmeans", sep = '_')
  
  ## Creating hierarchical cluster model, and assigning the result to the data used to create the tsne
  message(paste("[AP]\t\t-> hierarchical with factors k =", factors_to_search))
  fit_cluster_hierarchical <- hclust(dist(scale(d_tsne_1[datacols])))
  d_tsne_1$cluHierarchical = factor(cutree(fit_cluster_hierarchical, k=factors_to_search))
  names(d_tsne_1)[names(d_tsne_1) == 'cluHierarchical'] <- paste(prefix_clustering_columns, "cluHierarchical", sep = '_')
  
  # merge datasets to have them all in the same df
  message(paste("[AP]\tMerging results."))
  out_matrix_complete <- merge(x = df, y = d_tsne_1, by = 0, all.x = T)
  rownames(out_matrix_complete) <- out_matrix_complete$Row.names
  out_matrix_complete <- out_matrix_complete[, !(names(out_matrix_complete) %in% c("Row.names"))]
  
  # Plotting the cluster models onto t-SNE output
  # Now time to plot the result of each cluster model, based on the t-SNE map.
  
  if (do_plots) {
    message(paste("[AP]\tPlotting results"))
    # first produce plots
    plot_noclusters <- ggplot(d_tsne_1, aes_string(x=names(d_tsne_1)[1], y=names(d_tsne_1)[2])) +  
      geom_point(size=0.25) +
      guides(colour=guide_legend(override.aes=list(size=6))) +
      # xlab("") + ylab("") +
      ggtitle(plot_title) +
      theme_light(base_size=20) +
      theme(axis.text.x=element_blank(),
            axis.text.y=element_blank()) +
      scale_colour_brewer(palette = "Set2")
    
    plot_k <- plot_cluster(d_tsne_1, var_cluster = paste(prefix_clustering_columns, "cluKmeans", sep = '_'), 
                           palette = "Spectral", xname = datacols[1], yname = datacols[2], 
                           plot_title = paste0(plot_title, " k-means (k ", factors_to_search, ")") )  
    plot_h <- plot_cluster(d_tsne_1, var_cluster = paste(prefix_clustering_columns, "cluHierarchical", sep = '_'), 
                           palette = "Spectral", xname = datacols[1], yname = datacols[2], 
                           plot_title = paste0(plot_title, " hierarchical (k ", factors_to_search, ")" ) )
    
    plot_k_contour <- ggplot(d_tsne_1, aes_string(x=datacols[1], y=datacols[2], color=paste(prefix_clustering_columns, "cluKmeans", sep = '_') )) +
      geom_density_2d(alpha = 0.7) +
      geom_point(size = 1, alpha = 0.9) +
      guides(colour=guide_legend(override.aes=list(size=6))) +
      # geom_contour() +
      # xlab("") + ylab("") +
      ggtitle(paste0(plot_title, " k-means (k ", factors_to_search, ")")) +
      theme_light(base_size=20) +
      theme(
        # axis.text.x=element_blank(),
        # axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") + 
      scale_colour_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(factors_to_search) )
    
    plot_h_contour <- ggplot(d_tsne_1, aes_string(x=datacols[1], y=datacols[2], color=paste(prefix_clustering_columns, "cluHierarchical", sep = '_') )) +
      geom_density_2d(alpha = 0.7) +
      geom_point(size = 1, alpha = 0.9) +
      guides(colour=guide_legend(override.aes=list(size=6))) +
      # geom_contour() +
      # xlab("") + ylab("") +
      ggtitle(paste0(plot_title, " hierarchical (k ", factors_to_search, ")" )) +
      theme_light(base_size=20) +
      theme(
        # axis.text.x=element_blank(),
        # axis.text.y=element_blank(),
        legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.box = "horizontal") + 
      scale_colour_manual(values = colorRampPalette(brewer.pal(11, "Spectral"))(factors_to_search) )
    
    # then do plots in pdf
    pdf(plot_filename, height=9, width=12)
    print (plot_noclusters)
    print (plot_h)
    print (plot_k)
    print (plot_h_contour)
    print (plot_k_contour)
    
    print (
      grid.arrange(plot_k, plot_h,  
                   plot_k_contour, plot_h_contour, ncol=2)    
    )
    dev.off()
    
  } # if do plots
  
  return (out_matrix_complete)
}


###############################################################
#' @title Plot heatmaps
#' 
#' @author Andrea Calabria
#' @details version 0.1 (03 August 2019) 
#'
#' @rdname plotHeatmapsOvertime
#' @docType methods
#' @aliases plotHeatmapsOvertime
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#' @param plot_filename full path or output plot file name.
#' @param plot_title plot title to use.
#' @param seed_number seed
#' @param color_names_sorted_by_number Colors to use, a vector of named colors
#'
#' @return Plots
#' @usage todo
#' @description ToDo
#' 
###############################################################
plotHeatmapsOvertime <- function(df, 
                                 data_columns,
                                 plot_filename = "plot_file.heatmap.pdf",
                                 plot_title = "Heatmap", 
                                 color_names_sorted_by_number,
                                 seed_number = 123, height=9, width=12) {
  pdf(file = plot_filename, height=height, width=width)
  plot_main_title <- plot_title
  heatmap.2(as.matrix(df[data_columns]),
            key=TRUE,
            scale="none",
            Rowv=T, Colv=F,
            density.info="none", trace="none",
            labRow = FALSE,
            # cexRow=1.0, cexCol=1.0,
            col = col2hex(as.character(color_names_sorted_by_number)),
            main = plot_main_title,
            dendrogram="row",
            # colsep=seq(1,dim(patients_summary_onTime_bySeparatedCells_Tissue)[2]-1),
            # rowsep=seq(1,dim(patients_summary_onTime_bySeparatedCells_Tissue)[1]-1),
            # sepcolor="gray", sepwidth=c(0.005,0.005), margin=c(6, 15),
            # cellnote=patients_summary_onTime_bySeparatedCells_Tissue,
            # notecex=1.0,
            # notecol="black",
            na.color=par("bg")
  )
  heatmap.2(as.matrix(df[data_columns]),
            key=TRUE,
            scale="none",
            Rowv=T, Colv=T,
            density.info="none", trace="none",
            labRow = FALSE,
            # cexRow=1.0, cexCol=1.0,
            col = col2hex(as.character(color_names_sorted_by_number)),
            main = plot_main_title,
            dendrogram="both",
            # colsep=seq(1,dim(patients_summary_onTime_bySeparatedCells_Tissue)[2]-1),
            # rowsep=seq(1,dim(patients_summary_onTime_bySeparatedCells_Tissue)[1]-1),
            # sepcolor="gray", sepwidth=c(0.005,0.005), margin=c(6, 15),
            # cellnote=patients_summary_onTime_bySeparatedCells_Tissue,
            # notecex=1.0,
            # notecol="black",
            na.color=par("bg")
  )
  pheatmap(df[data_columns],
           show_colnames = T, show_rownames = F, 
           cluster_rows = T, cluster_cols = F,
           color = col2hex(as.character(color_names_sorted_by_number)),
           # cutree_rows = 5,
           main = plot_main_title)
  pheatmap(df[data_columns],
           show_colnames = T, show_rownames = F, 
           cluster_rows = T, cluster_cols = T,
           color = col2hex(as.character(color_names_sorted_by_number)),
           # cutree_rows = 5,
           main = plot_main_title)
  dev.off()
  
}

###############################################################
#' @title Plot profiles over time
#' 
#' @author Andrea Calabria
#' @details version 0.1 (03 August 2019) 
#'
#' @rdname plotProfilesBarPlotsPanels
#' @docType methods
#' @aliases plotProfilesBarPlotsPanels
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#' @param plot_filename full path or output plot file name.
#' @param plot_title plot title to use.
#' @param plot_footnote_details footnote to include in the plot; useful for comments.
#' @param height plot height of the putput PDF plot.
#' @param width plot width of the oputput PDF plot.
#' @param seed_number seed
#'
#' @return Plots
#' @usage todo
#' @description Requires specific columns!  "FollowUp", "NIS". If missing, an exception is raised.
#' 
###############################################################
plotProfilesBarPlotsPanels <- function(df_freqs, df_perc,
                                       data_columns,
                                       plot_filename = "plot_file.barplotpanel.pdf",
                                       plot_title = "Panel of profiles", 
                                       plot_footnote_details = "IS classified by relabeled flags, filtered IS not shared over time, relabeled in 3 classes and corrected.",
                                       scale_color_manual_colors_sortedbyname = c(),
                                       seed_number = 123, height=9, width=12) {
  do_plots <- TRUE
  required_cols <- c("FollowUp", "NIS")
  if (length(grep(FALSE, (c("FollowUp", "NIS") %in% c(colnames(df_freqs), colnames(df_perc)))) > 0)) { # some required columns are not there!
    message(paste("[AP]\tERROR: some required columns are not present in input DF. Required cols: ", paste(required_cols, collapse = ' - ')))
    do_plots <- FALSE
  }
  
  if (length(scale_color_manual)>0) {
    plot_freqs_bar <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      geom_text(aes(label=NIS), vjust=1.6, color="black", size=3.5, position = position_stack()) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_bar <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      geom_text(aes(label=round(NIS, 1)), vjust=1.6, color="black", size=3.5, position = position_stack()) +
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_freqs_line <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      # geom_smooth(stat = "smooth", position = "identity", span = 10, alpha = 0.4) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_line <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      # geom_smooth(stat = "smooth", position = "identity", span = 10, alpha = 0.4) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_freqs_line_smooth <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      # geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      # geom_smooth(method = "loess", formula = y ~ splines::ns(x, 3), position = "identity", span = 10, alpha = 0.4, level = 0.75) +
      geom_smooth(method = "loess", formula = y ~ splines::ns(x), position = "identity", alpha = 0.4, level = 0.75) +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_line_smooth <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
      geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
      scale_color_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_fill_manual(values = scale_color_manual_colors_sortedbyname) +
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
  } else {
    plot_freqs_bar <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(aes(label=NIS), vjust=1.6, color="black", size=3.5, position = position_stack()) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_bar <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class)) +
      geom_bar(stat = "identity", position = "stack") +
      geom_text(aes(label=round(NIS, 1)), vjust=1.6, color="black", size=3.5, position = position_stack()) +
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_freqs_line <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      # geom_smooth(stat = "smooth", position = "identity", span = 10, alpha = 0.4) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_line <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      # geom_smooth(stat = "smooth", position = "identity", span = 10, alpha = 0.4) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_freqs_line_smooth <-
      ggplot(df_freqs[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      # geom_smooth(stat = "smooth", position = "identity", span = 0.4) +
      # geom_smooth(method = "loess", formula = y ~ splines::ns(x, 3), position = "identity", span = 10, alpha = 0.4, level = 0.75) +
      geom_smooth(method = "loess", formula = y ~ splines::ns(x), position = "identity", alpha = 0.4, level = 0.75) +
      scale_x_continuous(breaks = seq(0, max(df_freqs$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in absolute numbers"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
    
    plot_percents_line_smooth <-
      ggplot(df_perc[data_columns], aes(x = FollowUp, y = NIS, fill = Class, color = Class), na.rm = T, se = TRUE) +
      geom_point() +
      # geom_smooth(stat = "smooth", position = "identity", span = 0.44) +
      geom_smooth(method = "loess", formula = y ~ splines::ns(x), stat = "smooth", position = "identity", alpha = 0.4, level = 0.75) +
      scale_x_continuous(breaks = seq(0, max(df_perc$FollowUp, na.rm = T), 6) ) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") + 
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      labs(title = paste0("IS in percentage by time"), x = "Months after gene therapy", y = "N. IS", colour = "Class")
  }
  
  if (do_plots) {
    pdf(file = plot_filename, height=height, width=width)
    grid.arrange(plot_freqs_bar, plot_percents_bar, 
                 plot_freqs_line, plot_percents_line, 
                 plot_freqs_line_smooth, plot_percents_line_smooth, 
                 ncol = 2, 
                 top = textGrob(plot_title, gp = gpar(fontsize = 22)),
                 bottom = textGrob(plot_footnote_details, gp = gpar(fontface = 3, fontsize = 18), hjust=1, x = 1)
    )  
    dev.off()
    png(file = gsub(".pdf", ".png", plot_filename), height=height, width=width, units = "in", res = 300)
    print(grid.arrange(plot_freqs_bar, plot_percents_bar, 
                       plot_freqs_line, plot_percents_line, 
                       plot_freqs_line_smooth, plot_percents_line_smooth, 
                       ncol = 2, 
                       top = textGrob(plot_title, gp = gpar(fontsize = 22)),
                       bottom = textGrob(plot_footnote_details, gp = gpar(fontface = 3, fontsize = 18), hjust=1, x = 1)
        )
    )
    dev.off()
  } # if (do_plots) 
  
}

###############################################################
#' @title Get sorted occurrences by column
#' 
#' @author Andrea Calabria
#' @details version 0.1 (03 August 2019) 
#'
#' @rdname getSortedCountsOfLevels
#' @docType methods
#' @aliases getSortedCountsOfLevels
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#'
#' @return df
#' @usage todo
#' @description ToDo
#' 
###############################################################
# per ciascuno dei valori (livelli) presenti nel dataframe intero, calcola occorrenze ordinate! colonna
getColumnSortedCountsOfLevels <- function(df, data_columns, levels_tosearch = NULL) {
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  this_df <- df[data_columns]
  this_df[is.na(this_df)] <- 0 # avoid here inside NA
  levels <- NULL
  if (length(levels_tosearch) == 0) {
    levels <- sort(unique(do.call(c, this_df)))  
  } else {
    levels <- levels_tosearch
  }
  message("[AP]\tRetrieved Levels:", paste(levels, collapse = ' - '))
  df_freq <- as.data.frame(apply(this_df, 2, function(x) { table(factor(x, levels = levels)) }))
  return (df_freq)
}

###############################################################
#' @title Get sharing levels of IS by levels (retrieved in the whole df) over columns
#' 
#' @author Andrea Calabria
#' @details version 0.1 (08 August 2019) 
#'
#' @rdname getSharedISnumber_byLevels
#' @docType methods
#' @aliases getSharedISnumber_byLevels
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#'
#' @return df
#' @usage todo
#' @description This is an extensino (and coding improvement) of the former getSharedISnumber in which I also look for sharing of IS belonging to each single level. Here the matrix of IS does not longer contain counts (freqs or n. reads) but LEVELS! A level is a flag or a category.
#
###############################################################
# per ciascuno dei valori (livelli) presenti nel dataframe intero, calcola lo sharing per colonna
getSharedISnumber_byLevels <- function(df, 
                                       data_columns, 
                                       compact = TRUE) {
  df <- df[data_columns]
  # in questa funzione, a partire dal df (RxC), per ogni colonna c  in C fissi le is di c (ovvero quelle >0 not NA) e conti quante sono condivise restituendo un vettore
  message("[AP]\tConverting input df by adding 0 to NA")
  df[is.na(df)] <- 0 # avoid here inside NA
  if (compact) {
    df <- compactDfByColumns(df)
  }
  levels <- sort(unique(do.call(c, df)))
  message("[AP]\tRetrieved Levels:", paste(levels, collapse = ' - '))
  
  message(paste("[AP]\tNow looping over columns to comput sharing results"))
  c_index <- 1
  slice_df_freqs_melt <- NULL
  for (col in colnames(df)) {
    message(paste("[AP]\t-> processing the colum\t", col, "\tposition", as.character(c_index), "of", as.character(ncol(df)), "\t[", as.character(round(c_index*100/ncol(df),2)), "%]"))
    slice_df <- df[which(df[col]>0),] # slice df
    for (lev in levels) {
      slice_df_lev <- slice_df[which(slice_df[col] == lev),] # slice df
      slice_df_lev_freqs <- getColumnSortedCountsOfLevels(df = slice_df_lev, data_columns = colnames(slice_df_lev), levels_tosearch = levels)
      slice_df_lev_freqs$SourceColumn <- col
      slice_df_lev_freqs$SourceLevel <- lev
      slice_df_lev_freqs$SharedLevel <- rownames(slice_df_lev_freqs)
      slice_df_lev_freqs_melt <- melt(data = slice_df_lev_freqs, id.vars = c("SourceLevel", "SourceColumn", "SharedLevel"), variable.name = "SharedColumn", value.name = "Count")
      slice_df_lev_freqs_melt <- slice_df_lev_freqs_melt[which(slice_df_lev_freqs_melt$Count > 0),]
      # place globally these results
      slice_df_freqs_melt <- rbind(slice_df_freqs_melt, slice_df_lev_freqs_melt)
    }
    c_index <- c_index + 1
  }
  
  return (slice_df_freqs_melt)
}

###############################################################
#' @title Plot streams of flag sharing
#' 
#' @author Andrea Calabria
#' @details version 0.1 (12 August 2019) 
#'
#' @rdname plotStreamsWithAlluvionalsAndBoxplots
#' @docType methods
#' @aliases plotStreamsWithAlluvionalsAndBoxplots
#'
#' @param df an input dataframe of IS matrix.
#' @param data_columns data columns (not annotations)
#' @param plot_filename full path or output plot file name.
#' @param plot_title plot title to use.
#' @param plot_footnote_details footnote to include in the plot; useful for comments.
#' @param height plot height of the putput PDF plot.
#' @param width plot width of the oputput PDF plot.
#' @param seed_number seed
#'
#' @return Plots
#' @usage todo
#' @description Requires specific columns!  "FollowUp", "NIS". If missing, an exception is raised.
#' 
###############################################################
plotStreamsWithAlluvionalsAndBoxplots <- function(df,
                                       data_columns,
                                       color_palette = c("blue4", "blue3", "blue", "deepskyblue", "cyan", "seagreen1", "limegreen", "darkgreen", "chartreuse1", "greenyellow", "gold", "orange", "chocolate1", "red", "darkred", "darkorchid4", "darkorchid1", "pink", "plum", "violet"),
                                       plot_filename = "plot_file.alluvionalbarplot.pdf",
                                       plot_title = "Streams of sharing levels", 
                                       plot_subtitle = "For each tp, the plots presents the source of IS by class (cell differentiation). IS count not unique.",
                                       plot_x_colname = "SharedTimePoint", plot_y_colname = "CountFiltered", plot_fill_colname = "SourceColumn", plot_facet_colname = "LabelClass",
                                       seed_number = 123, height=9, width=12, plot_wrap_ncol = 2) {
  do_plots <- TRUE
  require(ggalluvial)
  required_cols <- c(plot_x_colname, plot_y_colname, plot_fill_colname, plot_facet_colname)
  if (length(grep(FALSE, (required_cols %in% c(colnames(df)))) > 0)) { # some required columns are not there!
    message(paste("[AP]\tERROR: some specified columns are not present in input DF. "))
    do_plots <- FALSE
  }
  
  plot_bars <-
    ggplot(data = df[data_columns], aes_string(x = plot_x_colname, y = plot_y_colname, fill = plot_fill_colname)) +
      geom_bar(stat = "identity", position = "stack") +
      # scale_fill_brewer(type = "qual", palette = "Set3") +
      # scale_color_brewer(type = "qual", palette = "Set3") +
      scale_fill_manual(values = color_palette) + 
      scale_color_manual(values = color_palette) +
      facet_wrap( ~ LabelClass, ncol = plot_wrap_ncol) +
      # geom_text(aes(label=NIS), vjust=1.6, color="black", size=3.5, position = position_stack())+
      scale_x_continuous(breaks = seq(0, max(df[plot_x_colname], na.rm = T), 6) ) +
      theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
      theme(strip.text = element_text(size=16, colour="white")) +
      theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
      labs(list(title = plot_title, x = "Months after gene therapy", y = "N. IS", colour = "Source Time", subtitle = plot_subtitle))
  
  plot_alluvional <- 
    ggplot(data = df[data_columns], aes_string(x = plot_x_colname, y = plot_y_colname, alluvium = plot_fill_colname, fill = plot_fill_colname) ) +
      geom_alluvium(aes_string(fill = plot_fill_colname, colour = plot_fill_colname), alpha = .75, decreasing = FALSE) +
      scale_x_continuous(breaks = seq(0, max(df[plot_x_colname], na.rm = T), 6) ) +
      # scale_fill_brewer(type = "qual", palette = "Set3") +
      # scale_color_brewer(type = "qual", palette = "Set3") +
      scale_fill_manual(values = color_palette) + 
      scale_color_manual(values = color_palette) +
      facet_wrap(~ LabelClass, ncol = plot_wrap_ncol) +
      theme(strip.text.y = element_text(size = 16, colour = "blue", angle = 0)) +
      theme(strip.text = element_text(size=16, colour="white")) +
      theme(strip.background = element_rect(fill="darkblue", colour="white", size=1)) +
      theme(axis.text.x = element_text(size=16), axis.text.y = element_text(size=16), axis.title = element_text(size=16), plot.title = element_text(size=20)) +
      theme(legend.direction = "horizontal", legend.position = "bottom", legend.box = "horizontal") +
      # theme(axis.text.x = element_text(angle = -30, hjust = 0)) +
      labs(list(title = plot_title, x = "Months after gene therapy", y = "N. IS", colour = "Source Time", subtitle = plot_subtitle))
  
  if (do_plots) {
    pdf(file = plot_filename, height=height, width=width)
    print (plot_alluvional)
    print (plot_bars)
    dev.off()
  } # if (do_plots) 
  
}

###############################################################
#' @title Get color codes by time for TIGET list
#' 
#' @author Andrea Calabria
#' @details version 0.1 (12 August 2019) 
#'
#' @rdname getTigetColorSchema_byTimePoint
#' @docType methods
#' @aliases getTigetColorSchema_byTimePoint
#'
#' @param timepoints vector of timepoints (numbers or string)
#' @param unit days or months (default); converted in months
#'
#' @return list of colors (palette)
#' @usage todo
#' @description To have for each time point always the same color code, use this function
#' 
###############################################################
getTigetColorSchema_byTimePoint <- function(timepoints, unit = "months") {
  color_df <- data.frame(
    "TimePoint" = c("00", "01", "02", "03", "06", "09", "12", "18", "24", "30", "36", "42", "48", "54", "60", "66", "72", "78", "84", "90", "96"),
    "ColorCode" = c("gray70", "blue4", "blue3", "blue", "deepskyblue", "cyan", "seagreen1", "limegreen", "darkgreen", "chartreuse1", "greenyellow", "gold", "orange", "chocolate1", "red", "darkred", "darkorchid4", "darkorchid1", "pink", "plum", "violet")
  )
  color_df$TimePoint_int <- as.numeric(as.character(color_df$TimePoint))
  rownames(color_df) <- color_df$TimePoint_int
  
  in_timepoints <- NULL
  if (unit == "days") {
    in_timepoints <- round(as.numeric(as.character(timepoints))/30)
  } else if (unit == "months") {
    in_timepoints <- as.numeric(as.character(timepoints))
  } else {
    message(paste("[AP]\tError: Timepoints unit INVALID:", unit))
  }
  
  out_colors <- as.character(color_df[which(color_df$TimePoint_int %in% in_timepoints), c("ColorCode")])
  
  return (out_colors)
}


###############################################################
#' @title CIS with Grubbs
#' 
#' @author Andrea Calabria
#' @details version 0.1 (20 September 2019) 
#'
#' @rdname CISGrubbs
#' @docType methods
#' @aliases CISGrubbs
#'
#' @param df an input dataframe of IS matrix.
#' @param annotation_cols columns of annotation
#' @param genomic_annotation_genebased_file Annotation file: gene based genomic annotation file. It must be formatted UCSC based, see description for details.
#' @param grubbs_flankingene_bp Extra bp for flanking genes. Default 100000.
#' @param threshold_alpha P-value for filtering significant genes. Default 0.05
#' @param gene_name_col Column name in input df for column gene. Default = "GeneName",
#' @param chr_name_col Column name in input df for column chr. Default = "chr",
#' @param add_standard_padjust Do you want to compute standard p.adjust methods? default = TRUE
#' @param compactDfByRows If the matrix requires be compacted, set this to TRUE (Default = F)
#'
#' @return DF of CIS results
#' @usage todo
#' @description Do CIS analysis with Grubbs test for outliers (as EM and DC designed).
#' 
#' Gene name File formatting:
#'      name2 chrom strand min_txStart max_txEnd minmax_TxLen average_TxLen      name min_cdsStart max_cdsEnd minmax_CdsLen average_CdsLen
# 1     A1BG chr19      -    58858171  58864865         6694          6694 NM_130786     58858387   58864803          6416        6416.00
# 2 A1BG-AS1 chr19      +    58863335  58866549         3214          3214 NR_015380     58866549   58866549             0           0.00
# 3     A1CF chr10      -    52559168  52645435        86267         86267 NM_138933     52566488   52619700         53212       48635.50
# 4      A2M chr12      -     9220303   9268825        48522         48522 NM_000014      9220418    9268445         48027       46276.75
# 5  A2M-AS1 chr12      +     9217772   9220651         2879          2879 NR_026971      9220651    9220651             0           0.00
# 6    A2ML1 chr12      +     8975067   9029377        54310         43044 NM_144670      8975247    9027607         52360       41100.00
#' 
###############################################################
CISGrubbs <- function(df,
                      annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                      gene_name_col = "GeneName",
                      chr_name_col = "chr",
                      genomic_annotation_genebased_file, 
                      grubbs_flankingene_bp = 100000, 
                      threshold_alpha = 0.05,
                      add_standard_padjust = TRUE, 
                      compactDfByRows = FALSE
                      ) {
  
  library(sqldf)
  
  ok_checks_passed <- TRUE
  if (!file.exists(genomic_annotation_genebased_file)) {
    ok_checks_passed <- FALSE
  }
  
  if (ok_checks_passed) {
    if (compactDfByRows) {
      message(paste("[AP]\tCompacting input dataframe (NOTE: be sure that all data columns are real data columns, no 'all' cols accepted."))
      df <- compactDfByRows(df = df, annotation_columns = annotation_cols, data_columns = setdif(colnames(df), annotation_cols))
    }
    
    # annotations
    message(paste("[AP]\tCIS analysis started, files are ok. Importing data."))
    refgenes <- read.csv(file = genomic_annotation_genebased_file, 
                         header=TRUE, fill=T, sep='\t', 
                         check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
    
    df_bygene <- sqldf("select chr, integration_locus, GeneName, GeneStrand, count(distinct integration_locus) as n_IS_perGene, 
                          min(integration_locus) as min_bp_integration_locus, max(integration_locus) as max_bp_integration_locus, (max(integration_locus) - min(integration_locus)) as IS_span_bp, avg(integration_locus) as avg_bp_integration_locus, median(integration_locus) as median_bp_integration_locus, count(distinct strand) as distinct_orientations
                       from df 
                       where 1 
                       group by GeneName, chr")
    df_bygene$chr <- gsub("^", "chr", df_bygene$chr)
    # df_bygene$TotIS_asSumByGene <- sum(df_bygene$n_IS_perGene) # sum of IS
    
    df_bygene_withannotation <- sqldf("select df_bygene.chr, df_bygene.integration_locus, df_bygene.GeneName, df_bygene.n_IS_perGene, refgenes.average_TxLen, df_bygene.min_bp_integration_locus, df_bygene.max_bp_integration_locus, df_bygene.avg_bp_integration_locus, df_bygene.median_bp_integration_locus, df_bygene.distinct_orientations, df_bygene.IS_span_bp
                                      from df_bygene, refgenes
                                      where df_bygene.chr = refgenes.chrom and df_bygene.GeneStrand = refgenes.strand and df_bygene.GeneName = refgenes.name2 ")
    df_bygene_withannotation$TotIS_asDfRow <- nrow(df)
    df_bygene_withannotation$geneIS_frequency_byHitIS <- df_bygene_withannotation$n_IS_perGene / df_bygene_withannotation$TotIS_asDfRow
    
    # add columns as grubbs test requires
    message(paste("[AP]\tComputing tests."))
    extrabp <- grubbs_flankingene_bp
    n_elements <- nrow(df_bygene_withannotation)
    df_bygene_withannotation$raw_gene_integration_frequency <- df_bygene_withannotation$n_IS_perGene / df_bygene_withannotation$average_TxLen
    df_bygene_withannotation$integration_frequency_withtolerance <- (df_bygene_withannotation$n_IS_perGene / (df_bygene_withannotation$average_TxLen + extrabp)) * 1000
    df_bygene_withannotation$minus_log2_integration_freq_withtolerance <- -log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2)
    # average_minus_log2_integration_freq <- mean(df_bygene_withannotation$minus_log2_integration_freq_withtolerance)
    # stdev_minus_log2_integration_freq <- sd(df_bygene_withannotation$minus_log2_integration_freq_withtolerance)
    # df_bygene_withannotation$ratioZEM_minus_log2_integration_freq_withtolerance <- (average_minus_log2_integration_freq - df_bygene_withannotation$minus_log2_integration_freq_withtolerance) / stdev_minus_log2_integration_freq
    
    df_bygene_withannotation$zscore_minus_log2_integration_freq_withtolerance <- scale(-log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2))
    df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance <- -scale(-log(x = df_bygene_withannotation$integration_frequency_withtolerance, base = 2))
    df_bygene_withannotation$t_z_mlif <- sqrt( (n_elements * (n_elements - 2) * (df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance)^2) /
                                                 ( ((n_elements -1)^2) - (n_elements * (df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance)^2))
    )
    
    # df_bygene_withannotation$minus_log2_integration_freq <- -log(x = ((df_bygene_withannotation$raw_gene_integration_frequency)*1000), base = 2)
    # df_bygene_withannotation$z_minus_log_integration_freq <- scale(df_bygene_withannotation$minus_log2_integration_freq)
    # df_bygene_withannotation$t_z_mlif <- sqrt( (nrow(df_bygene_withannotation) * (nrow(df_bygene_withannotation) - 2) * (df_bygene_withannotation$z_minus_log_integration_freq)^2) /
    #                                              ( ((nrow(df_bygene_withannotation) -1)^2) - (nrow(df_bygene_withannotation) * (df_bygene_withannotation$z_minus_log_integration_freq)^2))
    #                                            )
    
    df_bygene_withannotation$tdist2t <- T.DIST.2T(df_bygene_withannotation$t_z_mlif, df = (nrow(df_bygene_withannotation) - 2))
    # df_bygene_withannotation$tdist <- dt(x = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2) * 2
    # df_bygene_withannotation$tdist <- (1 - dt(x = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2))
    df_bygene_withannotation$tdist_pt <- pt(q = df_bygene_withannotation$t_z_mlif, df = nrow(df_bygene_withannotation) - 2)
    
    df_bygene_withannotation$tdist_bonferroni_AC <- ifelse(df_bygene_withannotation$tdist2t * nrow(df_bygene_withannotation) > 1, 1, df_bygene_withannotation$tdist2t * nrow(df_bygene_withannotation))
    if (add_standard_padjust) {
      df_bygene_withannotation$tdist_bonferroni <- p.adjust(df_bygene_withannotation$tdist2t, method = "bonferroni", n = length(df_bygene_withannotation$tdist2t))
      df_bygene_withannotation$tdist_fdr <- p.adjust(df_bygene_withannotation$tdist2t, method = "fdr", n = length(df_bygene_withannotation$tdist2t))
      df_bygene_withannotation$tdist_benjamini <- p.adjust(df_bygene_withannotation$tdist2t, method = "BY", n = length(df_bygene_withannotation$tdist2t))
    }
    
    df_bygene_withannotation$tdist_positive_and_corrected <- ifelse((df_bygene_withannotation$tdist_bonferroni_AC < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                                    df_bygene_withannotation$tdist_bonferroni_AC, 
                                                                    NA)
    df_bygene_withannotation$tdist_positive <- ifelse((df_bygene_withannotation$tdist2t < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                      df_bygene_withannotation$tdist2t, 
                                                      NA)
    EM_correction_N <- length(df_bygene_withannotation$tdist_positive[!is.na(df_bygene_withannotation$tdist_positive)])
    df_bygene_withannotation$tdist_positive_and_correctedEM <- ifelse((df_bygene_withannotation$tdist2t * EM_correction_N < threshold_alpha & df_bygene_withannotation$neg_zscore_minus_log2_integration_freq_withtolerance > 0), 
                                                                       df_bygene_withannotation$tdist2t * EM_correction_N, 
                                                                       NA)
    
    return (df_bygene_withannotation)
    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files."))
  }  # if (ok_checks_passed)
  
}

###############################################################
#' @title MS Excel like T.DIST.2T
#' 
#' @author Andrea Calabria
#' @details version 0.1 (20 September 2019) 
#'
#' @rdname T.DIST.2T
#' @docType methods
#' @aliases T.DIST.2T
#'
#' @param x an input vector.
#' @param df degrees of freedom
#'
#' @return T.DIST.2T
#' @usage todo
#' @description ToDo
#' 
###############################################################

T.DIST.2T <- function(x, df) {
  return ( (1 - pt(x, df))*2 )
}


###############################################################
#' @title Aggregate Metadata
#' 
#' @author Andrea Calabria
#' @details version 0.1 (23 September 2019) 
#'
#' @rdname aggregateMetadata
#' @docType methods
#' @aliases aggregateMetadata
#'
#' @param af_df Association file (AF) as df.
#' @param groupby_fields Fields of AF to br grouped by in the SQL statement. Default = c("SubjectID", "CellMarker", "Tissue", "TimePoint").
#' @param timepoint_field Among the groupby_fields, select the time point. Default c("TimePoint"). If present, it will arrange this field as number to convert/pad. 
#' @param timepoint_field_to_string String col name for time point as converted in string. Default = "Followup_StrPad",
#' @param timepoint_field_to_int String col name for time point as converted in integer (rounded integer) Default = "Followup_Int",
#' @param groupby_outcolname Column name of the output column for aggregation of the groupby_fields. Default = "groupby_query"
#' @param colname_to_string Column name of the output column for timepoint conversion to months (as string). Default = FollowupMonths_StrPad
#' @param colname_to_string Column name of the output column for timepoint conversion to months (as string). Default = FollowupMonths_StrPad
#' 
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

aggregateMetadata <- function(af_df, 
                              groupby_fields = c("SubjectID", "CellMarker", "Tissue", "TimePoint"),
                              timepoint_field = c("TimePoint"),
                              timepoint_field_to_string = "Followup_StrPad",
                              timepoint_field_to_int = "Followup_Int",
                              groupby_outcolname = "groupby_query",
                              where_string = "1",
                              additional_select_querystring = "",
                              padspace = 2
) {
  require(sqldf)
  min_set_columns <- c("FusionPrimerPCRDate", "LinearPCRDate", "VCN", "DNAngUsed", "Kapa", "ulForPool")
  error <- ""
  ok_checks_passed <- TRUE
  if (length(setdiff(groupby_fields, colnames(af_df))) > 0) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Some columns in your groupby_fields are not present in AF columns.")
  } 
  if ( "" %in% colnames(af_df) ) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> You have an empty column in the AF dataframe. Remove it before continuing.")
  }  
  if (length(setdiff(min_set_columns, colnames(af_df))) > 0) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Some columns of the MINIMAL COLUMN SETis missing in input AF.\n\t\t\tMin set: ", paste(min_set_columns, collapse = ","), " # Missing set setdiff:", paste(setdiff(min_set_columns, colnames(af_df)), collapse = ","))
  }
  if (!(timepoint_field %in% groupby_fields)) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Timepoint in groupby_fields not exists as colname AF.")
  }  
  if (!(groupby_outcolname == "groupby_query")) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> groupby_outcolname = groupby_query, please change it to your grouping values.")
  }
  
  if (ok_checks_passed) {
    message(paste("[AP]\tAggregating metadata."))
    
    # groupby_fields <- c("SubjectID", "CellMarker", "Tissue", "TimePoint") # SubjectID, CellMarker, Tissue, TimePoint
    groupby_fields_string <- paste(groupby_fields, collapse = ", ")
    label_metadata_agg_sql_query <- paste0("SELECT ", groupby_fields_string, additional_select_querystring, ", 
                                           min(FusionPrimerPCRDate) as FusionPrimerPCRDate, min(LinearPCRDate) as LinearPCRDate, 
                                           AVG(VCN) as VCN, AVG(DNAngUsed) as Avg_DNAngUsed, AVG(Kapa) as Kapa, 
                                           sum(DNAngUsed) as DNAngUsed, sum(ulForPool) as ulForPool 
                                           FROM af_df
                                           WHERE ", where_string, "  
                                           GROUP BY ", groupby_fields_string, " 
                                           ORDER BY ", groupby_fields_string)
    # run SQL
    label_metadata_agg_sql_query_result <- sqldf(label_metadata_agg_sql_query)
    
    if (length(timepoint_field) > 0 ) {
      label_metadata_agg_sql_query_result$Followup_StrPad <- str_pad(round(as.numeric(label_metadata_agg_sql_query_result[,timepoint_field])), padspace, pad = "0")
      names(label_metadata_agg_sql_query_result)[names(label_metadata_agg_sql_query_result) == "Followup_StrPad"] <- timepoint_field_to_string
      label_metadata_agg_sql_query_result$Followup_StrPad <- str_pad(round(as.numeric(label_metadata_agg_sql_query_result[,timepoint_field])), padspace, pad = "0")
      names(label_metadata_agg_sql_query_result)[names(label_metadata_agg_sql_query_result) == "Followup_Int"] <- timepoint_field_to_int
      
      groupby_fields_withtstring <- gsub(timepoint_field, timepoint_field_to_string, groupby_fields)
      label_metadata_agg_sql_query_result$groupby_query <- apply(label_metadata_agg_sql_query_result[groupby_fields_withtstring], 1, paste, collapse = '_')
      rownames(label_metadata_agg_sql_query_result) <- label_metadata_agg_sql_query_result$groupby_query
      names(label_metadata_agg_sql_query_result)[names(label_metadata_agg_sql_query_result) == "groupby_query"] <- groupby_outcolname
      
    } else {
      label_metadata_agg_sql_query_result$groupby_query <- apply(label_metadata_agg_sql_query_result[groupby_fields], 1, paste, collapse = '_')
      rownames(label_metadata_agg_sql_query_result) <- label_metadata_agg_sql_query_result$groupby_query
      names(label_metadata_agg_sql_query_result)[names(label_metadata_agg_sql_query_result) == "groupby_query"] <- groupby_outcolname
    }
    
    return (label_metadata_agg_sql_query_result)
  } else {
    message(paste("[AP]\tERROR: Some problems with input files:", error))
    return (NULL)
  }  # if (ok_checks_passed)
  
}



###############################################################
#' @title Load gene annotated as onco or tumor suppressors
#' 
#' @author Andrea Calabria
#' @details version 0.1 (27 September 2019) 
#'
#' @rdname loadOncoTSgenes
#' @docType methods
#' @aliases loadOncoTSgenes
#'
#' @param onco_db_file Default = "/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/201806_uniprot-Proto-oncogene.tsv", 
#' @param tumsup_db_file Default = "/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/201806_uniprot-Tumor-suppressor.tsv",
#' @param species Default = "all" # alternatives: human, mouse
#' 
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

loadOncoTSgenes <- function(onco_db_file = "/Users/calabria.andrea/Dropbox (HSR Global)/TIGET/Workbench/isatk/script/R/publicdb/201806_uniprot-Proto-oncogene.tsv", 
                            tumsup_db_file = "/Users/calabria.andrea/Dropbox (HSR Global)/TIGET/Workbench/isatk/script/R/publicdb/201806_uniprot-Tumor-suppressor.tsv",
                            species = "all" # alternatives: human, mouse
                            ) {
  require(Hmisc)
  error <- ""
  ok_checks_passed <- TRUE
  if (!file.exists(onco_db_file) | !file.exists(tumsup_db_file)) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Input files are not existing or wrong path.")
  }
  if (!(species %in% c("all", "human", "mouse"))) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Specie not existing, choose one among: all, human, mouse")
  }
  
  if (ok_checks_passed) {
    message(paste("[AP]\tLoading annotated genes"))
    
    # acquire DB
    onco_db <- read.csv(file = onco_db_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    tumsup_db <- read.csv(file = tumsup_db_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    # onco_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/uniprot-oncogene_Mouse.tab", header=TRUE, fill=T, sep='\t', check.names = FALSE)
    # tumsup_db <- read.csv("/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/uniprot-tumorsuppressor_Mouse.tab", header=TRUE, fill=T, sep='\t', check.names = FALSE)
    oncots_df_touse <- NULL # output df
    if (species == "mouse") {
      specie <- "Mus musculus (Mouse)"
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Title case"))
      # subset data for this specie
      mouse_onco_db <- onco_db[which(onco_db$Organism == specie),]
      # mouse_onco_db <- onco_db
      # rownames(mouse_onco_db) <- mouse_onco_db$`Gene names  (primary )`
      mouse_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      # mouse_tumsup_db <- tumsup_db
      # rownames(mouse_tumsup_db) <- mouse_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      mouse_onco_db_genes <- unique(as.character(mouse_onco_db$`Gene names  (primary )`))
      mouse_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(mouse_onco_db[which(mouse_onco_db$Status == "reviewed" & !is.na(mouse_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      mouse_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(mouse_onco_db_genes_allnames)))) # merge datasets
      mouse_onco_db_genes_allnames_df <- data.frame("OncoGene" = mouse_onco_db_genes_allnames)
      rownames(mouse_onco_db_genes_allnames_df) <- mouse_onco_db_genes_allnames_df$OncoGene
      
      mouse_tumsup_db_genes <- unique(as.character(mouse_tumsup_db$`Gene names  (primary )`))
      mouse_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(mouse_tumsup_db[which(mouse_tumsup_db$Status == "reviewed" & !is.na(mouse_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      mouse_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(mouse_tumsup_db_genes_allnames))) ) # merge datasets
      mouse_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = mouse_tumsup_db_genes_allnames)
      rownames(mouse_tumsup_db_genes_allnames_df) <- mouse_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      mouse_oncotumsup_db_genes_allnames_df <- merge(x = mouse_onco_db_genes_allnames_df, y = mouse_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(mouse_oncotumsup_db_genes_allnames_df) <- mouse_oncotumsup_db_genes_allnames_df$Row.names
      names(mouse_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      mouse_oncotumsup_db_genes_allnames_df <- cbind(mouse_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(mouse_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3,
                  ifelse(.onco > 0, 1,
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      oncots_df_touse <- mouse_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
    }
    
    if (species == "human") {
      specie <- "Homo sapiens (Human)"
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Upper case."))
      # subset data for this specie
      human_onco_db <- onco_db[which(onco_db$Organism == specie),]
      # human_onco_db <- onco_db
      # rownames(human_onco_db) <- human_onco_db$`Gene names  (primary )`
      human_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      # human_tumsup_db <- tumsup_db
      # rownames(human_tumsup_db) <- human_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      human_onco_db_genes <- unique(as.character(human_onco_db$`Gene names  (primary )`))
      human_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(human_onco_db[which(human_onco_db$Status == "reviewed" & !is.na(human_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      # human_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(human_onco_db_genes_allnames)))) # merge datasets
      human_onco_db_genes_allnames <- gsub(";", "", unique(human_onco_db_genes_allnames)) # merge datasets
      human_onco_db_genes_allnames_df <- data.frame("OncoGene" = human_onco_db_genes_allnames)
      rownames(human_onco_db_genes_allnames_df) <- human_onco_db_genes_allnames_df$OncoGene
      
      human_tumsup_db_genes <- unique(as.character(human_tumsup_db$`Gene names  (primary )`))
      human_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(human_tumsup_db[which(human_tumsup_db$Status == "reviewed" & !is.na(human_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      # human_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(human_tumsup_db_genes_allnames))) ) # merge datasets
      human_tumsup_db_genes_allnames <- gsub(";", "", unique(human_tumsup_db_genes_allnames)) # merge datasets
      human_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = human_tumsup_db_genes_allnames)
      rownames(human_tumsup_db_genes_allnames_df) <- human_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      human_oncotumsup_db_genes_allnames_df <- merge(x = human_onco_db_genes_allnames_df, y = human_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(human_oncotumsup_db_genes_allnames_df) <- human_oncotumsup_db_genes_allnames_df$Row.names
      names(human_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      human_oncotumsup_db_genes_allnames_df <- cbind(human_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(human_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3, 
                  ifelse(.onco > 0, 1, 
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      oncots_df_touse <- human_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
    }
    
    if (species == "all") {
      message(paste("[AP]\t-> Specie selected:", species, " -> output gene names will be Title case"))
      # subset data for this specie
      # allspecies_onco_db <- onco_db[which(onco_db$Organism == specie),]
      allspecies_onco_db <- onco_db
      # rownames(allspecies_onco_db) <- allspecies_onco_db$`Gene names  (primary )`
      # allspecies_tumsup_db <- tumsup_db[which(tumsup_db$Organism == specie),]
      allspecies_tumsup_db <- tumsup_db
      # rownames(allspecies_tumsup_db) <- allspecies_tumsup_db$`Gene names  (primary )`
      
      # get gene list: only reviewd genes and all gene aliases (from all species, thus for mice you must do lowe case and capitalize the gene name)
      allspecies_onco_db_genes <- unique(as.character(allspecies_onco_db$`Gene names  (primary )`))
      allspecies_onco_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(allspecies_onco_db[which(allspecies_onco_db$Status == "reviewed" & !is.na(allspecies_onco_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      allspecies_onco_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(allspecies_onco_db_genes_allnames)))) # merge datasets
      allspecies_onco_db_genes_allnames_df <- data.frame("OncoGene" = allspecies_onco_db_genes_allnames)
      rownames(allspecies_onco_db_genes_allnames_df) <- allspecies_onco_db_genes_allnames_df$OncoGene
      
      allspecies_tumsup_db_genes <- unique(as.character(allspecies_tumsup_db$`Gene names  (primary )`))
      allspecies_tumsup_db_genes_allnames <- unique(unlist(strsplit(as.character(unique(as.character(allspecies_tumsup_db[which(allspecies_tumsup_db$Status == "reviewed" & !is.na(allspecies_tumsup_db$`Gene names`)), c("Gene names")]))), " ", fixed = TRUE), function(x) {c(x)}))
      allspecies_tumsup_db_genes_allnames <- gsub(";", "", unique(capitalize(tolower(allspecies_tumsup_db_genes_allnames))) ) # merge datasets
      allspecies_tumsup_db_genes_allnames_df <- data.frame("TumorSuppressor" = allspecies_tumsup_db_genes_allnames)
      rownames(allspecies_tumsup_db_genes_allnames_df) <- allspecies_tumsup_db_genes_allnames_df$TumorSuppressor
      
      # merge df
      allspecies_oncotumsup_db_genes_allnames_df <- merge(x = allspecies_onco_db_genes_allnames_df, y = allspecies_tumsup_db_genes_allnames_df, by = 0, all = T)
      rownames(allspecies_oncotumsup_db_genes_allnames_df) <- allspecies_oncotumsup_db_genes_allnames_df$Row.names
      names(allspecies_oncotumsup_db_genes_allnames_df) <- c("GeneName", "OncoGene", "TumorSuppressor")
      allspecies_oncotumsup_db_genes_allnames_df <- cbind(allspecies_oncotumsup_db_genes_allnames_df, data.frame(
        "Onco1_TS2" = apply(allspecies_oncotumsup_db_genes_allnames_df, 1, function(x) {
          .onco <- ifelse(!(is.na(x[2])), 1, 0)
          .tums <- ifelse(!(is.na(x[3])), 1, 0)
          ifelse( (.onco > 0 & .tums > 0 ), 3, 
                  ifelse(.onco > 0, 1, 
                         ifelse(.tums > 0, 2, NA)
                  )
          )
        })
      )
      )
      
      oncots_df_touse <- allspecies_oncotumsup_db_genes_allnames_df
      return (oncots_df_touse)
      
    } # if all
    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files:", error))
    return (NULL)
  }  # if (ok_checks_passed)
  
}


###############################################################
#' @title Load gene annotated as onco or tumor suppressors
#' 
#' @author Andrea Calabria
#' @details version 0.1 (14 October 2019) 
#'
#' @rdname loadOncoTSgenes_fromCancerMine
#' @docType methods
#' @aliases loadOncoTSgenes_fromCancerMine
#'
#' @param cancermine_collated_file Default = "/Users/calabria.andrea/Dropbox (FONDAZIONE TELETHON)/Project Tumour Development/source/publicdb/201806_uniprot-Proto-oncogene.tsv", 
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

loadOncoTSgenes_fromCancerMine <- function(cancermine_collated_file = "/Users/calabria.andrea/Dropbox (HSR Global)/TIGET/Workbench/isatk/script/R/publicdb/cancermine/cancermine_collated.201910.tsv") {
  require(Hmisc)
  error <- ""
  ok_checks_passed <- TRUE
  if (!file.exists(cancermine_collated_file)) {
    ok_checks_passed <- FALSE
    error <- paste(error, "\n[AP]\t\t\t-> Input files are not existing or wrong path.")
  }
  
  if (ok_checks_passed) {
    message(paste("[AP]\tLoading annotated genes"))
    
    # acquire DB
    cancermine_db <- read.csv(file = cancermine_collated_file, header=TRUE, fill=T, sep='\t', check.names = FALSE)
    
    cancermine_db$GeneName <- cancermine_db$gene_normalized
    cancermine_db$OncoTS <- cancermine_db$role
    cancermine_db$OncoTS <- gsub("Oncogene", "OncoGene", gsub("_", "", cancermine_db$OncoTS))
    
    return (cancermine_db)

    
  } else {
    message(paste("[AP]\tERROR: Some problems with input files:", error))
    return (NULL)
  }  # if (ok_checks_passed)
  
}


###############################################################
#' @title findLowerOutliersByColumns_byZScore
#' 
#' @author Andrea Calabria
#' @details version 0.1 (September 2019) 
#'
#' @rdname findLowerOutliersByColumns_byZScore
#' @docType methods
#' @aliases findLowerOutliersByColumns_byZScore
#'
#' @param df
#'
#' @return dataframe
#' @usage todo
#' @description ToDo
#' 
###############################################################

findLowerOutliersByColumns_byZScore <- function(df, 
                                                data_columns,
                                                n_sigma_threshold = 1.5 ) {
  # 1. do numbers of elements by columns and normalize values
  # 2. find elements boyond N sigma
  
  this_df <- df[data_columns]
  this_df[is.na(this_df)] <- 0
  df_stats <- data.frame("N_IS" = apply(this_df, 2, function(x) {length(x[x>0])} ))
  df_stats$zscore <- scale(df_stats$N_IS)
  # df_stats_desc <- describe(df_stats$N_IS)
  
  lower_outliers <- rownames(df_stats[which(df_stats$zscore < (-1*n_sigma_threshold)),])
  if (length(lower_outliers) > 0) {
    return (lower_outliers)
  } else {
    return (c())
  }
  
}


###############################################################
#' @title Relabel an input string df of flag frequencies based on transition matrix
#' 
#' @author Andrea Calabria
#' @details version 0.1.1 bug fixed (24 April 2020), source 0.1 (29 November 2019)
#'
#' @rdname relabelFlagString_byTransitionMatrix
#' @docType methods
#' @aliases relabelFlagString_byTransitionMatrix
#'
#' @param df_singlerow an input dataframe of occurrences by label: colnames are class labels (flags), values are frequences.
#' @param pairs_transition_matrix df reporting transition labels from (rows) to (columns). No default.
#' @param skip_root To avoid biases in transition evaluation, skipt root; this is a trick but it works because a root is never changing any other leaves, thus if not impacting the step of evaluation (if the matrix of trnasitions is properly set). Default = T.
#' @param quiet Avoid verbosity in output. Default = T.
#'
#' @return a dataframe or a matrix of occurrences
#' @description Similar to table(df) but works faster and by rows.
#' 
###############################################################
relabelFlagString_byTransitionMatrix <- function(df_singlerow, pairs_transition_matrix, skip_root = T, root_id = "1", quiet = T) {
  if (!quiet) {message(paste("[AP]\tRelabel Flags by Transition matrix.", "input string: freqs =", paste(df_singlerow, collapse = ' '), "- classes =", paste(colnames(df_singlerow), collapse = ' ')))}
  
  # Step 1: SORTING: sort input df by frequency
  df_singlerow_sorted <- sort(df_singlerow, decreasing = T) # sorted but now remove 0's
  df_singlerow_sorted <- df_singlerow_sorted[,df_singlerow_sorted>0]
  relabelingmap <- data.frame("relabeled" = colnames(df_singlerow_sorted), stringsAsFactors=FALSE) # the output
  rownames(relabelingmap) <- colnames(df_singlerow_sorted)
  df_singlerow_sorted_relabeled <- NULL
  
  # in case of multiple values at the same identical (max) frequency
  occurrences_of_max <- length(df_singlerow_sorted[df_singlerow_sorted==max(df_singlerow_sorted)]) # max(df_singlerow_sorted)
  do_all_comparisons_for_max_value <- F
  if (occurrences_of_max > 2) {
    do_all_comparisons_for_max_value <- T
  }
  
  # Step 2: EVALUATION
  # 2.1. definisci tutti i matches da sorted
  df_singlerow_sorted_pairs_ij <- sortedPairsOfAllCombinations_fromVector(colnames(df_singlerow_sorted), forwardOnly = T)
  todo_list <- colnames(df_singlerow_sorted)
  if (skip_root) {
    todo_list <- setdiff(todo_list, root_id)
  }
  already_done <- c()
  tree_root <- NULL # this variable is instantiatied if a pair of classess has the same frequences, thus their root is defined such that eventual next comparisons will inherit the first comparison output (the root)
  if (length(todo_list)>1) {
    while (length(todo_list) > 0) {
      i <- todo_list[1]
      todo_list <- setdiff(todo_list, i) # now remove i from todolist
      already_done <- union(already_done, i)
      if (!quiet) {message(paste("[AP]\tdebug:\tAnalyze i =", i))}
      # define the list of comparisons to inspect
      islice_df_singlerow_sorted_pairs_ij <- df_singlerow_sorted_pairs_ij[which(df_singlerow_sorted_pairs_ij$From == i),]
      rownames(islice_df_singlerow_sorted_pairs_ij) <- islice_df_singlerow_sorted_pairs_ij$To
      ni <- as.numeric(as.character(df_singlerow_sorted[i]))
      for (k in rownames(islice_df_singlerow_sorted_pairs_ij)) { # here you are keeping the sorting! else it would be enough to list all levels(factor("To))
        # j <- as.character(islice_df_singlerow_sorted_pairs_ij[k, "To"])
        j <- k
        nj <- as.numeric(as.character(df_singlerow_sorted[,j]))
        if (!quiet) {message(paste("[AP]\tdebug:\t\tj =", j))}
        # 2.2 do comparisons of the frequencies and evaluate
        # you can have ONLY 2 cases: > or = -> those will be your branches
        if (ni > nj) { # branch of ">"
          # first thing: initialize the output frquencies with this max class
          if (length(df_singlerow_sorted_relabeled) == 0) {
            df_singlerow_sorted_relabeled <- df_singlerow_sorted[1]
          } 
          # then evaluate...
          change_j_to_this <- pairs_transition_matrix[which(pairs_transition_matrix$From == i & pairs_transition_matrix$To == j), "ConvertTo"]
          relabelingmap[j, c("relabeled")] <- change_j_to_this
          if (!quiet) {message(paste("[AP]\tdebug:\tchange j to =", change_j_to_this)); print(df_singlerow_sorted_relabeled)}
          # Step 3: UPDATE - Recursion
          # update all data!
          # first update the frquencies with labels
          if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) { # if the new col exists, copy freqs, else add value
            df_singlerow_sorted_relabeled[change_j_to_this] <- df_singlerow_sorted_relabeled[change_j_to_this] + nj
          } else {
            df_singlerow_sorted_relabeled[change_j_to_this] <- nj
          } # if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) 
          # then update loop variables
          todo_list <- setdiff(todo_list, j)
          # TODO other updates?
        } else {
          if (ni == nj) { # branch of "=="
            change_j_to_this <- pairs_transition_matrix[which(pairs_transition_matrix$From == i & pairs_transition_matrix$To == j), "ConvertTo"]
            change_i_to_this <- pairs_transition_matrix[which(pairs_transition_matrix$From == j & pairs_transition_matrix$To == i), "ConvertTo"]
            relabelingmap[j, c("relabeled")] <- change_j_to_this
            relabelingmap[i, c("relabeled")] <- change_i_to_this
            if (!quiet) {message(paste("[AP]\tdebug:\tchange j to =", change_j_to_this)); }
            if (!quiet) {message(paste("[AP]\tdebug:\tchange i to =", change_i_to_this)); }
            
            # Step 3: UPDATE - Recursion
            # update all data!
            # initialize the output frquencies with the first class if not yet existing
            if (length(df_singlerow_sorted_relabeled) == 0) {
              df_singlerow_sorted_relabeled <- df_singlerow_sorted[1]
              df_singlerow_sorted_relabeled[1] <- 0
              # first update the frquencies with labels of j
              if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) { # if the new col exists, copy freqs, else add value
                df_singlerow_sorted_relabeled[change_j_to_this] <- df_singlerow_sorted_relabeled[change_j_to_this] + nj
              } else {
                df_singlerow_sorted_relabeled[change_j_to_this] <- nj
              } # if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) 
              # then update loop variables
              # if (!do_all_comparisons_for_max_value) { # if you need to do all comparisons, then DO NOT remov ehte comparison from this list
                todo_list <- setdiff(todo_list, j) 
                # }
              # first update the frquencies with labels of i
              if (change_i_to_this %in% colnames(df_singlerow_sorted_relabeled)) { # if the new col exists, copy freqs, else add value
                df_singlerow_sorted_relabeled[change_i_to_this] <- df_singlerow_sorted_relabeled[change_i_to_this] + ni
              } else {
                df_singlerow_sorted_relabeled[change_i_to_this] <- ni
              } # if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) 
              # then update loop variables
              todo_list <- setdiff(todo_list, i)  
              
            } else {
              # first update the frquencies with labels of j
              if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) { # if the new col exists, copy freqs, else add value
                df_singlerow_sorted_relabeled[change_j_to_this] <- df_singlerow_sorted_relabeled[change_j_to_this] + nj
              } else {
                df_singlerow_sorted_relabeled[change_j_to_this] <- nj
              } # if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) 
              # then update loop variables
              todo_list <- setdiff(todo_list, j)
              # first update the frquencies with labels of i, only if not already done
              if(!(i %in% already_done)) {
                if (change_i_to_this %in% colnames(df_singlerow_sorted_relabeled)) { # if the new col exists, copy freqs, else add value
                  df_singlerow_sorted_relabeled[change_i_to_this] <- df_singlerow_sorted_relabeled[change_i_to_this] + ni
                } else {
                  df_singlerow_sorted_relabeled[change_i_to_this] <- ni
                } # if (change_j_to_this %in% colnames(df_singlerow_sorted_relabeled)) 
                # then update loop variables
                todo_list <- setdiff(todo_list, i)
              } # if(!(i %in% already_done))
            } # if (ncol(df_singlerow_sorted_relabeled) == 0) 
            
            if (!quiet) {message(paste("[AP]\tdebug:\tAfter both updates:")); print(df_singlerow_sorted_relabeled)}
            
            # just as a tree update, not used yet
            tree_root <- min(as.numeric(change_i_to_this), as.numeric(change_j_to_this))
            
          } else { # IMPOSSIBLE CASES HERE!
            message(paste("[AP]\tERROR!!! A case NOT considered appeared here! ID row:\t", rownames(df_singlerow)))
            print (df_singlerow)
          } # if (ni == nj)
        } # if (ni > nj)

        # update done list already_done
        already_done <- union(already_done, j)
      } # for (k in rownames(islice_df_singlerow_sorted_pairs_ij))
    } #  while (length(todo_list)
    
    # map input col data into new data fields -> extend df
    input_col_order <- colnames(df_singlerow)
    missing_t_in_relabeled <- df_singlerow[setdiff(input_col_order, colnames(df_singlerow_sorted_relabeled))]
    missing_t_in_relabeled[missing_t_in_relabeled>0] <- 0
    df_singlerow_sorted_relabeled <- cbind(df_singlerow_sorted_relabeled, missing_t_in_relabeled)[input_col_order]
    
    return (list("frequencies" = df_singlerow_sorted_relabeled, 
                 "relabling_map" = relabelingmap,
                 "relabling_map_from" = as.numeric(rownames(relabelingmap)), #as.numeric(rownames(relabeling_row$relabling_map))
                 "relabling_map_to" = as.numeric(as.character(relabelingmap$relabeled))
                )
            )
  } else {
    return (list("frequencies" = df_singlerow, 
                 "relabling_map" = relabelingmap,
                 "relabling_map_from" = as.numeric(rownames(relabelingmap)), #as.numeric(rownames(relabeling_row$relabling_map))
                 "relabling_map_to" = as.numeric(as.character(relabelingmap$relabeled))
                )
            )
  }  # if (length(todo_list)>1)
  
}

###############################################################
#' @title Relabel an input string df of flag frequencies based on transition matrix
#' 
#' @author Andrea Calabria
#' @details version 0.1 (29 November 2019) 
#'
#' @rdname pairsFrom_TransitionMatrix
#' @docType methods
#' @aliases pairsFrom_TransitionMatrix
#'
#' @param df_singlerow an input dataframe of occurrences by label: colnames are class labels (flags), values are frequences.
#' @param transition_matrix df reporting transition labels from (rows) to (columns). No default.
#'
#' @return a dataframe or a matrix of occurrences
#' @description Similar to table(df) but works faster and by rows.
#' 
###############################################################
pairsFrom_TransitionMatrix <- function(df, all_as_chars = T, varnames = c("From", "To"), valuename = "ConvertTo", quiet = T) {
  if (!quiet) {message(paste("[AP]\tDo pais from the input Transition matrix."))}
  pairsTM <- as.data.frame(melt(as.matrix(df), varnames = varnames, value.name = valuename, na.rm = F))
  if (all_as_chars) {
    pairsTM$From <- as.character(pairsTM$From)
    pairsTM$To <- as.character(pairsTM$To)
    pairsTM$ConvertTo <- as.character(pairsTM$ConvertTo)
  }
  return (pairsTM)
}

###############################################################
#' @title Get all combinations of pairs from string (with no repetition and from sorted by input)
#' 
#' @author Andrea Calabria
#' @details version 0.1 (29 November 2019) 
#'
#' @rdname sortedPairsOfAllCombinations_fromVector
#' @docType methods
#' @aliases sortedPairsOfAllCombinations_fromVector
#'
#' @param invector vector of string
#' @param forwardOnly If forwardOnly == T, return only the combinations from the first element on, no repetitions of the first element as "To".
#'
#' @return a dataframe or a matrix of occurrences
#' @description Similar to table(df) but works faster and by rows.
#' 
###############################################################
sortedPairsOfAllCombinations_fromVector <- function(invector, forwardOnly = T, quiet = T) {
  if (!quiet) {message(paste("[AP]\tDo SORTED (FWD or not) of pais from vector string."))}
  allPairs <- NULL
  bag_of_input_vector <- invector # you will erose this variable
  if (forwardOnly) {
    for (i in bag_of_input_vector) {
      bag_of_input_vector <- bag_of_input_vector[! bag_of_input_vector %in% i]
      for (j in bag_of_input_vector) {
        if (i != j) {
          allPairs <- rbind(allPairs, data.frame("From" = i, "To" = j))
        }
      }
    }
  } else {
    for (i in bag_of_input_vector) {
      # bag_of_input_vector <- bag_of_input_vector[! bag_of_input_vector %in% i]
      for (j in bag_of_input_vector) {
        if (i != j) {
          allPairs <- rbind(allPairs, data.frame("From" = i, "To" = j))
        }
      }
    }
  }

  return (allPairs)
}

###############################################################
#' @title Pad occurrence DF with missing columns
#' 
#' @author Andrea Calabria
#' @details version 0.1 (24 April 2020) 
#'
#' @rdname padMissingColsOccurrences
#' @docType methods
#' @aliases padMissingColsOccurrences
#'
#' @param df_occurrences df of occurrences with colnames as string of numbers (not padded)
#' @param min_numer min number of col names to look for padding. Default: NA
#' @param max_numer man number of col names to look for padding. Default: NA
#'
#' @return a dataframe of padded occurrences (dim: rows as input df, cols >= input df)
#' @description To use for relabelFlagString_byTransitionMatrix Usage example: flag_matrix_relabeled_tracking <- padMissingColsOccurrences(df_occurrences = flag_matrix_relabeled_tracking, min_numer = min(labeling_to), max_numer = max(labeling_to))
#' 
###############################################################
padMissingColsOccurrences <- function(df_occurrences, min_numer = NA, max_numer = NA) {
  message(paste("[AP]\tPadding df with missing numerated columns. From cols:", paste(colnames(df_occurrences), collapse = ','), " to padded missing cols:"))
  actual_numeric_cols <- as.numeric(as.character(colnames(df_occurrences)))
  if (is.na(min_numer)) { min_numer <- min(actual_numeric_cols)}
  if (is.na(max_numer)) { max_numer <- max(actual_numeric_cols)}
  missing_numeric_cols <- sort(unique(setdiff(x = seq(min_numer, max_numer), y = actual_numeric_cols)))
  message(paste("[AP]\t\t", paste(missing_numeric_cols, collapse =' - ')))
  # now add pad cols to input df
  out_col_order <- as.character(sort(c(actual_numeric_cols, missing_numeric_cols)))
  empty_df_to_add <- as.data.frame(matrix(ncol = length(missing_numeric_cols), nrow = nrow(df_occurrences)))
  names(empty_df_to_add) <- as.character(missing_numeric_cols)
  empty_df_to_add[is.na(empty_df_to_add)] <- 0
  df_out <- cbind(df_occurrences, empty_df_to_add)
  return (df_out[out_col_order])
}


###############################################################
#' @title Read GREAT file of TSS distance by IS
#' 
#' @author Andrea Calabria
#' @details version 0.1 (used in data filtering), 2020-05-19
#'
#' @rdname readGreatISoutput
#' @docType methods
#' @aliases readGreatISoutput
#'
#' @param great_filename File name of GREAT output file for IS distance.
#' @param sep separator to read the file (tab).
#' @param skip number of rows to skip because comments (1).
#'
#' @return a df of IS annotated
#' @usage TODO
#' @note : 
#'
###############################################################
readGreatISoutput <- function(great_filename,
                              sep = '\t',
                              skip = 1
) {
  message(paste("[AP]\tAcquiring GREAT file of IS distance to genes, input file:",  great_filename))
  great_df <- read.csv(file = great_filename, header=F, fill=T, sep=sep, check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""), skip = skip)
  names(great_df) <- c("ISname", "GreatNote")
  great_df$GreatNote <- gsub("\\(|\\)", "", great_df$GreatNote)
  
  great_notes_df <- cbind( as.data.frame(t(as.data.frame(lapply(strsplit(as.character(great_df$GreatNote), " ", fixed = TRUE), function(x) {c(x)[1:2]})))) )
  names(great_notes_df) <- c("GeneName", "ISdistance_String")
  great_notes_df$ISdistance <- as.numeric(as.character(great_notes_df$ISdistance_String))
  great_notes_df <- cbind(great_notes_df, great_df)
  great_notes_df <- great_notes_df[which(!is.na(great_notes_df$ISdistance_String)),]
  
  return (great_notes_df)
}

###############################################################
#' @title Compute sharing between a source and a output
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-05-20
#'
#' @rdname sharingSourceToOutput_unsorted
#' @docType methods
#' @aliases readGreatISoutput
#'
#' @param df File name of GREAT output file for IS distance.
#' @param sep separator to read the file (tab).
#' @param skip number of rows to skip because comments (1).
#'
#' @return a df of IS annotated with sharing counts (not sums).
#' @usage TODO
#' @note : 
#'
###############################################################
sharingSourceToOutput_unsorted <- function(df,
                                           source_regex = "CD34",
                                           output_regex = "CD13|CD14|CD15|CD3|CD4|CD8|CD19|CD36|CD38|GLY|Plasma",
                                           annotation_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                                           # min_sharing_source_output = 2,
                                           return_allinputcols = TRUE) {
  message(paste("[AP]\tComputing sharing from a source (named with regex) and output/offspring (as regex)."))
  source_cols <- grep(source_regex, colnames(df), value = T)
  output_cols <- grep(output_regex, colnames(df), value = T)
  # df_paired <- df[c(annotation_cols, source_cols, output_cols)]
  # df_source <- df[c(annotation_cols, source_cols)]
  # df_output <- df[c(annotation_cols, output_cols)]
  df_paired <- df[c(source_cols, output_cols)]
  df_source <- df[c(source_cols)]
  df_output <- df[c(output_cols)]
  df_source[is.na(df_source)] <- 0
  df_output[is.na(df_output)] <- 0
  df_paired[is.na(df_paired)] <- 0
  df_source$Source_Sharing <- apply(df_source[source_cols], 1, function(x) {length(x[x>0])} )
  df_output$Output_Sharing <- apply(df_output[output_cols], 1, function(x) {length(x[x>0])} )
  df_paired$Paired_Sharing <- apply(cbind(df_source[c("Source_Sharing")], df_output[c("Output_Sharing")]), 1, function(x) {length(x[x>0])} )
  df_out <- NULL
  if (return_allinputcols) {
    df_out <- cbind(df, df_source[c("Source_Sharing")], df_output[c("Output_Sharing")], df_paired[c("Paired_Sharing")]) 
  } else {
    df_out <- cbind(df[annotation_cols], df_source, df_output, df_paired[c("Paired_Sharing")]) 
  }
  return (df_out)
}



###############################################################
#' @title Compute pairwise Fisher test on gene frequencies from CIS results
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-07-08
#'
#' @rdname compareGeneFrequency_Fisher
#' @docType methods
#' @aliases compareGeneFrequency_Fisher
#'
#' @param g1_cis_df CIS results df group 1.
#' @param g2_cis_df CIS results df group 2.
#' @param min_is_per_gene Min number of IS per gene to be reported as CIS gene.
#'
#' @return a df of gene frequences 
#' @usage TODO
#' @description This function is aimed at comparing the groups (G1 and G2) in terms of gene frequencies
#' by applying Fisher exact test on the confusion matrix 2x2 by gene (FDR corrected).
#' @note : 
#'
###############################################################
compareGeneFrequency_Fisher <- function(g1_cis_df, g2_cis_df, 
                                        min_is_per_gene = 3, 
                                        gene_set_method = "INTERSECTION",
                                        scary_genes_toannotate = c("Lmo2", "Smg6", "Mecom", "Mds", "Ccnd2"),
                                        remove_only_unbalanced_0 = TRUE ){
  message(paste("[AP]\tComparing G1 and G2 gene frequencies with Fisher test."))
  # filename_infix <- ".Intersection" # ".Intersection" ".Union"
  # gene_set_method <- "INTERSECTION" # "INTERSECTION" "UNION"
  # min_is_per_gene <- 3
  
  significance_threshold_minus_log_p <- -log(0.05, base = 10)
  annotation_threshold_ontots <- 0.5
  annotation_threshold_ontots_log <- -log(annotation_threshold_ontots, base = 10)
  annotation_cols_to_get <- c("GeneName", "OncoGene", "TumorSuppressor", "Onco1_TS2", "ClinicalRelevance", "DOIReference", "KnownGeneClass", "KnownClonalExpension", "CriticalForInsMut")
  annotation_cols_to_get_asSQLstring <- paste(annotation_cols_to_get, collapse = ", ")
  
  # dcast(data = gdf_iss_abundance_paperslice_melt, chr+ integration_locus + strand + GeneName + GeneStrand ~ SampleID, value.var = "Abundance", fun.aggregate = sum)
  # 
  # gdf_iss_abundance_paperslice_melt_cast <- dcast(data = gdf_iss_abundance_paperslice_melt, chr+ integration_locus + strand + GeneName + GeneStrand ~ SampleID, value.var = "Abundance", fun.aggregate = sum)
  # 
  # formula for RPKM: RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
  # add info of gene len normalization
  g1_cis_df$IS_per_kbGeneLen <- (g1_cis_df$n_IS_perGene / g1_cis_df$average_TxLen)*1000
  g1_cis_df$Sum_IS_per_kbGeneLen <- sum(g1_cis_df$IS_per_kbGeneLen)
  g1_cis_df$IS_per_kbGeneLen_perMDepth_TPM <- (g1_cis_df$IS_per_kbGeneLen / g1_cis_df$Sum_IS_per_kbGeneLen)*1000000
  g1_cis_df$KnownGeneClass <- ifelse(!is.na(g1_cis_df$Onco1_TS2), (ifelse((g1_cis_df$Onco1_TS2) == 1, "OncoGene", "TumSuppressor")), "Other")
  g1_cis_df$CriticalForInsMut <- ifelse(!is.na(g1_cis_df$KnownClonalExpension), "True", "False")
  g1_cis_df[is.na(g1_cis_df)] <- NA
  # write.table(x = g1_cis_df, file = paste(cis_results_folder, ProjectID, ".CIS.results.G1.extended.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = F, col.names = T)
  
  g2_cis_df$IS_per_kbGeneLen <- (g2_cis_df$n_IS_perGene / g2_cis_df$average_TxLen)*1000
  g2_cis_df$Sum_IS_per_kbGeneLen <- sum(g2_cis_df$IS_per_kbGeneLen)
  g2_cis_df$IS_per_kbGeneLen_perMDepth_TPM <- (g2_cis_df$IS_per_kbGeneLen / g2_cis_df$Sum_IS_per_kbGeneLen)*1000000
  g2_cis_df$KnownGeneClass <- ifelse(!is.na(g2_cis_df$Onco1_TS2), (ifelse((g2_cis_df$Onco1_TS2) == 1, "OncoGene", "TumSuppressor")), "Other")
  g2_cis_df$CriticalForInsMut <- ifelse(!is.na(g2_cis_df$KnownClonalExpension), "True", "False")
  g2_cis_df[is.na(g2_cis_df)] <- NA
  # write.table(x = g2_cis_df, file = paste(cis_results_folder, ProjectID, ".CIS.results.G1.extended.tsv", sep = ""), sep = "\t", quote = FALSE, row.names = F, col.names = T)
  
  g1_cis_df_forfreqplot <- g1_cis_df[which(g1_cis_df$n_IS_perGene >= min_is_per_gene),]
  g2_cis_df_forfreqplot <- g2_cis_df[which(g2_cis_df$n_IS_perGene >= min_is_per_gene),]
  
  study_mergegenename <- NULL
  if (gene_set_method == "UNION") {
    ### -------- version 1: Group 1 OR Group 2 genes ----------- equivalent to UNION ###
    union_of_genes <- union(as.character(g1_cis_df_forfreqplot[,c("GeneName")]), as.character(g2_cis_df_forfreqplot[,c("GeneName")]))
    new_overall_genet_cast_merge <- merge(x = g1_cis_df[which(g1_cis_df$GeneName %in% union_of_genes),], y = g2_cis_df[which(g2_cis_df$GeneName %in% union_of_genes),], by = annotation_cols_to_get)
    names(new_overall_genet_cast_merge) <- gsub("\\.x", "_G1", colnames(new_overall_genet_cast_merge))
    names(new_overall_genet_cast_merge) <- gsub("\\.y", "_G2", colnames(new_overall_genet_cast_merge))
    study_mergegenename <- new_overall_genet_cast_merge  
  } else if (gene_set_method == "INTERSECTION") {
    ### -------- version 2: Group 1 AND Group 2 genes ----------- equivalent to INTERSECTION ###
    # gene intersection
    overall_genet_cast_merge <- merge(x = g1_cis_df_forfreqplot, y = g2_cis_df_forfreqplot, by = annotation_cols_to_get)
    names(overall_genet_cast_merge) <- gsub("\\.x", "_G1", colnames(overall_genet_cast_merge))
    names(overall_genet_cast_merge) <- gsub("\\.y", "_G2", colnames(overall_genet_cast_merge))
    # use a backup
    study_mergegenename <- overall_genet_cast_merge
  } else {
    message(paste0("[AP]\tERROR: wrong method selected for gene_set_method."))
  }
  
  # check len and return results
  if (nrow(study_mergegenename) == 0) {
    message(paste0("[AP]\tWARNING: No genes in comomn resulted in G1 and G2 for frequency comparison. Rerutning NULL."))
    return (NULL)
  } else {
    ### --- now move ahead wuth the rest --- 
    # now do comparison by gene frequency
    study_mergegenename$scary <- ifelse(study_mergegenename$GeneName %in% scary_genes_toannotate, TRUE, FALSE)
    selected_col_for_fisher <- c("n_IS_perGene_G1", "TotIS_asDfRow_G1", "n_IS_perGene_G2", "TotIS_asDfRow_G2")
    # selected_col_for_fisher <- c("rd027_NISperGene", "rd027_TotalIS", "tm035_NISperGene", "tm035_TotalIS")
    ft_pval <- data.frame("FisherTest_pvalue" = apply(study_mergegenename[selected_col_for_fisher], 1, function(x) {
      .m <- matrix(c(x[1], x[2]-x[1], x[3], x[4]-x[3]), nrow = 2, dimnames = list(RD027 = c("ISofGene", "TotalIS"), TM035 = c("ISofGene", "TotalIS")))
      .ft <- fisher.test(.m)
      .ft$p.value
    })
    )
    study_mergegenename <- cbind(study_mergegenename, ft_pval)
    study_mergegenename$FisherTest_pvalue_significant <- ifelse(study_mergegenename$FisherTest_pvalue < 0.05, TRUE, FALSE)
    
    ### ---- to correct the pvalue, first try including or not the 0
    # get the min number of IS for 0 genes in the other sample
    mean_nis_pergene_a <- ceiling(mean(study_mergegenename$n_IS_perGene_G1[study_mergegenename$n_IS_perGene_G1>0]))
    mean_nis_pergene_c <- ceiling(mean(study_mergegenename$n_IS_perGene_G2[study_mergegenename$n_IS_perGene_G2>0]))
    study_mergegenename$to_exclude_from_test <- apply(study_mergegenename[selected_col_for_fisher], 1, function(x) {
      ifelse( (x[1] == 0 | x[3] == 0),
              ifelse( (((x[1] < mean_nis_pergene_a) & (x[3] == 0)) | ((x[1] == 0) & (x[3] < mean_nis_pergene_c))),
                      TRUE,
                      FALSE),
              FALSE
      )
    })
    
    if (remove_only_unbalanced_0) {
      study_mergegenename <- study_mergegenename[which(study_mergegenename$to_exclude_from_test == FALSE),] 
    }
    
    # corrected pvalue
    # study_mergegenename$FisherTest_pvalue_bonferroni <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "bonferroni", n = length(study_mergegenename$FisherTest_pvalue))
    study_mergegenename$FisherTest_pvalue_fdr <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "fdr", n = length(study_mergegenename$FisherTest_pvalue))
    study_mergegenename$FisherTest_pvalue_benjamini <- p.adjust(study_mergegenename$FisherTest_pvalue, method = "BY", n = length(study_mergegenename$FisherTest_pvalue))
    # add vars for volcano plot
    study_mergegenename$log_FC <- log(x = (study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2), base = 10)
    study_mergegenename$FC_G1G2 <- study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2
    study_mergegenename$log_pvalue <- -log(x = study_mergegenename$FisherTest_pvalue, base = 10)
    study_mergegenename$log_pvalue_fdr <- -log(x = study_mergegenename$FisherTest_pvalue_fdr, base = 10)
    # study_mergegenename$log_pvalue <- -log(x = study_mergegenename$FisherTest_pvalue, base = 10)
    study_mergegenename$log2_FC <- log(x = (study_mergegenename$geneIS_frequency_byHitIS_G1 / study_mergegenename$geneIS_frequency_byHitIS_G2), base = 2)
    study_mergegenename$log2_FC_TPM <- log(x = (study_mergegenename$IS_per_kbGeneLen_perMDepth_TPM_G1 / study_mergegenename$IS_per_kbGeneLen_perMDepth_TPM_G2), base = 2)
    
    return(study_mergegenename)
  } # if (nrow(study_mergegenename) == 0)
  
}

###############################################################
#' @title Parse read CIGAR and get info of chimera
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-10-05
#'
#' @rdname parseCigarForChimera
#' @docType methods
#' @aliases parseCigarForChimera
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return a df of sureads results 
#' @usage TODO
#' @description parse CIGAR string and get data for AAV studies (and general chimera). Definition of chimera (as subread): the ONLY portion of the read (identified by CIGAR) that flanks target genome alignment.
#' @note : 
#'
###############################################################
parseCigarForChimera <- function(df_alm, 
                                 col_read_name = "name",
                                 cols_to_include = c("start", "end", "gene_chr", "gene_annotationsource", "gene_elementtype", "gene_start", "gene_end", "gene_strand", "gene_details", "distance_to_gene", "sample", "sourcefile", "score"),
                                 key_cols = c("name", "cigar"),
                                 cigar_opts_to_keep = c('M', 'D', 'I'),
                                 vector_chr_string = "chrV",
                                 min_cigar_alm_width = 5, 
                                 quiet = T
                                 ){
  require(sqldf)
  require(GenomicAlignments)
  
  message(paste("[AP]\tParse CIGAR string to find query alignments (re-arrangements) and Chimera"))
  reads_id <- levels(factor(df_alm[, col_read_name]))
  df_query_summary <- NULL
  cigar_opts_to_keep <- paste0("('", paste(cigar_opts_to_keep, collapse = "', '"), "')")
  for (read_id in reads_id) {
    if (!(quiet)) {
      message(paste0("[AP]\t\t", read_id))
    }
    df_ir <- NULL
    slice_t <- df_alm[which(df_alm$name == read_id),]
    for (cigar_string in as.character(slice_t[,c("cigar")])) {
      ir_t <- cigarRangesAlongQuerySpace(cigar_string, 
                                         flag=NULL, before.hard.clipping=T, after.soft.clipping=F, 
                                         ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=FALSE, with.ops=T)
      ir_t_len <- as.numeric(width(cigarRangesAlongQuerySpace(cigar_string, flag=NULL, before.hard.clipping=T, after.soft.clipping=F, ops=CIGAR_OPS, drop.empty.ranges=FALSE, reduce.ranges=T, with.ops=T)))
      ir_t_chr <- slice_t[which(slice_t$cigar == cigar_string), "chr"]
      ir_t_strand <- slice_t[which(slice_t$cigar == cigar_string), "strand"]
      cigar_df <- as.data.frame(ir_t[[1]])
      names(cigar_df) <- c("read_alm_start", "read_alm_end", "read_alm_width", "cigar_flag")
      read_info <- data.frame( "name" = rep(read_id, nrow(cigar_df)),
                               "cigar" = rep(cigar_string, nrow(cigar_df)),
                               "chr" = ir_t_chr,
                               "strand" = ir_t_strand,
                               "read_len" = ir_t_len
      )
      bind_df <- cbind(read_info, cigar_df)
      # if read is reverse, invert iranges 
      if (levels(factor(bind_df$strand)) == '-') {
        bind_df$query_start <- apply( bind_df[c("read_len", "read_alm_start", "read_alm_end", "read_alm_width")], 1, function(x) {ifelse( x[4] > 0, 
                                                                                                                                          ir_t_len - x[3] + 1, 
                                                                                                                                          ir_t_len - x[2] + 1
        )
        }
        )
        bind_df$query_end <- ir_t_len - bind_df$read_alm_start + 1
      } else {
        bind_df$query_start <- bind_df$read_alm_start
        bind_df$query_end <- bind_df$read_alm_end
      }
      bind_df <- merge(x = bind_df, y = slice_t[which(slice_t$cigar == cigar_string), c(key_cols, cols_to_include)], by = key_cols, all.x = T)
      
      # merge df with all other cigar strings
      if (length(df_ir) == 0) {
        df_ir <- bind_df
      } else {
        df_ir <- rbind(df_ir, bind_df)
      } # if (length(df_ir) == 0)
    } # for (cigar_string in as.character(slice_t[,c("cigar")])) 
    
    # summary 
    read_query_summary <- sqldf( paste0("select * from df_ir where cigar_flag in ", cigar_opts_to_keep, " and read_alm_width >= ", min_cigar_alm_width, " order by query_start") )
    
    # create a string of read composition
    read_query_string_alm <- sqldf( paste0("select cigar, chr from df_ir where cigar_flag in ", cigar_opts_to_keep, " and read_alm_width >= ", min_cigar_alm_width, " group by cigar order by query_start") )
    read_query_string_alm_type <- gsub("chr", "", paste0(as.character(read_query_string_alm[,"chr"]), collapse = '-'))
    # do we have multi-mapping on target genome?
    read_query_string_alm_multimapping <- ifelse( (length(levels(factor( (read_query_string_alm[which( !(read_query_string_alm$chr %in% c(vector_chr_string)) ), "chr"]) ))) > 1), TRUE, FALSE)
    
    # ---- find chimera(s) ----
    # init chimera
    read_query_summary$chimera <- FALSE
    read_query_summary$junction <- FALSE
    read_query_summary$target_genome_position_wrt_junction <- NA
    read_query_summary$read_alm_type <- read_query_string_alm_type
    read_query_summary$target_genome_multimapping <- read_query_string_alm_multimapping
    
    vector_subread_chimera_cigar <- NULL
    alignment_on_target_genome <- F
    target_genome_position_wrt_junction <- NA
    target_genome_position_wrt_junction_found <- F
    if (nrow(read_query_summary) > 1) {
      for (i in seq(1, (nrow(read_query_summary)-1))) {
        # case A: V+G+
        # case B: G+V+
        # case C: G+V+G+ -> this case is included in A and B
        
        # case A: V+G+
        if (read_query_summary[i,"chr"] == vector_chr_string & read_query_summary[(i+1),"chr"] != vector_chr_string) {
          vector_subread_chimera_cigar <- c(vector_subread_chimera_cigar, as.character(read_query_summary[i,"cigar"]))
          if ( !target_genome_position_wrt_junction_found ) {
            target_genome_position_wrt_junction <- "R"
            target_genome_position_wrt_junction_found <- T
          } else {
            target_genome_position_wrt_junction <- paste0(target_genome_position_wrt_junction, "R")
          }
        } 
        # case B: G+V+
        if (read_query_summary[i,"chr"] != vector_chr_string & read_query_summary[(i+1),"chr"] == vector_chr_string) {
          vector_subread_chimera_cigar <- c(vector_subread_chimera_cigar, as.character(read_query_summary[(i+1),"cigar"]))
          if ( !target_genome_position_wrt_junction_found ) {
            target_genome_position_wrt_junction <- "L"
            target_genome_position_wrt_junction_found <- T
          } else {
            target_genome_position_wrt_junction <- paste0(target_genome_position_wrt_junction, "L")
          }
        } 
        # end if
      } # for (i in seq(1, (nrow(read_query_summary)-1)))   
    } else {
      vector_subread_chimera_cigar <- c()
    } # if (nrow(read_query_summary) > 1)
    read_query_summary$target_genome_position_wrt_junction <- target_genome_position_wrt_junction
    read_query_summary$junction <- ifelse(read_query_summary$cigar %in% vector_subread_chimera_cigar, TRUE, FALSE)
    # read_query_summary$chimera <- ifelse("TRUE" %in% as.character(read_query_summary$junction), TRUE, FALSE)
    read_query_summary$chimera <- ifelse(length(vector_subread_chimera_cigar) > 0, TRUE, FALSE)
    # now find the aav junction side, to best characterize the aav juction - integration. complex cases are here NOT managed (TODO)
    read_query_summary$vector_junction_locus <- ifelse(read_query_summary$junction == T, 
                                                       ifelse(read_query_summary$strand =='+',
                                                              ifelse(read_query_summary$target_genome_position_wrt_junction %in% c('L', 'LR'), 
                                                                     read_query_summary$start,
                                                                     read_query_summary$end), 
                                                              ifelse(read_query_summary$target_genome_position_wrt_junction %in% c('L', 'LR'), 
                                                                     read_query_summary$end,
                                                                     read_query_summary$start)),
                                                       NA)
    read_query_summary$query_order <- seq(1, nrow(read_query_summary), 1)
    # find locus
    locus_found_and_reported <- F
    indels <- NULL
    
    # read_query_summary$integration_locus <- apply(read_query_summary[c("chr", "strand", "target_genome_position_wrt_junction", "start", "end", "chimera")], 1, function(x) {
    #         .locus <- NA    
    #         if (x[1] != vector_chr_string & x[6] == "TRUE" & !is.na(x[3])) { 
    #           # simple cases L, R with - +
    #           if (x[2] == "+" & x[3] == "L") { .locus <- x[5]}
    #           if (x[2] == "+" & x[3] == "R") { .locus <- x[4]}
    #           if (x[2] == "-" & x[3] == "L") { .locus <- x[4]}
    #           if (x[2] == "-" & x[3] == "R") { .locus <- x[5]}
    #           # complex cases LR with +-
    #           if (!(locus_found_and_reported) & x[3] == "LR" & x[2] == "+") { 
    #             .locus <- x[5]
    #             locus_found_and_reported <- T
    #           } # LR+
    #           if (!(locus_found_and_reported) & x[3] == "LR" & x[2] == "-") { 
    #             .locus <- x[4]
    #             locus_found_and_reported <- T
    #           } # LR-
    #           
    #           # # just look at indels at AAV site with complex cases
    #           # if (locus_found_and_reported & x[3] == "LR" & x[2] == "+") { 
    #           #   .locus <- x[5]
    #           #   locus_found_and_reported <- T
    #           # } # LR+
    #           
    #         } # if (x[1] != vector_chr_string) 
    #         .locus
    #       } # function
    #     ) # apply
    
    integration_locus_tobind <- NA
    # slice_data <- sqldf(paste0('select distinct chr, strand, target_genome_position_wrt_junction, start, end, chimera from read_query_summary where chr not like "', vector_chr_string, '"'))
    slice_data <- sqldf(paste0('select distinct chr, strand, target_genome_position_wrt_junction, start, end, chimera from read_query_summary where 1'))
    locus_found_and_reported <- FALSE
    for (elem in seq(1, nrow(slice_data)) ) {
      .locus <- integration_locus_tobind
      # message(paste(slice_data[elem, "name"], slice_data[elem, "start"], slice_data[elem, "end"], slice_data[elem, "target_genome_position_wrt_junction"], "\n"))
      # message(paste("--> e", elem, as.character(locus_found_and_reported), "\n"))
      if ( slice_data[elem, "chr"] != vector_chr_string & slice_data[elem, "chimera"] == "TRUE" & !is.na(slice_data[elem, "target_genome_position_wrt_junction"]) ) { 
        # simple cases L, R with - +
        if (slice_data[elem, "strand"] == "+" & slice_data[elem, "target_genome_position_wrt_junction"] == "L") { .locus <- slice_data[elem, "end"]}
        if (slice_data[elem, "strand"] == "+" & slice_data[elem, "target_genome_position_wrt_junction"] == "R") { .locus <- slice_data[elem, "start"]}
        if (slice_data[elem, "strand"] == "-" & slice_data[elem, "target_genome_position_wrt_junction"] == "L") { .locus <- slice_data[elem, "start"] }
        if (slice_data[elem, "strand"] == "-" & slice_data[elem, "target_genome_position_wrt_junction"] == "R") { .locus <- slice_data[elem, "end"] }
        # complex cases LR with +-
        if (!(locus_found_and_reported) & slice_data[elem, "target_genome_position_wrt_junction"] == "LR" & slice_data[elem, "strand"] == "+") { 
          .locus <- slice_data[elem, "end"]
          locus_found_and_reported <- T
          # message(paste("\tFOUND:", as.numeric(.locus), as.character(locus_found_and_reported), "\n"))
        } # LR+
        if (!(locus_found_and_reported) & slice_data[elem, "target_genome_position_wrt_junction"] == "LR" & slice_data[elem, "strand"] == "-") { 
          .locus <- slice_data[elem, "start"]
          locus_found_and_reported <- T
          # message(paste("\tFOUND:", as.numeric(.locus), as.character(locus_found_and_reported), "\n"))
        } # LR-
        
        # # just look at indels at AAV site with complex cases
        # if (locus_found_and_reported & x[3] == "LR" & x[2] == "+") { 
        #   .locus <- x[5]
        #   locus_found_and_reported <- T
        # } # LR+
        integration_locus_tobind <- as.numeric(.locus)    
      } # if (x[1] != vector_chr_string) 
    } # for (elem in nrow(slice_data))
    
    read_query_summary$integration_locus <- ifelse(!(read_query_summary$chr %in% vector_chr_string), integration_locus_tobind, NA)
    
    # merge df with all other cigar strings
    if (length(df_query_summary) == 0) {
      df_query_summary <- read_query_summary
    } else {
      df_query_summary <- rbind(df_query_summary, read_query_summary)
    } # if (length(df_ir) == 0)
  } # for (header in reads_id) 
  
  return (df_query_summary)
}


plotRanges <- function(x, xlim=x, main=deparse(substitute(x)), col="black", sep=0.5, ...) {
  height <- 1
  # if (is(xlim, "IntegerRanges"))
  xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col=col, ...)
  title(main)
  axis(1)
}

extend_chromosomes = function(bed, chromosome, prefix = "zoom_") {
  zoom_bed = bed[bed[[1]] %in% chromosome, , drop = FALSE]
  zoom_bed[[1]] = paste0(prefix, zoom_bed[[1]])
  rbind(bed, zoom_bed)
}



###############################################################
#' @title Plot chimera in CIRCOS
#' 
#' @author Andrea Calabria
#' @details version 0.1, 2020-11-19
#'
#' @rdname plot_circos_rearrangement_imp
#' @docType methods
#' @aliases plot_circos_rearrangement_imp
#'
#' @param df_alm input df with alignment results
#' @param min_cigar_alm_width minimum alignment subread widt for each CIGAR string. default = 5. This value will be applied to ALL tags (indels included)
#'
#' @return plot
#' @description parse CIGAR string and get data for AAV studies (and general chimera). Definition of chimera (as subread): the ONLY portion of the read (identified by CIGAR) that flanks target genome alignment.
#' @note : 
#' @usage 
#' plot_circos_rearrangement_imp(df = alm_reads_for_rearrangements_withstatsbyread, sample_read = "m64047_200620_005306/19269857/ccs", 
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.circos.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), vector_chr = "chrV", zoomed_chr_index = 23, vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband")
###############################################################
plot_circos_rearrangement_imp <- function(df, sample_read, outfile_png, outfile_pdf,
                                          cols_to_search = c("chr", "start", "end", "ratio_bp_aligned_on_raw"), 
                                          vector_chr = "chrV", 
                                          zoomed_chr_index = 23, 
                                          vector_cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
                                          species = "hg19",
                                          bp_res= 300,
                                          color_vector_fwd="orange",
                                          color_vector_rev="green",
                                          color_target_rev="violet",
                                          color_target_fwd="blue") {
  if (FALSE %in% (cols_to_search %in% colnames(df)) ) {
    message(paste0("[AP]\tERROR: Columns in df are different from expected: ", paste(cols_to_search, collapse = ',')))
    return (NULL)
  } else {
    # slice data
    t <- sample_read
    slice_t <- df[which(df$name == t),]
    slice_t <- slice_t[order(slice_t$query_start),]
    bed_1 <- slice_t[(1:(nrow(slice_t)-1)), cols_to_search]
    names(bed_1) <- c("chr", "start", "end", "value1")
    bed_2 <- slice_t[(2:(nrow(slice_t))), cols_to_search]
    names(bed_2) <- c("chr", "start", "end", "value1")
    
    # orientations
    bed_all_stranded <- slice_t[c("chr", "start", "end", "strand")]
    bed_all_stranded$strand <- ifelse(bed_all_stranded$strand == '+', 1, 0)
    names(bed_all_stranded) <- c("chr", "start", "end", "value1")
    
    # get start-end of each read
    bed_all_long <- slice_t[c("chr", "start", "end", "strand", "query_start", "query_end")]
    outbed_line_sorted_start <- NULL
    if (nrow(bed_all_long) > 1) {
      for (i in seq(1, (nrow(bed_all_long)-1))) {
        if (bed_all_long[i,"strand"] == '+') {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("end")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } else {
          from_chr <- bed_all_long[i,c("chr")]
          from_pos <- bed_all_long[i,c("start")]
          if (bed_all_long[i+1,"strand"] == '+') {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("start")]
          } else {
            to_chr <- bed_all_long[i+1,c("chr")]
            to_pos <- bed_all_long[i+1,c("end")]
          } # if (bed_all_long[i+1,"strand"] == '+')
        } # if (bed_all_long[i,"strand"] == '+')
        if (length(outbed_line_sorted_start) > 0) {
          outbed_line_sorted_start <- rbind(outbed_line_sorted_start, data.frame("from_chr" = from_chr,
                                                                                 "from_pos" = from_pos,
                                                                                 "to_chr" = to_chr,
                                                                                 "to_pos" = to_pos))  
        } else {
          outbed_line_sorted_start <- data.frame("from_chr" = from_chr,
                                                 "from_pos" = from_pos,
                                                 "to_chr" = to_chr,
                                                 "to_pos" = to_pos)
        } # if (length(outbed_line_sorted_start) > 0)
        
      } # for (i in seq(1, nrow(bed_all_long)))
    } # if (nrow(bed_all_long) > 1)
    outbed_line_sorted_start$value1 <- 1
    bed_1_line <- outbed_line_sorted_start[grep("from", colnames(outbed_line_sorted_start))]
    names(bed_1_line) <- c("chr", "start")
    bed_1_line$end <- bed_1_line$start
    bed_1_line$value1 <- 1
    bed_2_line <- outbed_line_sorted_start[grep("to", colnames(outbed_line_sorted_start))]
    names(bed_2_line) <- c("chr", "start")
    bed_2_line$end <- bed_2_line$start
    bed_2_line$value1 <- 1
    
    target_chr <- gsub("chr", "", setdiff( levels(factor(slice_t$chr)), vector_chr))
    # if 0 -> only AAV (chrV), else target too
    if (length(target_chr) > 0) {
      # go into genomics
      human_cytoband <- read.cytoband(species = species)$df
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(cytoband_rbind)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # acquire target chr
      if (target_chr %in% c("X", "Y")) {
        if (target_chr == "X") {target_chr <- grep("X", chromosome)}
        if (target_chr == "Y") {target_chr <- grep("Y", chromosome)}
      } else {
        target_chr <- as.numeric(target_chr)
      }
      # determine range of zoom
      xrange = c(cytoband$chr.len, cytoband$chr.len[vector_chr])
      normal_chr_index = target_chr
      zoomed_chr_index = zoomed_chr_index
      sector.width = c(xrange[normal_chr_index] / sum(xrange[normal_chr_index]), 
                       xrange[zoomed_chr_index] / sum(xrange[zoomed_chr_index])) 
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 300)
      extended <-extend_chromosomes(cytoband_df, vector_chr)
      chr_target_ff<- bed_all_stranded[bed_all_stranded$chr!=vector_chr,]$chr
      s <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$start
      e <- bed_all_stranded[bed_all_stranded$chr==chr_target_ff,]$end
      m <- (s+e)/2
      r <- extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,]
      increase<-slice_t$read_len[1]+1000
      r$V2 <- m-(increase/2)
      r$V3 <- m+(increase/2)
      extended[(extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m) ,] <- r
      itr5 <- extended[extended$V4=="ITR.5p" & extended$V1==vector_chr,]$V2
      itr3 <- extended[extended$V4=="ITR.3p" & extended$V1==vector_chr,]$V3
      extended2 <- extended[((extended$V1==chr_target_ff) & (extended$V2<m) & (extended$V3>m)) | (extended$V1!=chr_target_ff & extended$V1!=vector_chr)
                            | ((extended$V1==vector_chr & (extended$V2>=itr5 & extended$V3<=itr3))),]
      extended2$V6 <-''
      extended2[extended2$V1==vector_chr,]$V6 <- extended2[extended2$V1==vector_chr,]$V4
        
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
    
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      l <- seq(m-5000,m+5000,600)
      for (subread in seq(1, nrow(bed_all_stranded))) {
          circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                              track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                              panel.fun = function(region, value, ...) {
                                                  x1=region$start
                                                  y1=region$end
                                                  if(value==1){
                                                    v<-value
                                                    if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                      v<-v+2
                                                      circos.arrow(x1,y1, col = f(v), 
                                                                   border = 1, arrow.head.length = cm_x(0.2) ) 
                                                  }
                                                  else{
                                                    v<-value
                                                    if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                      v<-v+2
                                                      circos.arrow(x1, y1, col = f(v), 
                                                                   border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                                  }
                                                 
                                                  
                                                  i = getI(...)
                                                  cell.xlim = get.cell.meta.data("cell.xlim")
                                                  # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                              })
      

      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.initializeWithIdeogram(extended2, 
                                    chromosome.index = c(chromosome[target_chr], vector_chr),
                                    sector.width = sector.width, tickLabelsStartFromZero = FALSE, major.by = bp_res)
    
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                              track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                              panel.fun = function(region, value, ...) {
                                                x1=region$start
                                                y1=region$end
                                                if(value==1){
                                                  v<-value
                                                  if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                    v<-v+2
                                                  circos.arrow(x1,y1, col = f(v), 
                                                               border = 1, arrow.head.length = cm_x(0.2) ) 
                                                }
                                                else{
                                                  v<-value
                                                  if (bed_all_stranded[subread,]$chr!=vector_chr)
                                                    v<-v+2
                                                  circos.arrow(x1, y1, col = f(v), 
                                                               border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                                }
                                                
                                                
                                                i = getI(...)
                                                cell.xlim = get.cell.meta.data("cell.xlim")
                                                # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                              })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      
    } else {
      vector_cytoband <- read.csv(file = vector_cytoband_file, header=F, fill=T, check.names = FALSE, sep = '\t')
      # cytoband_rbind <- rbind(human_cytoband, vector_cytoband)
      cytoband <- read.cytoband(vector_cytoband)
      
      cytoband_df = cytoband$df
      chromosome = cytoband$chromosome
      # PNG
      png(file = outfile_png, height=5, width=5, units = "in", res = 300)
      circos.par("start.degree" = 90)
      
      itr5 <- cytoband_df[cytoband_df$V4=="ITR.5p" & cytoband_df$V1==vector_chr,]$V2
      
      itr3 <- cytoband_df[cytoband_df$V4=="ITR.3p" & cytoband_df$V1==vector_chr,]$V3 
      
      cytoband_df <- cytoband_df[(cytoband_df$V2>=itr5 & cytoband_df$V3<=itr3),]
      
      circos.initializeWithIdeogram(cytoband = cytoband_df,
                                    chromosome.index = c(vector_chr), 
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
      # PDF
      pdf(file = outfile_pdf, height=5, width=5)
      circos.par("start.degree" = 90)
      circos.initializeWithIdeogram(cytoband = cytoband_df, chromosome.index = c(vector_chr),
                                    tickLabelsStartFromZero = FALSE, major.by = bp_res)
      f = colorRamp2(breaks = c(0, 1,2,3), colors = c(color_vector_rev, color_vector_fwd,color_target_rev,color_target_fwd))
      for (subread in seq(1, nrow(bed_all_stranded))) {
        circos.genomicTrackPlotRegion(bed_all_stranded[subread,], stack = TRUE, 
                                      track.height = 0.05, bg.col = NA, bg.border = "gray90",
                                      panel.fun = function(region, value, ...) {
                                        x1=region$start
                                        y1=region$end
                                        if(value==1){
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1,y1, col = f(v), 
                                                       border = 1, arrow.head.length = cm_x(0.2) ) 
                                        }
                                        else{
                                          v<-value
                                          if (bed_all_stranded[subread,]$chr!=vector_chr)
                                            v<-v+2
                                          circos.arrow(x1, y1, col = f(v), 
                                                       border = 1, arrow.position = "start", arrow.head.length = cm_x(0.2))
                                        }
                                        
                                        
                                        i = getI(...)
                                        cell.xlim = get.cell.meta.data("cell.xlim")
                                        # circos.lines(cell.xlim, c(i, i), lty = 2, col = "#FFFFFF")
                                      })
      }
      # circos.genomicLink(bed_1, bed_2, col = rand_color(nrow(bed_1), transparency = 0.5))
      circos.genomicLink(bed_1_line, bed_2_line, col = "black", directional = 1)
      dev.off()
      circos.clear()
    }  # if (length(target_chr) > 0) 
  } # if (FALSE %in% (cols_to_search %in% colnames(df)) )
  
} # function



###############################################################
#' @title Plot rearrangements linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_rearrangement
#'
#' @param df_alm input df with alignment results
#' @param reads_to_plot vector with names of the reads to plot
#' @param cytoband_file path to cytoband file NOTE: the file must have the following
#'                      header: chr	start	end	name	gieStain 
#' @param vector_length vector genome length
#' @return 
#' @description Plot the vector portions of the read
#' @note : 
#' @usage 
#' plot_linear_rearrangement(df = alm_reads_for_rearrangements_withstatsbyread, reads_to_plot = reads_to_plot<-c("m64047_200324_183710/83362987/ccs" ,"m64047_200324_183710/102237007/ccs"),
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' vector_chr = "chrV", cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
#' annotation_file="source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.bed")
plot_linear_rearrangement <- function(df_alm, reads_to_plot, outfile_png, outfile_pdf,
                                      cytoband_file, annotation_file, vector_chr = "chrV",
                                      vector_length=6000,
                                      height=8, width=8, res = 300,
                                      offset_arrow=110
                                      ){
  library(karyoploteR)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 15
  plot.params$data1height <- 10
  plot.params$data2height<-500
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  increase1<-0.2/length(reads_to_plot)
  increase2<-0.22/length(reads_to_plot)
  increase3<-0.3/length(reads_to_plot)
  
  i_r0<-0.1
  i_r1<-0.15
  y0<-0.1
  y1<-y0+increase1
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                     'target_genome_position_wrt_junction')]
    
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    for (i in 1:nrow(x))
    {
      if (x[i,]$strand=='+')
      {       
        i_col<-"orange"
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-"green"
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.05,lwd=15,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.05,lwd=13,lend='butt',ljoin='mitre', data.panel = 2)
      print(i_r0)
      print(i_r1)
      
      i_r0<-i_r0+increase1
      i_r1<-i_r1+increase1
    }
    
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    if (!is.na(x[1,]$vector_junction_locus))
    {
      j<-x[1,]$vector_junction_locus
     # if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
      #  j<-j+offset_arrow
      #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
        #j<-j-offset_arrow
      kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
    }
    if (nrow(x) > 1)
    {
      for (i in 2:nrow(x))
      {
        
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y0+0.0025*length(reads_to_plot) ,y1=y1,length=0.1, lty=2, data.panel = 2)
        y0<-y0+increase2
        y0_line<-y0_line+increase2
        y1<-y1+increase2
        if (!is.na(x[i,]$vector_junction_locus))
        {
          j<-x[i,]$vector_junction_locus
          #if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
           # j<-j+offset_arrow
          #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
           # j<-j-offset_arrow
          print(x[i,]$vector_junction_locus)
          kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
        }
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end-(offset_arrow/2) else x[i,]$start+(offset_arrow/2)
        
      }
    }
    y0<-y0+increase1
    y1<-y1+increase1
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(y0,y0),lty=2,data.panel = 2,col='grey',lwd=2)
    i_r0<-i_r0+increase3
    i_r1<-i_r1+increase3
    y0<-y0+increase3
    y1<-y1+increase3
    y_text<-y0-increase1
    plotted_segs<-plotted_segs+1
    y0_line<-y0_line+increase1+increase3+increase1/((30/length(reads_to_plot))*plotted_segs)
  }
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- 15
  plot.params$data1height <- 10
  plot.params$data2height<-500
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  increase1<-0.2/length(reads_to_plot)
  increase2<-0.22/length(reads_to_plot)
  increase3<-0.3/length(reads_to_plot)
  
  i_r0<-0.1
  i_r1<-0.15
  y0<-0.1
  y1<-y0+increase1
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    for (i in 1:nrow(x))
    {
      if (x[i,]$strand=='+')
      {       
        i_col<-"orange"
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-"green"
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.05,lwd=15,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.05,lwd=13,lend='butt',ljoin='mitre', data.panel = 2)
      print(i_r0)
      print(i_r1)
      
      i_r0<-i_r0+increase1
      i_r1<-i_r1+increase1
    }
    
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    if (!is.na(x[1,]$vector_junction_locus))
    {
      j<-x[1,]$vector_junction_locus
      # if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
      #  j<-j+offset_arrow
      #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
      #j<-j-offset_arrow
      kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
    }
    if (nrow(x) > 1)
    {
      for (i in 2:nrow(x))
      {
        
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y0+0.0025*length(reads_to_plot) ,y1=y1,length=0.1, lty=2, data.panel = 2)
        y0<-y0+increase2
        y0_line<-y0_line+increase2
        y1<-y1+increase2
        if (!is.na(x[i,]$vector_junction_locus))
        {
          j<-x[i,]$vector_junction_locus
          #if (x[1,]$target_genome_position_wrt_junction=="R" & x[1,]$strand=='+')
          # j<-j+offset_arrow
          #else if (x[1,]$target_genome_position_wrt_junction=="L" & x[1,]$strand=='-')
          # j<-j-offset_arrow
          print(x[i,]$vector_junction_locus)
          kpLines(kp,chr = vector_chr,x=c(j,j), y=c(y0_line-increase3/(30/length(reads_to_plot)),y0_line+increase3/(18/length(reads_to_plot))),col='red', lwd=5, data.panel = 2)
        }
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end-(offset_arrow/2) else x[i,]$start+(offset_arrow/2)
        
      }
    }
    y0<-y0+increase1
    y1<-y1+increase1
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(y0,y0),lty=2,data.panel = 2,col='grey',lwd=2)
    i_r0<-i_r0+increase3
    i_r1<-i_r1+increase3
    y0<-y0+increase3
    y1<-y1+increase3
    y_text<-y0-increase1
    plotted_segs<-plotted_segs+1
    y0_line<-y0_line+increase1+increase3+increase1/((30/length(reads_to_plot))*plotted_segs)
  }
  dev.off()
}




###############################################################
#' @title Plot Read linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_read
#'
#' @param df_alm input df with alignment results
#' @param read_to_plot name of the read to plot
#' 
#' @return 
#' @description Plot the different alignments against the read
#' @note : 
#' @usage 
#' plot_linear_read(df,read_to_plot = "m64047_200324_183710/102237007/ccs",
#'        outfile_png = "read.png",outfile_pdf = "read.pdf")
plot_linear_read <- function(df_alm, read_to_plot, outfile_png, outfile_pdf,
                                      vector_chr = "chrV",
                                      height=8, width=8, res = 300,
                                      offset_arrow=110,
                                      color_vector_fwd="orange",
                                      color_vector_rev="green", 
                                      color_target_fwd="blue",
                                      color_target_rev="violet"
){
  library(karyoploteR)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction')
  i_r0<-1
  i_r1<-0.95
  for (i in 1:nrow(x))
  {
    if (x[i,]$strand=='+')
    {       
      ar_start<-x[i,]$start
      ar_end<-x[i,]$end-offset_arrow
      if (ar_start>=ar_end)
        ar_start<-ar_end-1
      if (x[i,]$genome==vector_chr)
        i_col<-color_vector_fwd
      else
        i_col<-color_target_fwd
    }
    else
    {        
      ar_start<-x[i,]$end
      ar_end<-x[i,]$start+offset_arrow
      if (ar_start<=ar_end)
        ar_start<-ar_end+1
      if (x[i,]$genome=='chrV')
        i_col<-color_vector_rev
      else
        i_col<-color_target_rev
    }
    kpArrows(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.1,lwd=15,lend='butt',ljoin='mitre')
    kpArrows(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.1,lwd=13,lend='butt',ljoin='mitre')
    if (x[i,]$junction==TRUE)
    {
        if(x[i,]$target_genome_position_wrt_junction=='R')
          j<-x[i,]$end
        else
          j<-x[i,]$start
     
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0,i_r0-2*(i_r0-i_r1)),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-0.1
    i_r1<-i_r1-0.1
  }
  legend(x = "bottomleft", fill = c(color_vector_fwd, color_vector_rev,color_target_fwd,color_target_rev), title=c('Strand and Genome'),legend = legend(x = "bottomleft", fill = c("orange", "green","blue","violet"), title=c('Strand and Genome'),legend = c("Vector (+)", "Vector (-)", "Target (+)","Target (-)")))
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction')
  i_r0<-1
  i_r1<-0.95
  for (i in 1:nrow(x))
  {
    if (x[i,]$strand=='+')
    {       
      ar_start<-x[i,]$start
      ar_end<-x[i,]$end-offset_arrow
      if (ar_start>=ar_end)
        ar_start<-ar_end-1
      if (x[i,]$genome==vector_chr)
        i_col<-"orange"
      else
        i_col<-"blue"
    }
    else
    {        
      ar_start<-x[i,]$end
      ar_end<-x[i,]$start+offset_arrow
      if (ar_start<=ar_end)
        ar_start<-ar_end+1
      if (x[i,]$genome=='chrV')
        i_col<-"green"
      else
        i_col<-"violet"
    }
    kpArrows(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', r0=i_r0, r1=i_r1,length=0.1,lwd=15,lend='butt',ljoin='mitre')
    kpArrows(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col, r0=i_r0, r1=i_r1,length=0.1,lwd=13,lend='butt',ljoin='mitre')
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0,i_r0-2*(i_r0-i_r1)),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-0.1
    i_r1<-i_r1-0.1
  }
  legend(x = "bottomleft", fill = c(color_vector_fwd, color_vector_rev,color_target_fwd,color_target_rev), title=c('Strand and Genome'),legend = legend(x = "bottomleft", fill = c("orange", "green","blue","violet"), title=c('Strand and Genome'),legend = c("Vector (+)", "Vector (-)", "Target (+)","Target (-)")))
  dev.off()
}




###############################################################
#' @title Plot rearrangements linearized
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_rearrangement
#' @docType methods
#' @aliases plot_linear_rearrangement
#'
#' @param df_alm input df with alignment results
#' @param reads_to_plot vector with names of the reads to plot
#' @param cytoband_file path to cytoband file NOTE: the file must have the following
#'                      header: chr	start	end	name	gieStain 
#' @param vector_length vector genome length
#' @param annotation_file bed file with annotation
#' @param color_bed_file bed file with rgb colors, include all intervals, see example one
#' @param lwd_black width of the external arrow (this is used to make the arrow border, should be always >= lwd_color)
#' @param lwd_color width of the arrow
#' @param offset_arrow once you have set the two lwd pars you can use this parameter to be sure the arrow end in the right position
#' to do that change this parameter according to red lines which indicates the junction locus. You probably need to change this parameters for different
#' set of reads to plot
#' @return 
#' @description Plot the vector portions of the read
#' @note : 
#' @usage 
#' plot_linear_rearrangement(df = alm_reads_for_rearrangements_withstatsbyread, reads_to_plot = reads_to_plot<-c("m64047_200324_183710/83362987/ccs" ,"m64047_200324_183710/102237007/ccs"),
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' outfile_pdf = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.pdf", sep = ""), 
#' vector_chr = "chrV", offset_arrow=110, lwd_black = 160, lwd_color=148,
#' cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
#' annotation_file="source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.bed",
#' )
plot_linear_rearrangement_impr <- function(df_alm, reads_to_plot, outfile_png, outfile_pdf,
                                      cytoband_file, annotation_file,color_bed_file, vector_chr = "chrV",
                                      vector_length=6000,
                                      height=8, width=8, res = 300,
                                      offset_arrow=110, lwd_black = 160, lwd_color=148,
                                      color_vector_fwd="orange", color_vector_rev="green",
                                      ideogram_height = 5, aav_colored_thickness = 0.3
){
  library(karyoploteR)
  library(bezier)
  library(ggplot2)
  all_slice<-df_alm[is.element(df_alm$name, reads_to_plot),]
  all_slice<-all_slice[all_slice$chr==vector_chr,]
  lwd_black<-lwd_black/nrow(all_slice)
  lwd_color<-lwd_color/nrow(all_slice)
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  
  
 
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- ideogram_height
  plot.params$data1height <- 15
  plot.params$data2height<-200
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  for(i in seq(1:(nrow(bed_color)-1)))
  {
    rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
    r<-rgb_string[1]
    g<-rgb_string[2]
    b<-rgb_string[3]
    kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i+1,]$V2,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  }
  i<-i+1
  rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
  r<-rgb_string[1]
  g<-rgb_string[2]
  b<-rgb_string[3]
  kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i,]$V3,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  
  increase<-0.65/nrow(all_slice)
  
  i_r0<-0.1
  i_r1<-i_r0+increase
  y0<-0.1
  y1<-y0+increase
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  y_r0_ar<-i_r0
  y_r1_ar<-i_r1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    #y_r0_ar<-y_r0_ar+(0.12/length(reads_to_plot))*(plotted_segs-1)+0.01
    #y_r1_ar<-y_r1_ar+0.007*(plotted_segs-1)+0.01
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    for (i in 1:nrow(x))
    {
      
      if (x[i,]$strand=='+')
      {       
        i_col<-color_vector_fwd
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-color_vector_rev
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      #draw segments
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', length=0.05,lwd=lwd_black,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col,length=0.05,lwd=lwd_color,lend='butt',ljoin='mitre', data.panel = 2)
      if (!is.na(x[i,]$vector_junction_locus))
      {
        j<-x[i,]$vector_junction_locus
        kpLines(kp,chr = vector_chr,x=c(j,j), y=c(i_r0-increase/3,i_r0+increase/3),col='red', lwd=5, data.panel = 2)
      }
      #draw connecting curves
      if (i>1){
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        t<-seq(0,1,length=100)
        bz<-bezier(t,p=list(c(arrow_start,arrow_start,arrow_start,arrow_end,arrow_end,arrow_end),seq(i_r0-increase,i_r0, length=6)))
        kpLines(kp,chr="chrV",x=bz[,1],y=bz[,2],data.panel = 2, lty=1)
        #kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y_r0_ar ,y1=y_r0_ar+increase1-0.01,length=0.1, lty=2, data.panel = 2)
        
        
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end else x[i,]$start
        
        
      }
      i_r0<-i_r0+increase
      i_r1<-i_r1+increase
      
    }
    
    
    
    
    
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(i_r0,i_r0),lty=2,data.panel = 2,col='grey',lwd=2)
    y_text<-i_r0+increase/2
    i_r0<-i_r0+increase*1.5
    i_r1<-i_r1+increase*1.5
    
    
    
  }
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  
  
  annot<-read.table(file = annotation_file, sep = '\t', header = FALSE)
  colnames(annot)<-c('chr','start','end','feature')
  annot$value<-0.1
  custom.genome <- toGRanges(data.frame(chr=c(vector_chr), start=c(1), end=c(vector_length)))
  custom.cytobands <- toGRanges(cytoband_file)
  annot_gr<-toGRanges(annot)
  
  plot.params <- getDefaultPlotParams(plot.type=3)
  plot.params$ideogramheight <- ideogram_height
  plot.params$data1height <- 15
  plot.params$data2height<-200
  kp <- plotKaryotype(genome = custom.genome,cytobands = custom.cytobands, plot.type = 3, plot.params = plot.params, )
  kpText(kp,data=annot_gr,labels = annot$feature,data.panel = 1)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  for(i in seq(1:(nrow(bed_color)-1)))
  {
    rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
    r<-rgb_string[1]
    g<-rgb_string[2]
    b<-rgb_string[3]
    kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i+1,]$V2,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  }
  i<-i+1
  rgb_string<-strsplit(bed_color[i,]$rgb," ",fixed=TRUE)[[1]]
  r<-rgb_string[1]
  g<-rgb_string[2]
  b<-rgb_string[3]
  kpRect(kp,chr="chrV",x0=bed_color[i,]$V2,x1=bed_color[i,]$V3,y0=0.5,y1=0.5+aav_colored_thickness,col=rgb(r,g,b,maxColorValue = 255),data.panel = 1)
  
  increase<-0.65/nrow(all_slice)
  
  i_r0<-0.1
  i_r1<-i_r0+increase
  y0<-0.1
  y1<-y0+increase
  y0_line<-y0
  x_text<-3000
  y_text<-0.03
  plotted_segs<-1
  y_r0_ar<-i_r0
  y_r1_ar<-i_r1
  for (read in 1:length(reads_to_plot))
  {
    slice_t <- df_alm[which(df_alm$name == reads_to_plot[read]),]
    slice_t <- slice_t[order(slice_t$query_start),]
    x<-slice_t[slice_t$chr==vector_chr,c('chr','start','end','strand','vector_junction_locus',
                                         'target_genome_position_wrt_junction')]
    
    #y_r0_ar<-y_r0_ar+(0.12/length(reads_to_plot))*(plotted_segs-1)+0.01
    #y_r1_ar<-y_r1_ar+0.007*(plotted_segs-1)+0.01
    kpText(kp,x=x_text,y = y_text,chr=vector_chr,labels = c(reads_to_plot[read]),data.panel = 2)
    arrow_start <- if (x[1,]$strand=='+') x[1,]$end-(offset_arrow/2) else x[1,]$start+(offset_arrow/2)
    for (i in 1:nrow(x))
    {
      
      if (x[i,]$strand=='+')
      {       
        i_col<-color_vector_fwd
        ar_start<-x[i,]$start
        ar_end<-x[i,]$end-offset_arrow
        if (ar_start>=ar_end)
          ar_start<-ar_end-1
      }
      else
      {        
        i_col<-color_vector_rev
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start+offset_arrow
        if (ar_start<=ar_end)
          ar_start<-ar_end+1
      }
      #draw segments
      kpArrows(kp,chr=vector_chr,x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0, col='black', length=0.05,lwd=lwd_black,lend='butt',ljoin='mitre', data.panel = 2)
      kpArrows(kp,chr=vector_chr,x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0, col=i_col,length=0.05,lwd=lwd_color,lend='butt',ljoin='mitre', data.panel = 2)
      if (!is.na(x[i,]$vector_junction_locus))
      {
        j<-x[i,]$vector_junction_locus
        kpLines(kp,chr = vector_chr,x=c(j,j), y=c(i_r0-increase/3,i_r0+increase/3),col='red', lwd=5, data.panel = 2)
      }
      #draw connecting curves
      if (i>1){
        arrow_end <- if (x[i,]$strand=='+') x[i,]$start else x[i,]$end
        t<-seq(0,1,length=100)
        bz<-bezier(t,p=list(c(arrow_start,arrow_start,arrow_start,arrow_end,arrow_end,arrow_end),seq(i_r0-increase,i_r0, length=6)))
        kpLines(kp,chr="chrV",x=bz[,1],y=bz[,2],data.panel = 2, lty=1)
        #kpArrows(kp,chr=vector_chr,x0=arrow_start,x1=arrow_end, y0=y_r0_ar ,y1=y_r0_ar+increase1-0.01,length=0.1, lty=2, data.panel = 2)
        
        
        arrow_start <- if (x[i,]$strand=='+') x[i,]$end else x[i,]$start
        
        
      }
      i_r0<-i_r0+increase
      i_r1<-i_r1+increase
      
    }
    
    
    
    
    
    kpLines(kp,chr=vector_chr,x=c(0,6000),y=c(i_r0,i_r0),lty=2,data.panel = 2,col='grey',lwd=2)
    y_text<-i_r0+increase/2
    i_r0<-i_r0+increase*1.5
    i_r1<-i_r1+increase*1.5
    
    
    
  }
dev.off()
}

###############################################################
#' @title Plot Read linearized with bed colors
#' 
#' @author Carlo Cipriani
#' @details version 0.1, 2021-06-3
#'
#' @rdname plot_linear_read_colorized
#' @docType methods
#' @aliases plot_linear_read
#'
#' @param df_alm input df with alignment results
#' @param read_to_plot name of the read to plot
#' @param space_btw_rects space between two rects, default = 0.05
#' @param rect_thickness thickness of a rect, default = 0.05
#' @param cex_legend scaling factor of legend, default = 0.5 
#' @return 
#' @description Plot the different alignments against the read using the value in the bed as colors
#' @note : 
#' @usage 
#' plot_linear_read(df,read_to_plot = "m64047_200324_183710/102237007/ccs",
#'        outfile_png = "read.png",outfile_pdf = "read.pdf")
plot_linear_read_colorized <- function(df_alm, read_to_plot, outfile_png, outfile_pdf,
                             vector_chr = "chrV",
                             height=8, width=8, res = 300,
                             offset_arrow=110,
                             color_bed_file,
                             color_target_fwd="blue",
                             color_target_rev="violet",
                             space_btw_rects=0.05,
                             rect_thickness = 0.05,
                             cex_legend = 0.5
){
  library(karyoploteR)
  bed_color <- read.table(file = color_bed_file, sep = '\t', header = FALSE)
  sr<-bed_color[1,]
  sr$V1<-"start"
  sr$V2<-0
  sr$V3<-0
  se<-bed_color[nrow(bed_color),]
  se$V1<-"end"
  se$V2<-se$V3+1
  se$V3<-se$V3+1
  bed_color<-rbind(sr,bed_color,se)
  bed_color['rgb']<-gsub(",","",bed_color$V9)
  
  png(file = outfile_png, height=height, width=width, units = "in", res = res)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction','start','end')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction','orig_start','orig_end')
  space <- space_btw_rects
  thickness <- rect_thickness
  i_r0<-1
  i_r1<-i_r0 - thickness
  for (i in 1:nrow(x))
  {
    ar_start<-x[i,]$start
    ar_end<-x[i,]$end
    ltype = 1
    if(x[i,]$genome==vector_chr)
    {
      if(x[i,]$strand=='-')
      {
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start
      }
      
      range_start<-x[i,]$orig_start
      range_new_end<-x[i,]$orig_start
      range_end<-x[i,]$orig_end
      for (jj in seq(1:(nrow(bed_color)-1))){
        if  (range_start>=bed_color[jj,]$V2 & range_start<=bed_color[jj+1,]$V2)
        {
          rgb_string<-strsplit(bed_color[jj,]$rgb," ",fixed=TRUE)[[1]]
          r<-rgb_string[1]
          g<-rgb_string[2]
          b<-rgb_string[3]
          i_col<-rgb(r,g,b,maxColorValue = 255)
          range_new_end<-bed_color[jj+1,]$V2-1
          span<-range_new_end-range_start
          if (range_new_end>=range_end)
          {
            
            range_new_end<-range_end
            span<-range_new_end-range_start
            kpRect(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
            ar_start<-ar_start+span
            break
          }
          else
          {
            
            if(x[i,]$strand=='+')
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start+span, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
              ar_start<-ar_start+span
            }
            
            else
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start-span, y0=i_r0 ,y1=i_r0-thickness, col=i_col, lty=ltype)
              ar_start<-ar_start-span
            }
            
          }
          
          range_start<-bed_color[jj+1,]$V2
          
        }
      }
    }
    else
    {
      if (x[i,]$strand=='+')
      {
        i_col = color_target_fwd
        ltype=1
      }
      else
      {
        i_col = color_target_rev
        ltype=2
      }
      kpRect(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
    }
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0, i_r0-thickness),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-thickness
    i_r1<-i_r1-thickness
    i_r0<-i_r0-space
    i_r1<-i_r1-space
  }
  legend(x = "bottomleft", lty = c(1,2), legend=c("+", "-"), title=c('Target Strand'), cex=cex_legend)
  dev.off()
  
  pdf(file = outfile_pdf, height=height, width=width)
  slice_t <- df_alm[which(df_alm$name == read_to_plot),]
  slice_t <- slice_t[order(slice_t$query_start),]
  end<-slice_t$query_end[nrow(slice_t)]
  custom.genome <- toGRanges(data.frame(chr=c("read"), start=c(1), end=c(end)))
  plot.params <- getDefaultPlotParams(plot.type=5)
  plot.params$ideogramheight <- 5
  plot.params$data2height<-100
  kp <- plotKaryotype(genome = custom.genome, plot.type = 5, plot.params = plot.params)
  kpAddMainTitle(kp, main=read_to_plot)
  
  
  x<-slice_t[,c('name','query_start','query_end','strand','chr',
                'junction','target_genome_position_wrt_junction','start','end')]
  x$name<-'read'
  colnames(x)<-c('chr','start','end','strand','genome','junction','target_genome_position_wrt_junction','orig_start','orig_end')
  space <- space_btw_rects
  thickness <- rect_thickness
  i_r0<-1
  i_r1<-i_r0 - thickness
  for (i in 1:nrow(x))
  {
    ar_start<-x[i,]$start
    ar_end<-x[i,]$end
    ltype = 1
    if(x[i,]$genome==vector_chr)
    {
      if(x[i,]$strand=='-')
      {
        ar_start<-x[i,]$end
        ar_end<-x[i,]$start
      }
      
      range_start<-x[i,]$orig_start
      range_new_end<-x[i,]$orig_start
      range_end<-x[i,]$orig_end
      for (jj in seq(1:(nrow(bed_color)-1))){
        if  (range_start>=bed_color[jj,]$V2 & range_start<=bed_color[jj+1,]$V2)
        {
          rgb_string<-strsplit(bed_color[jj,]$rgb," ",fixed=TRUE)[[1]]
          r<-rgb_string[1]
          g<-rgb_string[2]
          b<-rgb_string[3]
          i_col<-rgb(r,g,b,maxColorValue = 255)
          range_new_end<-bed_color[jj+1,]$V2-1
          span<-range_new_end-range_start
          if (range_new_end>=range_end)
          {
            
            range_new_end<-range_end
            span<-range_new_end-range_start
            kpRect(kp,chr='read',x0=ar_start, x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
            ar_start<-ar_start+span
            break
          }
          else
          {
            
            if(x[i,]$strand=='+')
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start+span, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
              ar_start<-ar_start+span
            }
            
            else
            {
              kpRect(kp,chr='read',x0=ar_start, x1=ar_start-span, y0=i_r0 ,y1=i_r0-thickness, col=i_col, lty=ltype)
              ar_start<-ar_start-span
            }
            
          }
          
          range_start<-bed_color[jj+1,]$V2
          
        }
      }
    }
    else
    {
      if (x[i,]$strand=='+')
      {
        i_col = color_target_fwd
        ltype=1
      }
      else
      {
        i_col = color_target_rev
        ltype=2
      }
      kpRect(kp,chr='read',x0=ar_start,x1=ar_end, y0=i_r0 ,y1=i_r0-thickness, col=i_col,lty=ltype)
    }
    if (x[i,]$junction==TRUE)
    {
      if(x[i,]$target_genome_position_wrt_junction=='R')
        j<-x[i,]$end
      else
        j<-x[i,]$start
      
      kpLines(kp,chr = 'read',x=c(j,j), y=c(i_r0, i_r0-thickness),col='red', lwd=5, data.panel = 2)
      
    }
    i_r0<-i_r0-thickness
    i_r1<-i_r1-thickness
    i_r0<-i_r0-space
    i_r1<-i_r1-space
  }
  legend(x = "bottomleft", lty = c(1,2), legend=c("+", "-"), title=c('Target Strand'), cex=cex_legend)
  dev.off()
  
}




###############################################################
#' @title Acquire CIS from folder
#' 
#' @author Andrea
#' @details version 0.1, 2021-10-30
#'
#' @rdname acquire_CIS_infolder
#' @docType methods
#' @aliases acquire_CIS_infolder
#'
#' @param folder input folder
#' @param prefix file suffix, typically trial ID as patient code prefix
#' @param suffix file suffix, tsv
#' @param sep file separator
#' @return 
#' @description Plot the vector portions of the read
#' @note : 
#' @usage 
#' plot_linear_rearrangement(df = alm_reads_for_rearrangements_withstatsbyread, reads_to_plot = reads_to_plot<-c("m64047_200324_183710/83362987/ccs" ,"m64047_200324_183710/102237007/ccs"),
#' outfile_png = paste0(dest_dir, "results/", ProjectID, ".plot_alm_reads_stats.linearized.InVivo.V-V-V-V-X.01.png", sep = ""), 
#' vector_chr = "chrV", cytoband_file = "source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.cytoband",
#' annotation_file="source/metadata/CAG_Tomato_withBackbone/AAV-CAG-tdTomato.withBackbone.bed")
################################################################
acquire_CIS_infolder <- function(folder, 
                                 prefix,
                                 suffix = ".tsv",
                                 study = "GTtrial",
                                 sep = "\t"
                                 ){
  cislist <- list.files(folder, pattern=paste0(prefix, ".*", "\\", suffix, "$"))
  message(paste0("[AP]\tAcquire CIS from all patients in the folder: ", folder))
  out_df <- NULL
  for (i in cislist) {
    message(paste0("[AP]\t\t-> Processing file: ", i))
    patient_id <- strsplit(i, "\\.")[[1]][1]
    this_cis <- read.csv(file = paste0(folder, i), header=T, fill=T, check.names = FALSE, sep = sep)
    this_cis$SubjectID <- patient_id
    this_cis$ProjectID <- study
    if (length(out_df)>0) {
      out_df <- rbind(out_df, this_cis)
    } else {
      out_df <- this_cis
    } # if
  } # for
  return ( out_df )
}



###############################################################
#' @title Acquire sharing data from ISAnalytics prototype (v0.4 2015)
#' 
#' @author Andrea
#' @details version 0.1, 2021-12-10
#'
#' @rdname getCD34outputFromSharing
#' @docType methods
#' @aliases getCD34outputFromSharing
#'
#' @param file_stats stats file, the output of descriptive stats from the new ISAnalytics package
#' @param trial_wd file suffix, typically trial ID as patient code prefix
#' @param nint_padding padding timepoint, number of 0s, default = 2
#' @param input_results_folder input file folder
#' @param file_pattern_search file pattern to search. Default: ".filter3.purity.stem_population.sharing.ISA.4.tsv"
#' @param patients_to_exclude if any, vector of patients to exclude. 
#' @param colname_to_convertto_cellmarker if any, vector of patients to exclude. 
#' @param quiet print verbose if True.
#' @return shared34_stats as dataframe of results
#' @description Plot the vector portions of the read
#' @note : IMPORTANT: file names MUST start with patient ID. Strict requirement.
#' @usage 
################################################################
getCD34outputFromSharing <- function(file_stats,
                                     trial_wd, 
                                     nint_padding = 2,
                                     input_results_folder,
                                     file_pattern_search = ".filter3.purity.stem_population.sharing.ISA.4.tsv",
                                     patients_to_exclude = c(),
                                     colname_to_convertto_cellmarker = "SuperGroup",
                                     quiet = F
) {
  shared34_stats <- NULL
  if (!quiet) {message(paste("[AP]\tImporting stats"))}
  stats_summary_df_full_descriptivestats_complete_postQF <- read.csv(file = file_stats, header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
  stats_summary_df_full_descriptivestats_complete_postQF$NumIS <- stats_summary_df_full_descriptivestats_complete_postQF$nIS
  if (!("TimePoint" %in% colnames(stats_summary_df_full_descriptivestats_complete_postQF)) & ("TimepointMonths" %in% colnames(stats_summary_df_full_descriptivestats_complete_postQF)) ) {
    stats_summary_df_full_descriptivestats_complete_postQF$TimePoint <- stats_summary_df_full_descriptivestats_complete_postQF$TimepointMonths
  } else {
    message(paste("[AP]\tERROR:\tNo columns for TIMEPOINT"))
  }
  
  stats <- stats_summary_df_full_descriptivestats_complete_postQF[which(!(is.na(stats_summary_df_full_descriptivestats_complete_postQF$NumIS)) & 
                                                                          stats_summary_df_full_descriptivestats_complete_postQF$NumIS > 5 ), ]
  if (colname_to_convertto_cellmarker %in% colnames(stats) & !("CellMarker" %in% colnames(stats))) {
    stats$CellMarker <- stats[,colname_to_convertto_cellmarker]
  }
  stats$gdf <- apply(stats[c("SubjectID", "CellMarker", "Tissue", "TimePoint")], 1, function(x) {paste(x[1], x[2], x[3], str_pad(round(as.numeric(x[4])), nint_padding, pad = "0"), sep = '_')})
  rownames(stats) <- stats$gdf
  stats$groupby_marker_tissue_fumonths_wet <- stats$gdf
  
  label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet <- stats
  # update metadata
  label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$MarkerTissue <- apply(label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet[c("CellMarker", "Tissue")], 1, function(x) {paste(x[1], x[2], sep = '_')}) 
  label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$FollowUp <- as.numeric(as.character(label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$TimePoint))
  # erythroid_markers <- c("CD36", "GLY")
  # label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$BloodLineage <- ifelse(
  #   label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$CellMarker %in% erythroid_markers, "Erythroid", "Other"
  # )
  
  # input_results_folder <- paste0(basedir, "11.stem_population/202110")
  infiles_fullpath <- list.files(paste0(input_results_folder), pattern = file_pattern_search, full.names = T)
  infiles_name <- list.files(paste0(input_results_folder), pattern = file_pattern_search, full.names = F)
  if (!quiet) {message(paste0("[AP]\tFiles to parse:\n", paste(infiles_name, collapse = "\n\t--> ")))}
  for (infile in infiles_fullpath) {
    patient_id <- as.character(as.data.frame(t(as.data.frame(lapply(strsplit(as.character(gsub(paste0(file_pattern_search), "", infile)), ".", fixed = TRUE), function(x) {c(x)}))))[,4])
    message(paste("[AP]\tProcessing patient", patient_id, " file:", infile))
    if (!(patient_id %in% patients_to_exclude)) {
      annotation_cols <- c("IS_code", "seqnames", "start", "end", "width", "strand", "closest_gene", "strand_gene", "in.progenitor", "sharing" )
      # get sharing IS
      patientSharedIS <- read.csv(file = paste0(infile), 
                                  header=T, fill=T, sep='\t',
                                  check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
      # rownames(patient_stats) <- gsub("_SLiM", "", gsub("_LAM", "", gsub(paste(patient_id, "_", sep = ''), "", rownames(patient_stats))))
      patientSharedISslice <- patientSharedIS[which(patientSharedIS$in.progenitor == "TRUE" & patientSharedIS$sharing > 1),] # select only sharing from CD34 BM
      patientSharedISslice_stats <- data.frame(
        "SharedCD34IS_NumIS" = apply(patientSharedISslice[setdiff(colnames(patientSharedISslice), annotation_cols)], 2, function(x) {length(x[!is.na(x)])})
        # "" = apply(gdf_iss_nona[setdiff(colnames(gdf_iss_nona), id_cols)], 2, function(x) {sum(x, na.rm = T)}),
        # "Hindex" = c(diversity(t(gdf_iss_nona[setdiff(colnames(gdf_iss_nona), id_cols)])))
      )
      # slice metadata and adapt names to matrix
      patientSharedISmeta <- label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet[which(label_metadata_agg_metadata_groupby_marker_tissue_fumonths_wet$SubjectID == patient_id ), ]
      rownames(patientSharedISmeta) <- gsub("_SLiM", "", gsub("_LAM", "", gsub(paste(patient_id, "_", sep = ''), "", rownames(patientSharedISmeta))))
      # acquire cumulative sharing IS
      patientSharedIScum <- getBinaryCumulative_ProgressiveUnion_byColumnsName(df = patientSharedISslice, metadata_df = patientSharedISmeta[intersect(colnames(patientSharedISslice), rownames(patientSharedISmeta)),], key_field = "MarkerTissue", starting_data_col_index = 11, number_of_last_cols_to_remove = 0, decreasing_sort = F)
      # do some stats
      patientSharedIScum_stats <- data.frame(
        "Cumulative_NumIS" = apply(patientSharedIScum, 2, function(x) {length(x[x>0])})
        # "" = apply(gdf_iss_nona[setdiff(colnames(gdf_iss_nona), id_cols)], 2, function(x) {sum(x, na.rm = T)}),
        # "Hindex" = c(diversity(t(gdf_iss_nona[setdiff(colnames(gdf_iss_nona), id_cols)])))
      )
      
      # merge results shared IS
      patient_iss_stats <- NULL
      patient_iss_stats <- merge(x = patientSharedISmeta, y = patientSharedISslice_stats, by = 0, all.y = TRUE)
      rownames(patient_iss_stats) <- patient_iss_stats$Row.names
      patient_iss_stats <- patient_iss_stats[, !(names(patient_iss_stats) %in% c("Row.names"))]
      patient_iss_stats <- patient_iss_stats[which(!is.na(patient_iss_stats$SubjectID)),]
      
      # # now with cumulative
      patient_iss_stats <- merge(x = patient_iss_stats, y = patientSharedIScum_stats, by = 0, all.y = TRUE)
      rownames(patient_iss_stats) <- patient_iss_stats$Row.names
      patient_iss_stats <- patient_iss_stats[, !(names(patient_iss_stats) %in% c("Row.names"))]
      # add number of all progenitors shared
      nIS_in_progenitor <- nrow(patientSharedISslice)
      patient_iss_stats$Progenitor_nIS <- nIS_in_progenitor
      # do some ratios
      patient_iss_stats$PercNIS_onOverallSharedCD34BM <- patient_iss_stats$NumIS * 100 / patient_iss_stats$Progenitor_nIS
      patient_iss_stats$PercCumulativeNIS_onOverallSharedCD34BM <- patient_iss_stats$Cumulative_NumIS * 100 / patient_iss_stats$Progenitor_nIS
      
      # reassign rown names to assemble the final dataset
      # rownames(patient_iss_stats) <- patient_iss_stats$groupby_marker_tissue_fumonths_wet
      rownames(patient_iss_stats) <- patient_iss_stats$gdf_names
      if (length(shared34_stats) > 0) {
        shared34_stats <- rbind(shared34_stats, patient_iss_stats)
      } else {
        shared34_stats <- patient_iss_stats
      } # if (length(shared34_stats) > 0)
    } else {
      message(paste("[AP]\tThis patient", patient_id, " will be skept since in the list of the excluded patients."))
    } # if (!(patient_id %in% patients_to_exclude))
    
  } # for (infile in infiles_fullpath)
  
  shared34_stats[is.na(shared34_stats)] <- NA
  rownames(shared34_stats) <- shared34_stats$gdf
  return (shared34_stats)
}

###############################################################
#' @title Reannotate chr R from graphDB data
#' 
#' @author Andrea
#' @details version 0.1, 2021-12-22
#'
#' @rdname reannotate_Repeat_IS
#' @docType methods
#' @aliases reannotate_Repeat_IS
#'
#' @param df graph db IS data
#' @param id_coils columns of annotations
#' @return df of IS reannotated
#' @description Change chrR annotations
#' @note : 
#' @usage 
################################################################
reannotate_Repeat_IS <- function(df, id_cols = c("Chr", "Pos", "Strand")) {
  # slice matrix with only repeats
  slice_R <- df[which(df[id_cols[1]] == "R" & df[id_cols[2]] == -1),]
  slice_notR <- df[which( !(df[id_cols[1]] == "R" & df[id_cols[2]] == -1) ),]
  slice_R[id_cols[2]] <- seq(1,nrow(slice_R))
  return( rbind(slice_notR, slice_R) )
}

###############################################################
#' @title Rewrite flag matrixes to analyze source of IS
#' 
#' @author Andrea
#' @details version 0.1, 2021-12-22
#'
#' @rdname rewriteFlagMatrix_forSourceIS
#' @docType methods
#' @aliases rewriteFlagMatrix_forSourceIS
#'
#' @param filepath path flag matrix file
#' @param id_coils columns of annotations
#' @param score_to_isolate int to isolate in the matrix df
#' @return df of IS isolated
#' @description Todo
#' @note : 
#' @usage 
################################################################
rewriteFlagMatrix_forSourceIS <- function(filepath, 
                                          id_cols = c("chr", "integration_locus", "strand", "GeneName", "GeneStrand"),
                                          score_to_isolate = 1,
                                          SubjectID,
                                          CellMarker = "Multi",
                                          Tissue = "Tissue"
                                          ) {
  mfl_multi <- read.csv(file = filepath, header=TRUE, fill=T, sep='\t', check.names = FALSE, na.strings = c("NONE", "NA", "NULL", "NaN", ""))
  mfl_multi[mfl_multi > score_to_isolate] <- NA
  mfl_multi[mfl_multi < score_to_isolate] <- NA
  
  mfl_multi_slice <- compactDfByRows(df = mfl_multi, data_columns = colnames(mfl_multi), annotation_columns = c())
  mfl_multi_slice_molten <- melt(as.matrix(mfl_multi), na.rm = T, value.name = "MultiCode", varnames = c("IS_ID", "TimePoint"))
  mfl_multi_slice_molten_details <- as.data.frame(t(as.data.frame(lapply(strsplit(as.character(mfl_multi_slice_molten$IS_ID), "_", fixed = TRUE), function(x) {c(x)}))))
  names(mfl_multi_slice_molten_details) <- id_cols
  mfl_multi_slice_molten_details$SubjectID <- "MLD02"
  mfl_multi_slice_molten_details$CellMarker <- "Multi"
  mfl_multi_slice_molten_details$Tissue <- "Tissue"
  mfl_multi_slice_molten <- cbind(mfl_multi_slice_molten_details, mfl_multi_slice_molten[-c(grep("IS_ID", colnames(mfl_multi_slice_molten)))])
  # write.table(mfl_multi_slice_molten, file = paste0("/Users/calabria.andrea/Dropbox (HSR Global)/Clinical Trial MLD/analyses/13.hsc_profiles/202112/MyBTEry34_allMarkers_noLy_BM_SC3/MLD02.hlfu_excl34MyBTEry.flag_matrix_timeShared_relabeled_correctedDC.forSourceIS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = '')
  
  return( mfl_multi_slice_molten )
}


