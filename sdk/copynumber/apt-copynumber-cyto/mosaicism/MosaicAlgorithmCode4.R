
gains.boundaries  <- c(  0.08764945,  0.15380349,  0.21465931,  0.27100300)
losses.boundaries <- c( -0.08293345, -0.17551812, -0.28048196, -0.40165383)

debug_current_chr      <<- "ERROR-should-be-set"
debug_current_basename <<- "ERROR-should-be-set"
debug_current_filename <<- "ERROR-should-be-set"

debug_write_segs_cnt <<- 0

write_segs <- function (segs,f_suffix,msg="") {
  debug_write_segs_cnt <<- debug_write_segs_cnt + 1
  filename <- paste(debug_current_basename,
                    "--",
                    formatC(debug_write_segs_cnt,format="d",width=2,flag="0"),
                    "-",
                    f_suffix,
                    sep="")
  write.table(segs,file=filename,row.names=FALSE,sep="\t")
  cat("### wrote segs:",filename,msg,"\n")
}

##
## x are log ratios
## gains.boundaries and losses.bounaries tell the algorithm where to make cut-offs (data derived fixed constants are above)
## bandwidth tells how many markers to average in the running average (note the boundary adjustments
## have been calibrated to the 6000 marker bandwidth here)
##

mosaic.segmentation.algorithm <- function(x, gains.boundaries, losses.boundaries,bandwidth=6000){

  # This function determines which class the mosaic falls into
  get.mos.prop <- function(most.extreme,boundaries){
    if (abs(most.extreme) >= abs(boundaries[4])){
      return(1)
    } else if (abs(most.extreme) >= abs(boundaries[3])){
      return(0.7)
    } else if (abs(most.extreme) >= abs(boundaries[2])){
      return(0.5)
    } else if (abs(most.extreme) >= abs(boundaries[1])){
      return(0.3)
    } else {
      return(0)
    }
  }

  # This function determines how to adjust the segment boundaries
  GetBoundaryAdjustment <- function(prop){
    if (prop == -1){
      return(2165)
    } else if (prop == -0.7){
      return(1777)
    } else if (prop == -0.5){
      return(1054)
    } else if (prop == -0.3){
      return(-1780)
    } else if (prop == 0.3){
      return(-2128)
    } else if (prop == 0.5){
      return(209)
    } else if (prop == 0.7){
      return(909)
    } else if (prop == 1.0){
      return(1342)
    }
  }

  ### new subroutine
  join.overlaps <- function(x) {
    write_segs(x,"join0")   # mosaicism:debugging

    # first find the overlaps
    segment.label <- rep(0,nrow(x))
    cur.label <- 1
    segment.label[1] <- cur.label
    row.ind <- 1
    while (row.ind < nrow(x)){
      cur.seg <- x[row.ind,]
      next.seg <- x[row.ind+1,]
      if (cur.seg[3] == next.seg[3]){
        # adjoining segments
        segment.label[row.ind +1] <- cur.label
      } else {
        cur.label <- cur.label + 1
        segment.label[row.ind +1] <- cur.label
      }
      row.ind <- row.ind +1
    }

    # now go through and join up adjacent regions
    cleaned.segmentation <- NULL

    start <- 1
    end <- start
    while (end < nrow(x)){
      end <- end+1
      if (segment.label[start] != segment.label[end]){
        cur.seg <- c(x[start,1],x[end-1,2],x[start,3])
        cleaned.segmentation <- c(cleaned.segmentation,list(cur.seg))
        start <- end
      }
    }
    cur.seg <- c(x[start,1],x[end,2],x[start,3])
    cleaned.segmentation <- c(cleaned.segmentation,list(cur.seg))
##    browser()

    # cleaned.segmentation is a list
    # rbind puts it onto a frame.
    x <- do.call("rbind",cleaned.segmentation)

    write_segs(x,"join1")   # mosaicism:debugging

    ## now check for any overlap segments
    ## at this point only those differing in call should exist
    cleaned.segmentation <- NULL
    row.ind <- 1
    cur.seg <- x[row.ind,]
    next.seg <- cur.seg
    while (row.ind < nrow(x)){
      cur.seg <- x[row.ind,]
      next.seg <- x[row.ind+1,]
      if (cur.seg[2] >= next.seg[1]){
        # overlaping segments
        if (cur.seg[1] < next.seg[1]){
          cur.seg[2] <- next.seg[1] -1
          cleaned.segmentation <- c(cleaned.segmentation,list(cur.seg))
        }
      } else {
        # no overlap
        cleaned.segmentation <- c(cleaned.segmentation,list(cur.seg))
      }
      row.ind <- row.ind +1
    }
    cleaned.segmentation <- c(cleaned.segmentation,list(next.seg))

    x <- do.call("rbind",cleaned.segmentation)
    write_segs(x,"join2")   # mosaicism:debugging
    x
  }

  # This function cleans out really small segments
  # then adjusts the boundaries to get final segmentation
  clean.segments <- function(x, min.size=20){
    write_segs(x,"clean1")
    cleaned.segmentation <- NULL

    # X: [ start, end, ???, ... ]
    # First get rid of really small segments
    region.size <- x[,2] - x[,1]
    x[region.size < min.size,3] <- 0
    ##print(x)

    Are.zero <- x[,3] == 0

    row.ind <- 1
    while (row.ind <= nrow(x)){
      cat("row.ind=",row.ind,"\n")
      cat("cleaned.segmentation=",paste(cleaned.segmentation,sep="\t"),"\n")
      #
      if (!Are.zero[row.ind]){
        cleaned.segmentation <- c(cleaned.segmentation,list(x[row.ind,]))
        row.ind <- row.ind + 1
      } else {
        modup <- 0
        while(Are.zero[row.ind + (modup + 1)] & (row.ind + (modup)) < nrow(x)){
          modup <- modup+1
        }
        cleaned.segmentation <- c(cleaned.segmentation,list(c(x[row.ind,1],x[row.ind+(modup),2],0)))
        row.ind <- row.ind + modup + 1
      }
    }
    x <- do.call("rbind",cleaned.segmentation)

    write_segs(x,"clean2")   # mosaicism:debugging

    ##print(x)
    # Now the hard part, fixing the boundaries
    # first we enlarge or shrink the regions where a CN change is detected
    if (nrow(x) > 1){
      row.ind <- 1
      if (x[1,3] != 0){
        end.adj <- GetBoundaryAdjustment(x[row.ind,3])
        x[row.ind,2] <- x[row.ind,2] - end.adj
      }
      row.ind <- 2
      while(row.ind < nrow(x)){
        if (x[row.ind,3] != 0){
          end.adj <- GetBoundaryAdjustment(x[row.ind,3])
          x[row.ind,1] <- max(1,x[row.ind,1] + end.adj)
          x[row.ind,2] <- x[row.ind,2] - end.adj
        }
        row.ind <- row.ind + 1
      }
      if (x[row.ind,3] != 0){
        end.adj <- GetBoundaryAdjustment(x[row.ind,3])
        x[row.ind,1] <- x[row.ind,1] + end.adj
      }
      ##print(x)
      ## now adjust out any gaps
      if (x[1,3] == 0){
        x[1,2] <- x[2,1]-1
      }

      row.ind <- 2
      while (row.ind < nrow(x)){
        if (x[row.ind,3] ==0){
          x[row.ind,1] <- x[row.ind-1,2]+1
          x[row.ind,2] <- x[row.ind+1,1]-1
        }
        row.ind <- row.ind + 1
      }
      if (x[row.ind,3] == 0){
        x[row.ind,1] <- x[row.ind-1,2]+1
      }
    }

    write_segs(x,"clean3")   # mosaicism:debugging

    #### This is the NEW bit. Corrects issue with regions that get grown
    #### It also checks for segments that go outside the allowable boundaries

    outsidebounds <- any(x[,1] < 1 | x[,1] > x[nrow(x),2] | x[,2] < 1 | x[,2] > x[nrow(x),2])

    #### Now it checks for errounous CN=0 regions, removes them and then joins overlapping segments
    Errounous <- x[,1] > x[,2] & x[,3] == 0.0
    if (any(Errounous) | outsidebounds){
      x[x[,1] < 1,1] <- 1
      x[x[,1] > x[nrow(x),2] ,1] <-  x[nrow(x),2]
      x[x[,2] < 1,2] <- 1
      x[x[,2] > x[nrow(x),2] ,2] <-  x[nrow(x),2]

      ## get rid of the incorrect segments
      x <- x[!Errounous,,drop=FALSE]
      ##print(x)
      x <- join.overlaps(x)
    }
    write_segs(x,"clean4")   # mosaicism:debugging
    x
  }

  ## This function produces an intial segmentation
  produce.segmentation <- function(x, gains.boundaries, losses.boundaries){

    segment.set <- NULL

    # the default values:
    # gains.boundaries <- c(0.08764945, 0.15380349, 0.21465931, 0.27100300)
    # losses.boundaries <- c(-0.08293345, -0.17551812, -0.28048196, -0.40165383)

    Neutral.zone.low <- losses.boundaries[1] # -0.08293345 (1-based)
    Neutral.zone.high <- gains.boundaries[1] #  0.08764945

    segment.start.ind <- 1
    if ((Neutral.zone.low < x[1]) & (x[1] <= Neutral.zone.high)){
      ## starting in a copy number neutral state
      i <- 1
    } else {
      ## starting in a non-neutral CN state
      cur.max <- x[1]
      i <- 1
      if (x[i] < Neutral.zone.low){
        ## It is a loss region
        cur.max <- x[i]
        while (x[i] < Neutral.zone.low & i < length(x)){
          i <- i + 1
          if (x[i] < cur.max){
            cur.max <- x[i]
          }
        }
     ##   cat(segment.start.ind,i-1,-get.mos.prop(cur.max,losses.boundaries),"\n")
        if (i < length(x)){
          segment.set <- c(segment.set,list(c(segment.start.ind,i-1,-get.mos.prop(cur.max,losses.boundaries))))
        } else {
          segment.set <- c(segment.set,list(c(segment.start.ind,i,-get.mos.prop(cur.max,losses.boundaries))))
        }
        segment.start.ind <- i
      } else {
        ## It is a gain
        segment.start.ind <- i
        cur.max <- x[i]
        while (x[i] > Neutral.zone.high & i < length(x)){
          i <- i + 1
          if (x[i] > cur.max){
            cur.max <- x[i]
          }
        }
    ##    cat(segment.start.ind,i-1,get.mos.prop(cur.max,gains.boundaries),"\n")
        if (i < length(x)){
          segment.set <- c(segment.set,list(c(segment.start.ind,i-1,get.mos.prop(cur.max,gains.boundaries))))
        } else {
          segment.set <- c(segment.set,list(c(segment.start.ind,i,get.mos.prop(cur.max,gains.boundaries))))
        }
        segment.start.ind <- i
      }
    }

    last_segment_start_ind <- segment.start.ind
    while (i < length(x)) {
      if (segment.start.ind != last_segment_start_ind) {
        cat("### PSeg: new seg: ",
            paste(segment.set[length(segment.set)],collapse="\t",sep="\t"),
            "\n",sep="")
        last_segment_start_ind <- segment.start.ind
      }

      cat("### PSeg:",i,x[i],segment.start.ind,"\n",sep="\t")

      if ((Neutral.zone.low < x[i]) & (x[i] <= Neutral.zone.high)){
        i <- i +1
      } else {
        segment.set <- c(segment.set,list(c(segment.start.ind,i-1,0)))
        segment.start.ind <- i
        if (x[i] < Neutral.zone.low){
          ## It is a loss region
          cur.max <- x[i]
          while (x[i] < Neutral.zone.low & i < length(x) ){
            i <- i + 1
            if (x[i] < cur.max){
              cur.max <- x[i]
            }
          }
       ##   cat(segment.start.ind,i-1,-get.mos.prop(cur.max,losses.boundaries),"\n")
          if (i < length(x)){
            segment.set <- c(segment.set,list(c(segment.start.ind,i-1,-get.mos.prop(cur.max,losses.boundaries))))
          } else {
            segment.set <- c(segment.set,list(c(segment.start.ind,i,-get.mos.prop(cur.max,losses.boundaries))))
          }
          segment.start.ind <- i
        } else {
          ## It is a gain
          segment.start.ind <- i
          cur.max <- x[i]
          while (x[i] > Neutral.zone.high & i < length(x) ){
            i <- i + 1
            if (x[i] > cur.max){
            cur.max <- x[i]
          }
          }
      ##    cat(segment.start.ind,i-1,get.mos.prop(cur.max,gains.boundaries),"\n")
          if (i < length(x)){
            segment.set <- c(segment.set,list(c(segment.start.ind,i-1,get.mos.prop(cur.max,gains.boundaries))))
          } else {
            segment.set <- c(segment.set,list(c(segment.start.ind,i,get.mos.prop(cur.max,gains.boundaries))))
          }
          segment.start.ind <- i
        }

      }
    }
    if (segment.start.ind < length(x)){
     ## cat(segment.start.ind,i,0,"\n")
      segment.set <- c(segment.set,list(c(segment.start.ind,i,0)))
    }
    segment.set <- do.call("rbind",segment.set)

    # return this
    segment.set
  }

  compute.confidences <- function(x,segment.set,gains.boundaries, losses.boundaries){

    confidences <- rep(0,nrow(segment.set))

    for (i in 1:nrow(segment.set)){
      cur.seg <- segment.set[i,]
      if (cur.seg[3] > 0){
        confidences[i] <- mean(runmed(x[cur.seg[1]:cur.seg[2]],251) > gains.boundaries[1])
      } else if (cur.seg[3] < 0){
        confidences[i] <- mean(runmed(x[cur.seg[1]:cur.seg[2]],251) < losses.boundaries[1])
      } else {
        confidences[i] <- mean(runmed(x[cur.seg[1]:cur.seg[2]],251) < losses.boundaries[1] | runmed(x[cur.seg[1]:cur.seg[2]],251) > gains.boundaries[1])
      }

    }

    cbind(segment.set,confidences)

  }



  #smooth out log ratios using long span moving average  x[i] = average(x[i-k/2], ..., x[i+ k/2])
  # Note that this "runmean" function handles endpoints (ie 1...k/2, n-k/2) as follows
  # for i=1...k2  x[i] = average(x[1], ... x[i+k/2])
  # for i=(n-k/2)...n  x[i] = average(x[i-k/2], ....,x[n])

  x.orig <- x
  x <- runmean(x,k=bandwidth)

  debug_write_segs_cnt <<- 0
  # mosaicism:debugging
  # means are ok.
  # filename <- paste(debug_current_basename,".runmean",sep="")
  # write.table(x,file=filename,row.names=FALSE,quote=FALSE,sep="\t")
  # cat("### wrote means:",filename,"\n")

  # Now get the segmentation for the mosaicism
  segment.set <- produce.segmentation(x,gains.boundaries, losses.boundaries)
  write_segs(segment.set,"after_prod")   # mosaicism:debugging

  segment.set <- clean.segments(segment.set,20)
  write_segs(segment.set,"after_clean")   # mosaicism:debugging

  segment.set <- compute.confidences(x.orig,segment.set,gains.boundaries, losses.boundaries)
  write_segs(segment.set,"after_conf")   # mosaicism:debugging

  segment.set
}
