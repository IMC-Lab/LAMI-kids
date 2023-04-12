## This file contains utilities for creating
## spatial grids of rectangles and hexagons
library(tidyverse)

## Given hexagon centers (bin_x, bin_y) and distances
## between adjacent hexagons (xsize, ysize), returns
## a dataframe with all vertices of each hexagon
vertices.rect <- function(grid_idx, bin_x, bin_y, xsize, ysize) {
    if (length(bin_x) != length(bin_y)) stop('x and y must be the same length')
    
    data.frame(i=1:4) %>%
        mutate(mx=ifelse(i<=2, 1, -1),
               my=ifelse(i==1|i==4, 1, -1)) %>%
        group_by(i, mx, my) %>%
        expand(nesting(grid_idx=grid_idx, bin_x=bin_x, bin_y=bin_y)) %>%
        mutate(vx=bin_x + mx*xsize/2,
               vy=bin_y + my*ysize/2) %>%
        group_by(grid_idx, bin_x, bin_y) %>%
        summarize(vi=list(i), vx=list(vx), vy=list(vy)) %>%
        return()
}

## Create a grid of spatial locations determined by
## the minimum, maximum, and step size for the x and y axes.
## This grid will be replicated for all existing grouping variables
## in the data frame df.
grid.rect <- function(df=data.frame(), xmin=-600, xmax=600, xsize=10,
                      ymin=-300, ymax=300, ysize=10) {
    ## make sure we include the maxima
    if ((xmax-xmin) %% xsize != 0) xmax <- xmax + xsize
    if ((ymax-ymin) %% ysize != 0) ymax <- ymax + ysize
    
    df %>%
        expand(bin_x=seq(xmin, xmax, xsize),
               bin_y=seq(ymin, ymax, ysize)) %>%
        mutate(grid_idx=row_number()) %>%
        left_join(., vertices.rect(.$grid_idx, .$bin_x, .$bin_y, xsize, ysize)) %>%
        return()
}

## Create a 2D histogram of the data binned over grid.rect(...)
bin.rect <- function(df, xmin=-600, xmax=600, xsize=10,
                     ymin=-300, ymax=300, ysize=10) {
    grid <- grid.rect(df, xmin=xmin, xmax=xmax, xsize=xsize,
                      ymin=ymin, ymax=ymax, ysize=ysize)
    
    ## find break points from bin centers
    x_bins <- unique(grid$bin_x)
    y_bins <- unique(grid$bin_y)
    x_breaks <- x_bins[-length(x_bins)] + xsize/2
    y_breaks <- y_bins[-length(y_bins)] + ysize/2
    
    ## bin df according to the grid
    df %>%
        mutate(bin_x=x_bins[findInterval(fix_x, x_breaks) + 1],
               bin_y=y_bins[findInterval(fix_y, y_breaks) + 1]) %>%
        group_by(bin_x, bin_y, .add=TRUE) %>%
        summarize(count=n()) %>%
        right_join(grid) %>%
        mutate(count=replace_na(count, 0)) %>%
        return()
}

## Given hexagon centers (bin_x, bin_y) and distances
## between adjacent hexagons (xsize, ysize), returns
## a dataframe with all vertices of each hexagon
vertices.hex <- function(grid_idx, bin_x, bin_y, xsize, ysize) {
    if (length(bin_x) != length(bin_y)) stop('x and y must be the same length')
    
    data.frame(i=1:6) %>%
        mutate(angle_deg=60*i - 30,
               angle_rad=pi/180*angle_deg) %>%
        group_by(i, angle_deg, angle_rad) %>%
        expand(nesting(grid_idx=grid_idx, bin_x=bin_x, bin_y=bin_y)) %>%
        mutate(vx=bin_x+xsize/sqrt(3)*cos(angle_rad),
               vy=bin_y+ysize*2/3*sin(angle_rad)) %>%
        group_by(grid_idx, bin_x, bin_y) %>%
        summarize(vi=list(i), vx=list(vx), vy=list(vy)) %>%
        return()
}

## Create a hexagonal grid over the domain (xmin, xmax) and
## range (ymin, ymax). xsize and ysize are the spacing between
## adjacent hexagons in each row/column.
##
## By default, ysize is set so that all hexagons are regular
grid.hex <- function(df=data.frame(), xmin=-600, xmax=600, xsize=10,
                     ymin=-300, ymax=300, ysize=xsize*sqrt(3)/2,
                     rounded=FALSE, shift_idx=c(FALSE, TRUE)) {
    grid <- grid.rect(df=df, xmin=xmin, xmax=xmax, xsize=xsize,
                      ymin=ymin, ymax=ymax, ysize=ysize) %>%
        select(-vi, -vx, -vy)
    shift_rows <- unique(grid$bin_y)[shift_idx]

    ## add in some "missing" hexes on the left side
    if (rounded) {
        grid <- grid %>%
            bind_rows(expand(df, bin_x=xmin-xsize,
                             bin_y=shift_rows) %>%
                      mutate(grid_idx=row_number()+max(grid$grid_idx))) %>%
            arrange(bin_x, bin_y) %>%
            mutate(grid_idx=row_number())
    }
    
    ## shift every other row and calculate vertices
    grid %>%
        mutate(bin_x=bin_x + xsize/2*(bin_y %in% shift_rows)) %>%
        left_join(., vertices.hex(.$grid_idx, .$bin_x, .$bin_y, xsize, ysize)) %>%
        return()
}

## Calculate the squared euclidean distance between two or more points
dist_sq <- function(x1, y1, x2, y2) {
    return((x1-x2)**2 + (y1-y2)**2)
}

## Create a 2D histogram of df over the grid grid.hex(df, ...)
bin.hex <- function(df, ...) {
    grid <- grid.hex(df, ...)
    g <- grid.hex(ungroup(df), ...)
    
    df %>%
        rowwise() %>%
        mutate(grid_idx=g$grid_idx[which.min(dist_sq(fix_x, fix_y,
                                                     g$bin_x, g$bin_y))]) %>%
        group_by(grid_idx, .add=TRUE) %>%
        summarize(count=n()) %>%
        right_join(grid) %>%
        mutate(count=replace_na(count, 0)) %>%
        return()
}
