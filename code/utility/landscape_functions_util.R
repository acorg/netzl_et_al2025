padding <- 0.5
opacity_val <- 0.9
# Functions to remove buttons
addObject3js <- function(
    data3js,
    object,
    number_of_ids = 1
){
  
  # Generate an object ID
  if(is.null(data3js$lastID)){ data3js$lastID <- 0 }
  object$ID <- max(data3js$lastID) + seq_len(number_of_ids)
  
  # If object is interactive and highlighted add a reference to itself to
  # it's highlight group by default
  if(!is.null(object$properties$interactive)){
    object$group <- object$ID
  }
  
  # Add the object to the plot data
  data3js$plot[[length(data3js$plot)+1]] <- object
  
  # Update the ID of the last object added
  data3js$lastID <- object$ID
  
  # Return the new data
  data3js
  
}

remove_buttons <- function(data3js){
  
  new_data3js = data3js
  
  new_data3js = data3js
  
  new_data3js[['lastID']] = 0
  new_data3js[['plot']] = list()
  
  N = data3js[['lastID']] 
  
  
  
  
  for (i in 1:N)
  {
    obj = data3js[['plot']][[i]]
    
    
    
    if ('toggle' %in% names(obj[['properties']])){
      obj[['properties']][['toggle']] <- NULL
    }
    
    new_data3js = addObject3js(new_data3js,obj)
    
    
  }
  
  
  
  return (new_data3js)
  
}


base_plot_data3js <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                              add_border = TRUE, add_axis = TRUE){
  
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  z_coords <- rep(0.02, length(highlighted_ags))
  ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
  ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
  ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
  labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  border_col <- "grey50"
    
  z_lims <- c(0,10)
  axis_at <- seq(z_lims[1], z_lims[2],2)
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "log"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){
    
    axis_labels <- 2^axis_at*10
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      # labeloffset = 0.11,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}

plot_idvl_landscapes_from_list <- function(data3js, idvl_landscapes, sr_colors){
  
  for (x in seq_along(idvl_landscapes)) {
    
    surface_options <- list()
    surface_options$col.surface = sr_colors[x]
    surface_options$col.surface.grid = adjustcolor(
      "grey",
            red.f = 0.25,
            green.f = 0.25,
            blue.f = 0.25
    )
    surface_options$opacity.surface = 0.2
    
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = idvl_landscapes[[x]],
      toggle = x,
      options = surface_options,
      crop2chull = FALSE,
      grid_spacing = 0.5,
      padding = padding
    )
    
  }
  
  return(data3js)
}


plot_landscapes_from_list <- function(data3js, titertables_groups, lndscp_fits,map, gmt_data, highlighted_ags,
                                      ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1, hide_buttons = TRUE, add_ag_label = FALSE, lndscp_colors,
                                      show_gmts = TRUE,
                                      show.individual.surfaces = FALSE,
                                      options.individual.surfaces = list(opacity.surface.grid = 0.4, 
                                                                         opacity.surface = 0.2, col.surface = "grey70", col.surface.grid = "grey70")){
  
  if(alternative_ba5){
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2], agCoords(map, optimization_number = opti_nr)[agNames(map) %in% "BA.4/BA.5", 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags)), 12) / 5
    # text_x <- c(c(x_coords[1:4] + ag_point_size[1:4]*0.15),c(x_coords[5:6] - ag_point_size[5:6]*0.25))
    #  text_y <- c(y_coords[1:4], c(y_coords[5:6] - ag_point_size[5:6]*0.2))
    text_x <- c(x_coords[1:6] - ag_point_size[1:6]*0.2)
    text_y <- c(y_coords[1:6] - ag_point_size[1:6]*0.2)
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]], "BA.4/BA.5(2)")
    
  } else {
    x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
    y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
    z_coords <- rep(0.02, length(highlighted_ags))
    ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
    #  text_x <- c(c(x_coords[1:3] + ag_point_size[1:3]*0.15),c(x_coords[4:5] - ag_point_size[4:5]*0.2))
    #  text_y <- c(y_coords[1:3], c(y_coords[4:5] - ag_point_size[4:5]*0.2))
    text_x <- c(x_coords[1:5] - ag_point_size[1:5]*0.2)
    text_y <- c(y_coords[1:5] - ag_point_size[1:5]*0.2)
  }
  
  if(add_ag_label){
    text_plot <- c(ag_plot_names[agNames(map)[agNames(map) %in% highlighted_ags]])
  } else {
    text_plot <- rep("", length(highlighted_ags))
  }
  
  if(length(lndscp_fits) > 1){
    min_offset <- -0.1
    max_offset <- 0.1
  } else {
    min_offset <- 0
    max_offset <- 0
  }
  offset <- seq(from = min_offset, to = max_offset, by = (max_offset - min_offset)/length(lndscp_fits))
  
  
  for (i in seq_along(lndscp_fits)) {
    
    # message(i)
    srg <- titertables_groups$sr_group[i]
    lndscp_fit <- lndscp_fits[[i]]
    
    coords <- cbind(x_coords, y_coords)
    
    coords <- coords[!is.na(x_coords),]
    # Add titers
    gmts <- filter(gmt_data, sr_group == srg)
    
    gmts <- gmts[match(rownames(coords), gmts$ag_name),]
    
    for (j in seq_len(nrow(coords))) {
      
      if(show_gmts){
        
        data3js <- r3js::lines3js(
          data3js,
          x = rep(coords[j, 1], 2),
          y = rep(coords[j, 2], 2),
          z = c(0, gmts$gmt[j]),
          col = "grey50",
          highlight = list(col = "red"),
          interactive = FALSE,
          #   toggle = sprintf("Titers, %s, %s", arm, visit),
          geometry = TRUE,
          opacity = 0.7,
          lwd = 0.2 #was 0.4
        )
        
        data3js <- r3js::points3js(
          data3js,
          x         = coords[j, 1],# + offset[i],
          y         = coords[j, 2],
          z         = gmts$gmt[j],
          #  size      = 0.7, #was 2
          size      = 0.9, #was 0.9
          #  col       = "grey50",
          col  = lndscp_colors[as.character(srg), 'Color'],
          highlight = list(col = "red"),
          #     label     = gmts$variant[j],
          #   toggle    = sprintf("Titers, %s, %s", arm, visit),
          opacity   = 1 # was 1
          
        )
      }
      # text_x <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 1] + agSize(map)[agNames(map) %in% highlighted_ags[1:4]]*0.025), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 1]- agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
      #  text_y <- c(c(agCoords(map)[agNames(map) %in% highlighted_ags[1:4], 2]), c(agCoords(map)[agNames(map) == "BA.4/BA.5", 2] - agSize(map)[agNames(map) == "BA.4/BA.5"]*0.07))
      
      # set points and coordinates of highlighted ags
      
      
      data3js <- r3js::text3js(
        data3js,
        x          = text_x,
        y          = text_y,
        z          = z_coords,
        text       = text_plot,
        toggle     = "Labels",
        size       = c(rep(14*0.02, length(text_x))), #agSize(map)[agNames(map) %in% highlighted_ags]*0.02,
        alignment  = "right"
      )
      
    }
    
    
    # Add landscapes
    data3js <- lndscp3d_surface(
      data3js = data3js,
      object = lndscp_fit,
      # zlim = c(0, 10),
      crop2chull = FALSE,
      # crop2base = TRUE,
      toggle = sprintf("Landscape, %s", srg),
      grid_spacing = 0.5,
      padding = padding,
      options = list(
        col.surface = lndscp_colors[as.character(srg), 'Color'],
        # opacity.surface = 0.5
        opacity.surface = opacity_val
      )
    )
    
    fit <- lndscp_fit
    if (show.individual.surfaces) {
      for (i in seq_len(nrow(fit$titers))) {
        individual_fit <- fit
        individual_fit$titers <- fit$titers[i, ]
        individual_fit$logtiters <- fit$logtiters[i, ]
        individual_fit$logtiters.upper <- individual_fit$logtiters.upper[i, 
        ]
        individual_fit$logtiters.lower <- individual_fit$logtiters.lower[i, 
        ]
        individual_fit$lessthans <- individual_fit$lessthans[i, 
        ]
        individual_fit$morethans <- individual_fit$morethans[i, 
        ]
        individual_fit$fitted.values <- NULL
        individual_fit$residuals <- NULL
        individual_fit$residuals.lessthan <- NULL
        individual_fit$residuals.morethan <- NULL
        if (!is.null(individual_fit$cone)) {
          individual_fit$cone$cone_coords <- individual_fit$cone$cone_coords[i, 
                                                                             , drop = F]
          individual_fit$cone$cone_heights <- individual_fit$cone$cone_heights[i]
        }
        data3js <- lndscp3d_surface(data3js = data3js, object = individual_fit, 
                                    crop2chull = FALSE, grid_spacing = 0.5, 
                                    padding = padding,
                                    options = options.individual.surfaces, toggle = "Individual surfaces")
      }
    }
    
    
  }
  
  if(hide_buttons){
    data3js <- remove_buttons(data3js)
  }
  
  
  
  return(data3js)
}



# sams landscape functions to add landscape from lndscp fits list
get_titertable <- function(data, group) {
  
  
  mean_titer <- function(titers) {
    if(length(titers) == 0) {
      "*"
    } else if(length(titers) >1) {
      mean_t <- mean(log2(as.numeric(titers)/10), na.rm = TRUE)
      as.character(2^mean_t*10)
    } else {
      titers
    }
  }
  
  data %>% 
    select(
      ag_name,
      sr_name,
      titer
    ) %>%
    mutate(
      titer = replace(titer, is.na(titer), "*")
    ) %>%
    pivot_wider(
      id_cols = sr_name,
      names_from = ag_name,
      values_from = titer,
      values_fn = mean_titer
    ) %>% 
    as.matrix() -> titermatrix
  
  attr(titermatrix, "sr_group") <- group$sr_group
  rownames(titermatrix) <- titermatrix[,"sr_name"]
  titermatrix <- titermatrix[,-1]
  
  #  print(titermatrix)
  #  print(titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F])
  #  titermatrix[titermatrix[,"BA.4/BA.5"] != "*",,drop=F]
  #  titermatrix[titermatrix[,"BA.4/BA.5(2)"] != "*",,drop=F]
  
  return(titermatrix)
  
}


plot_single_landscape_panel_webshot <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                                sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE,
                                                delete_html = TRUE, save_name = "temp"){
  
  
  to_save <- file.path(paste0(save_name, ".html"))
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  webshot(url=to_save,file = png_save)
  temp_plot <- readPNG(png_save)
  
  qplot(c(1:10),c(1:10), geom="blank") +
    annotation_custom(rasterGrob(temp_plot, height = unit(0.7, "npc")), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
    annotate(geom="text", x=label_x_pos, y=label_y_pos, label=label,size= label_size, hjust = 0) + 
    annotate(geom="text", x=label_x_pos, y=sr_group_y_pos, label=sr_group_label,size= sr_group_size, hjust = 0) +
    theme_void() -> p
  
  if(show_border) {
    p + theme(panel.border = element_rect(color = "grey50",
                                          fill = NA,
                                          size = 0.5))-> p
  }
  
  if(delete_html){
    if (file.exists(to_save)) {
      #Delete file if it exists
      file.remove(to_save)
    }
    if (file.exists(png_save)) {
      #Delete file if it exists
      file.remove(png_save)
    }
  }
  
  return(p) 
}


plot_single_landscape_panel <- function(landscape, label, label_size = 10, label_x_pos = 2, label_y_pos = 9,
                                        sr_group_label = "", sr_group_y_pos = 0, sr_group_size = 3, show_border = FALSE,
                                        delete_html = TRUE, save_name = "temp"){
  
  
  to_save <- file.path(paste0(save_name, ".html"))
  png_save <- gsub(".html", ".png", to_save)
  saveWidget(landscape, to_save, selfcontained = FALSE)
  
}



combine_landscape_and_calculated_gmt <- function(lndscp_fits, gmt_data, sr_group_fields = 1){
  
  lndscp_gmts <- lapply(lndscp_fits, function(x){
    
    gmts <- data.frame(logtiter = x$fitted.value,
                       ag_name = names(x$fitted.value),
                       sr_group = paste0(str_split(rownames(x$logtiters)[1], "_")[[1]][sr_group_fields], collapse = "_"))
    return(gmts)
  })
  
  lndscp_gmts <- do.call(rbind, lndscp_gmts)
  
  
  ## copmare lndscp gmts and claculated gmts
  comb_gmt <- rbind(lndscp_gmts %>%
                      mutate(Data = "Fitted Landscape GMT"),
                    gmt_data %>%
                      mutate(logtiter = gmt) %>%
                      select(sr_group, ag_name, logtiter) %>%
                      unique() %>%
                      mutate(Data = "Calculated GMT"))
  
  
  return(comb_gmt)
}



plot_lndscp_calculated_gmt_lineplot <- function(comb_data, ag_order = c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2", "BA.1", "BA.2"), lower_lim = -1, plot_colors = NULL){
  
  comb_data %>%
    ggplot(aes(x = ag_name, y = logtiter, color = Data, fill = Data)) + 
    geom_line(aes(group = Data), position = position_dodge(width = 0.3)) + 
    geom_point(shape = 21, color = "grey20", position = position_dodge(width = 0.3)) + 
    # scale_fill_manual(values = plot_colors,
    #                   name = "Serum group") +
    # scale_color_manual(values = plot_colors,
    #                    name = "Serum group") +
    scale_x_discrete(name = "Variant",
                     limits = ag_order) + 
    scale_y_continuous(limits = c(lower_lim, NA),
                       labels = function(x) round(2^x*10,2), 
                       breaks = c(-3:10),
                       name = "GMT") + 
    facet_wrap(~sr_group,
               labeller = label_wrap_gen(multi_line = TRUE),
               ncol = 6) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) -> p
  
  if(!is.null(plot_colors)){
    p <- p +
      scale_color_manual(values = plot_colors) + 
      scale_fill_manual(values = plot_colors)
  }
  
  return(p)
  
}

plot_gmt_diff <- function(gmt_diff, ag_order = c("D614G", "B.1.1.7", "B.1.351", "B.1.617.2", "BA.1", "BA.2")){
  
  gmt_diff <- gmt_diff %>%
    filter(Data != "Calculated GMT")
  
  gmt_diff %>%
    ggplot(aes(x = ag_name, y = gmt_diff, group = sr_group)) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point() +
    scale_x_discrete(name = "Variant",
                       limits = ag_order) +
    scale_y_continuous(name = "GMT difference",
                     #  limits = c(-2, 2),
                       breaks = seq(-6, 2, 0.5),
                       labels = function(x) ifelse(x < 0,paste0("-", round(2^abs(x), 1), "x"), paste0(round(2^x, 1), "x"))) +
    facet_wrap(~sr_group, ncol = 6) + 
    theme_bw() +
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          axis.text.y = element_text(size = 7)) -> p_gmt_diff
  
  return(p_gmt_diff)
  
}



