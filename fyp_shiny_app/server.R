# -------------------------------------------------------------------
# The function that reads an individual simulation of a given mutant
get_sim <- function(AAname, simNum) {
    if (AAname %in% main_df$name) {
        fileName <- paste(AAname, ".csv", sep = "")
        path <- paste(data_directory, fileName, sep = "")
        
        df <- read.csv(path)
        
        if (simNum == 1) {
            return(df[c(10001:100000), ])
        } else if (simNum == 2) {
            return(df[c(110001:200000), ])
        } else if (simNum == 3) {
            return(df[c(209999:299998), ])
        }
    }
    else {
        return(0)
    }
}

# Function that gives all simulations concatenated
get_all_sim <- function(AAname, data_directory) {
    if (AAname %in% main_df$name) {
        fileName <- paste(AAname, ".csv", sep = "")
        path <- paste(data_directory, fileName, sep = "")
        df <- read.csv(path)
        
        
        return(rbind(df[c(c(10001:100000), c(110001:200000), c(209999:299998)), ]))
    }
    else {
        return(0)
    }
}


create_compare_average_distance_plot <- function(AAname) {
    df <- get_all_sim(AAname, "../data/extracted_data/")
    first_letter_distance <-
        read.csv(paste(
            paste(
                "../data/population_samples/",
                substr(AAname, 1, 1),
                sep = ""
            ),
            "_first.csv",
            sep = ""
        ))
    
    second_letter_distance <-
        read.csv(paste(
            paste(
                "../data/population_samples/",
                substr(AAname, 2, 2),
                sep = ""
            ),
            "_second.csv",
            sep = ""
        ))
    p <-
        ggplot(
            data = melt(
                list(
                    AAname = df,
                    first_AA_average = first_letter_distance,
                    second_AA_average = second_letter_distance
                ),
                id.vars = "time",
                measure.vars = "distance",
                value.name = "distance"
            ),
            aes(x = distance, color = L1)
        ) +
        labs(title = "Cavity Size Comparison",
             subtitle = str_wrap("Comparison of cavity size between the mutant and averages for each of its amino acid",60)) +
        geom_density(size = 1.5, alpha = 0.6) +
        scale_x_continuous(
            name = "Distance",
            limits = c(0.6, 2.0),
            breaks = seq(0.6, 2.0, 0.2)
        ) +
        scale_y_continuous(name = "Density",
                           breaks = seq(0, 10, 2)) +
        scale_colour_discrete(name = "",
                              labels = c(
                                  paste(AAname, "mutant", sep = " "),
                                  "110th Average",
                                  "129th Average"
                              )) +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.title.x = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.title.y = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.text.x = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.subtitle = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.title = element_text(
                size = 16,
                face = "bold",
                colour = "black"
            )
        )
    
    return(p)
}

create_compare_average_benzene_plot <-
    function(AAname, threshold = 8) {
        first_letter_benzene <-
            read.csv(paste(
                paste(
                    "../data/population_samples/",
                    substr(AAname, 1, 1),
                    sep = ""
                ),
                "_first_benzene.csv",
                sep = ""
            ))
        second_letter_benzene <-
            read.csv(paste(
                paste(
                    "../data/population_samples/",
                    substr(AAname, 2, 2),
                    sep = ""
                ),
                "_second_benzene.csv",
                sep = ""
            ))
        df_benzene <- load_benzene_data(AAname, threshold)
        p <- ggplot(
            data = melt(
                list(
                    AAname = df_benzene,
                    first_AA_average = first_letter_benzene,
                    second_AA_average = second_letter_benzene
                ),
                id.vars = "Time",
                measure.vars = "Distance",
                value.name = "distance"
            ),
            aes(x = distance, color = L1)
        ) +
            scale_y_continuous(name = "Density") +
            labs(title = "Benzene binding Comparison",
                 subtitle = str_wrap("Comparison of benzene distance from the cavity between the mutant and averages for each of its amino acid", 60)) +
            geom_density(size = 1.5, alpha = 0.6) +
            scale_x_continuous(
                name = "Distance",
                limits = c(0, 30),
                breaks = seq(0, 30, 5)
            ) +
            theme_minimal() +
            theme(
                axis.text.y = element_blank(),
                axis.title.x = element_text(
                    size = 14,
                    face = "bold",
                    colour = "black"
                ),
                axis.title.y = element_text(
                    size = 14,
                    face = "bold",
                    colour = "black"
                ),
                axis.text.x = element_text(
                    size = 12,
                    face = "bold",
                    colour = "black"
                ),
                plot.subtitle = element_text(
                    size = 12,
                    face = "bold",
                    colour = "black"
                ),
                plot.title = element_text(
                    size = 16,
                    face = "bold",
                    colour = "black"
                )
            )
        
        return(p)
    }

create_collage <- function(AAname, threshold = 8) {
    p1 <- create_compare_average_distance_plot(AAname)
    
    p2 <- create_compare_average_benzene_plot(AAname, threshold)
    
    ggarrange(
        p1,
        p2,
        ncol = 2,
        nrow = 1,
        common.legend = TRUE,
        legend = "bottom"
    )
}

load_benzene_data <- function(aa, threshold) {
    directory_name <-
        paste(
            paste(
                "../data/bnz_distance/BNZ-COM-BNZ-AM-1-W15-110_129-",
                aa,
                sep = ""
            ),
            ".csv",
            sep = ""
        )
    df <- read.csv(directory_name)
    names(df) <- c("Time", "Distance")
    df$isInside <- df$Distance < threshold
    return(df)
}

benzene_timeseries_plot <- function(df, threshold) {
    ggplot(data = df, aes(x = Time, y = Distance)) +
        geom_line(color = "darkgrey") +
        geom_hline(
            yintercept = threshold,
            linetype = "dashed",
            size = 1.2,
            color = "red"
        ) +
        theme_minimal() +
        theme(
            axis.title.x = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.title.y = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.text.x = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.subtitle = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.title = element_text(
                size = 16,
                face = "bold",
                colour = "black"
            )
        )
}

distance_density_plot <- function(df, aminoAcids) {
    df$simulation <- as.factor(df$simulation)
    p <- ggplot(data = df, aes(x = distance, color = simulation)) +
        geom_density(size = 1.5, alpha = 0.6) +
        scale_x_continuous(
            name = "Distance",
            limits = c(0.6, 2.0),
            breaks = seq(0.6, 2.0, 0.2)
        ) +
        scale_y_continuous(name = "Density") +
        theme_minimal() +
        theme(
            axis.text.y = element_blank(),
            axis.title.x = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.title.y = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            ),
            axis.text.x = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.subtitle = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            plot.title = element_text(
                size = 16,
                face = "bold",
                colour = "black"
            ),
            legend.text = element_text(
                size = 12,
                face = "bold",
                colour = "black"
            ),
            legend.title = element_text(
                size = 14,
                face = "bold",
                colour = "black"
            )
        )
    return(p)
}

identify_colour <- function(sd, mean, value) {
    difference <- abs(mean - value)
    if (difference < sd) {
        return("green")
    } else if (difference < 1.5 * sd) {
        return("yellow")
    } else if (difference < 2 * sd) {
        return("orange")
    } else{
        return("red")
        
    }
}

# Define a server for the Shiny app
function(input, output) {
    
    ranges <- reactiveValues(x = NULL, y = NULL)
    # Fill in the spot we created for a plot
    output$AAplot <- renderPlot({
        amino_acid <- paste(input$firstAA, input$secondAA, sep = "")
        p1 <-
            distance_density_plot(get_all_sim(amino_acid, "../data/extracted_data/"),
                                  amino_acid)
        p1 +
            coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    })
    
    output$benzene_plot <- renderPlot({
        amino_acid <- paste(input$firstAA, input$secondAA, sep = "")
        df <- load_benzene_data(amino_acid, input$benzene_limit)
        insidePercentage <-
            length(df$isInside[df$isInside == TRUE]) / 5000
        output$insidePercentage <- renderText({
            
        })
        output$insidePercentage_valuebox <- renderValueBox({
            valueBox(
                paste(round(insidePercentage, digits = 2) * 100, "%", sep = ""),
                "Benzene bind",
                color = identify_colour(
                    sd(main_df$isInsidePercent),
                    mean(main_df$isInsidePercent),
                    insidePercentage
                )
            )
        })
        benzene_timeseries_plot(df, input$benzene_limit)
    })
    
    output$overview_title <-
        renderText({
            paste(paste(input$firstAA, input$secondAA, sep = ""),
                  " Profile",
                  sep = " ")
        })
    
    output$pi_valuebox <- renderValueBox({
        piValue <-
            main_df$pi[main_df$name == paste(input$firstAA, input$secondAA, sep = "")]
        
        valueBox(round(piValue, digits = 2),
                 "pi value",
                 color = identify_colour(sd(main_df$pi), mean(main_df$pi), piValue))
    })
    
    output$hydrophobicity_valuebox <- renderValueBox({
        hydrophobicity <-
            main_df$hydrophobicity[main_df$name == paste(input$firstAA, input$secondAA, sep = "")]
        
        valueBox(
            round(hydrophobicity, digits = 2),
            "hydrophobicity",
            color = identify_colour(
                sd(main_df$hydrophobicity),
                mean(main_df$hydrophobicity),
                hydrophobicity
            )
        )
    })
    
    output$mass_valuebox <- renderValueBox({
        mass <-
            main_df$mass[main_df$name == paste(input$firstAA, input$secondAA, sep = "")]
        
        valueBox(value = round(mass, digits = 2), subtitle = "mass",
                 color = identify_colour(sd(main_df$mass), mean(main_df$mass), mass))
    })
    
    output$var_score_valuebox <- renderValueBox({
        var_score <-
            main_df$var_score[main_df$name == paste(input$firstAA, input$secondAA, sep = "")]
        
        valueBox(round(var_score, digits = 2),
                 "var_score",
                 color = identify_colour(
                     sd(main_df$var_score),
                     mean(main_df$var_score),
                     var_score
                 ))
    })
    
    output$vdw_valuebox <- renderValueBox({
        vdwValue <-
            main_df$vdw_volume[main_df$name == paste(input$firstAA, input$secondAA, sep = "")]
        
        valueBox(round(vdwValue, digits = 2),
                 "Van der Waals volume",
                 color = identify_colour(sd(main_df$vdw_volume), mean(main_df$vdw_volume), vdwValue))
    })
    
    
    output$message <- renderText({shiny::tags$p("some text", style = "color:red; font-size: 200%;")})
    
    output$collage <-
        renderPlot({
            create_collage(AAname = paste(input$firstAA, input$secondAA, sep = ""))
        })
    
    
    # When a double-click happens, check if there's a brush on the plot.
    # If so, zoom to the brush bounds; if not, reset the zoom.
    observeEvent(input$plot1_dblclick, {
        brush <- input$plot1_brush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
            ranges$y <- c(brush$ymin, brush$ymax)
            
        } else {
            ranges$x <- NULL
            ranges$y <- NULL
        }
    })
    
    output$contents <- renderTable({
        
        # input$file1 will be NULL initially. After the user selects
        # and uploads a file, head of that data file by default,
        # or all rows if selected, will be shown.
        
        req(input$file1)
        
        # when reading semicolon separated files,
        # having a comma separator causes `read.csv` to error
        tryCatch(
            {
                df <- read.csv(input$file1$datapath,
                               header = input$header,
                               sep = input$sep,
                               quote = input$quote)
            },
            error = function(e) {
                # return a safeError if a parsing error occurs
                stop(safeError(e))
            }
        )
        
        if(input$disp == "head") {
            return(head(df))
        }
        else {
            return(df)
        }
        
    })
    
    runjs('
        var el2 = document.querySelector(".skin-blue");
        el2.className = "skin-blue sidebar-mini";
        ')
}