tags <- shiny::tags

header <- dashboardHeader()

anchor <- tags$a(
    tags$img(
        src = 'logo.png',
        height = '35',
        width = '35'
    ),
    'Mutant Explorer',
    style = "color:white"
)

header$children[[2]]$children <- tags$div(anchor,
                                          class = 'name')

dashboardPage(
    skin = "blue",
    header,
    ## Sidebar content
    dashboardSidebar(
        disable = TRUE,
        collapsed = TRUE,
        sidebarMenu(
            menuItem(
                "Mutant Profile",
                tabName = "mutant_profile",
                icon = icon("dashboard")
            ),
            menuItem("Compare", tabName = "compare", icon = icon("eye")),
            menuItem("Prediction", tabName = "prediction", icon = icon("road")),
            menuItem("About", tabName = "about", icon = icon("info"))
            
        )
    ),
    
    ## Body content
    dashboardBody(
        #useShinyjs(),
        tabItems(
            # First tab content
            tabItem(
                tabName = "mutant_profile",
                fluidRow(
                    box(
                        title = "Controls",
                        status = "info",
                        solidHeader = TRUE,
                        width = 3,
                        height = 300,
                        fluidRow(column(
                            12,
                            selectInput("firstAA", "110 Amino Acid:",
                                        choices =
                                            aa_list)
                        ), ),
                        fluidRow(column(
                            12,
                            selectInput("secondAA", "129 Amino Acid:",
                                        choices =
                                            aa_list)
                        ), ),
                        fluidRow(column(
                            12,
                            sliderInput(
                                "benzene_limit",
                                label = "Benzene binding threshold:",
                                min = 5,
                                max = 30,
                                value = 10,
                            )
                        ), ),
                        
                    ),
                    box(
                        title = "Overview",
                        status = "primary",
                        solidHeader = TRUE,
                        width = 9,
                        height = 300,
                        fluidRow(
                            valueBoxOutput("pi_valuebox"),
                            valueBoxOutput("hydrophobicity_valuebox"),
                            valueBoxOutput("mass_valuebox"),
                            valueBoxOutput("var_score_valuebox"),
                            valueBoxOutput("insidePercentage_valuebox"),
                            valueBoxOutput("vdw_valuebox")
                        )
                    )
                    
                ),
                fluidRow(
                    box(
                        title = "Distance Plot",
                        status = "primary",
                        solidHeader = TRUE,
                        width = 6,
                        plotOutput(
                            "AAplot",
                            height = 400,
                            dblclick = "plot1_dblclick",
                            brush = brushOpts(id = "plot1_brush",
                                              resetOnNew = TRUE)
                        )
                    ),
                    box(
                        title = "Benzene Plot",
                        status = "primary",
                        solidHeader = TRUE,
                        width = 6,
                        plotOutput("benzene_plot", height = 400)
                    ),
                    
                ),
                fluidRow(
                    box(
                        title = "Compare to Average",
                        status = "primary",
                        solidHeader = TRUE,
                        width = 12,
                        plotOutput("collage", height = 400)
                    ),
                )
            ),
            
            # Second tab content
            tabItem(tabName = "compare",
                    fluidRow(
                        box(
                            title = "First Mutant Options",
                            status = "info",
                            solidHeader = TRUE,
                            width = 3,
                            height = 300,
                            fluidRow(column(
                                12,
                                selectInput("compare_firstAA_1", "110 Amino Acid:",
                                            choices =
                                                aa_list)
                            ), ),
                            fluidRow(column(
                                12,
                                selectInput("compare_secondAA_1", "129 Amino Acid:",
                                            choices =
                                                aa_list)
                            ), ),
                            fluidRow(column(
                                12,
                                sliderInput(
                                    "compare_benzene_threshold_1",
                                    label = "Benzene binding threshold:",
                                    min = 5,
                                    max = 30,
                                    value = 10,
                                )
                            ), ),
                            
                        ),
                        box(
                            title = "Overview",
                            status = "primary",
                            solidHeader = TRUE,
                            width = 6,
                            height = 300,
                            
                        ),
                        box(
                            title = "Second Mutant Options",
                            status = "info",
                            solidHeader = TRUE,
                            width = 3,
                            height = 300,
                            fluidRow(column(
                                12,
                                selectInput("compare_firstAA_2", "110 Amino Acid:",
                                            choices =
                                                aa_list)
                            ), ),
                            fluidRow(column(
                                12,
                                selectInput("compare_secondAA_2", "129 Amino Acid:",
                                            choices =
                                                aa_list)
                            ), ),
                            fluidRow(column(
                                12,
                                sliderInput(
                                    "compare_benzene_threshold_2",
                                    label = "Benzene binding threshold:",
                                    min = 5,
                                    max = 30,
                                    value = 10,
                                )
                            ), ),
                            
                        )
                    )),
            
            # Second tab content
            tabItem(tabName = "prediction",
                    h2("Prediction")),
            tabItem(tabName = "about",
                    h2("About"))
        )
    )
)
