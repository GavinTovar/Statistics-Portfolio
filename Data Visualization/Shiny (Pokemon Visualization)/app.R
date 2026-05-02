### (ST 537) Animation/Interaction Self Directed Study
# setwd("C:/Users/Ghcto/OneDrive/Desktop/Spring 2024/ST 537/Module 9/Assignment 4")
# setwd("C:/Users/Ghcto/OneDrive/Desktop/School/Spring 2024/ST 537/Module 9")
Pokemon <- read.csv("Pokemon.csv")
library(shiny)
library(bslib)
library(tidyverse)
library(colorspace)

### Initial Cleaning
# Say's "Blastoise" is a type of Pokemon and can correct it
Pokemon[which(Pokemon$type1 == "Blastoise"),]
Pokemon$type1[16] <- "Water"
Pokemon$type2[16] <- ""

# Noticed "grass" misspelled for one observation
Pokemon[which(Pokemon$type1 == "Graass"),]
Pokemon$type1[979] <- "Grass"

# No Generation 0 wanted
Pokemon <- Pokemon %>%
  filter(generation != 0)

# Making Column Names Better
colnames(Pokemon) <- c("Number", "Name", "MainType", "SecondaryType", 
                       "Total", "HP", "Attack", "Defense", "SP.Attack", "SP.Defense", "Speed",
                       "Generation", "Legendary")

##### Server Construction
server <- function(input, output){
  
  # Plot 1 Settings
  selectedData <- reactive({
  Pokemon.div <- tapply(X = Pokemon[, input$col], INDEX = if(input$type){
    Pokemon$MainType} else{Pokemon$SecondaryType}, FUN = mean)
  Pokemon.div <- data.frame(types = names(Pokemon.div), mean.stat = Pokemon.div)
  Pokemon.div$stat.Z <- (Pokemon.div$mean.stat - mean(Pokemon.div$mean.stat)) / sd(Pokemon.div$mean.stat)
  Pokemon.div$avg <- ifelse(Pokemon.div$stat.Z > 0, "Above Average", "Below Average")
  Pokemon.div <- Pokemon.div[order(Pokemon.div$stat.Z),]
  Pokemon.div$types <- factor(Pokemon.div$types, levels = Pokemon.div$types)
  Pokemon.div
  })
  
  output$plot <- renderPlot({
  ggplot(selectedData(), aes(x = types, y = stat.Z, label = stat.Z)) + 
    geom_bar(aes(fill = avg), stat = 'identity', width = 0.5) +
      list(if(input$col.pal == "Default"){
    scale_fill_manual(name = "",
                      labels = c("Above Average", "Below Average"),
                      values = c("Above Average" = "#048c22", "Below Average" = "#f54040"))
        }else{
          scale_fill_discrete_diverging(palette = input$col.pal,
                                          name = "")
        }) +
    scale_x_discrete(name = ifelse(input$type, "Main Types", "Secondary Types")) +
    labs(y = "Standard Deviations",
         title = paste("Deviations from Pokemon Average", input$col,"Statistic for Each Main Type")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = c(0.85, 0.23),
          legend.key.size = unit(1.25, 'cm'),
          legend.background = element_blank(),
          text = element_text(size = 18))
  })
  
  # Plot 2 Settings
  selectedData2 <- reactive({
    if(input$xcol2 == "Stats"){
      Pokemon %>%
        select(c("HP", "Attack", "Defense", "SP.Attack", "SP.Defense", "Speed", input$ycol2)) %>%
        pivot_longer(1:6, names_to = "Stats", values_to = "Values")
    }else{
      if(input$ycol2 == "Stats"){
      Pokemon %>%
        select(c("HP", "Attack", "Defense", "SP.Attack", "SP.Defense", "Speed", input$xcol2)) %>%
        pivot_longer(1:6, names_to = "Stats", values_to = "Values")
    }else{
      if(input$xcol2 == "Generation" & input$ycol2 == "SecondaryType"){
      Pokemon %>%
        filter(SecondaryType != "") %>%
        count(SecondaryType, Generation)
    }else{
      if(input$ycol2 == "Generation" & input$xcol2 == "SecondaryType"){ 
      Pokemon %>%
        filter(SecondaryType != "") %>%
        count(SecondaryType, Generation)
    }else{
      if(input$xcol2 == "Generation" & input$ycol2 == "MainType"){ 
      Pokemon %>%
        count(MainType, Generation)
    }else{
      if(input$ycol2 == "Generation" & input$xcol2 == "MainType"){
      Pokemon %>%
        count(MainType, Generation)
    }else{
      if(input$ycol2 == "MainType" & input$xcol2 == "SecondaryType"){
      Pokemon %>%
        filter(SecondaryType != "") %>%
        count(MainType, SecondaryType)
    }else{
      if(input$ycol2 == "SecondaryType" & input$xcol2 == "MainType"){
      Pokemon %>%
        filter(SecondaryType != "") %>%
        count(MainType, SecondaryType)
    }else
      if(all(sapply(list(input$xcol2, input$ycol2), function(x) x == "Generation"))){
        Pokemon %>%
          count(Generation, Generation)
    }else
      if(all(sapply(list(input$xcol2, input$ycol2), function(x) x == "MainType"))){
        Pokemon %>%
          count(MainType, MainType)
    }else
      if(all(sapply(list(input$xcol2, input$ycol2), function(x) x == "SecondaryType"))){
        Pokemon %>%
          filter(SecondaryType != "") %>%
          count(SecondaryType, SecondaryType)
    }else{
      Pokemon %>%
        select(c("HP", "Attack", "Defense", "SP.Attack", "SP.Defense", "Speed")) %>%
        pivot_longer(1:6, names_to = "Stats", values_to = "Values")
      }
    }
  }
}
      }
    }
  }
}
})
  
  output$plot2 <- renderPlot({
    ggplot(data = selectedData2(), aes_string(x = input$xcol2, y = input$ycol2)) +
      geom_tile(aes(fill = 
      if(input$xcol2 == "Stats"){
        Values
      }else
        if(input$ycol2 == "Stats"){
          Values
        }else{
          n
        }), color = "white", lwd = 0.75, linetype = 1) +
      scale_fill_continuous_sequential(palette = input$col.pal2,
                                       name = if(input$xcol2 == "Stats"){"Values"}
                                       else if(input$ycol2 == "Stats"){"Values"}
                                       else("Count")) +
      labs(title = paste("Heatmap of", input$ycol2, "versus", input$xcol2)) +
      list(if(input$xcol2 == "Generation"){
        scale_x_continuous(name = "Generation",
                           breaks = seq(1, 8, by = 1))
      }else{
        scale_x_discrete(guide = guide_axis(n.dodge = 2))
      }) +
      coord_fixed() +
      theme_classic() +
      theme(text = element_text(size = 18))
  })
  
}

##### UI Construction

# Variable options for plot 1
vars <- colnames(Pokemon[, 5:11])
# Color Palette for plot 1
coolers <- c("Default","Blue-Red", "Blue-Red 2", "Blue-Red 3", "Red-Green", "Purple-Green", "Purple-Brown", 
             "Green-Brown", "Blue-Yellow 2", "Blue-Yellow 3", "Green-Orange", "Cyan-Magenta", 
             "Tropic", "Broc", "Cork", "Vik", "Berlin", "Lisbon", "Tofino")

# Variable options for plot 2
vars2 <- c("Stats", "MainType", "SecondaryType", "Generation")
# Color Palette for plot 2
coolers2 <- c("Grays", "Light Grays", "Blues 2", "Blues 3", "Purples 2", "Purples 3", "Reds 2",
             "Reds 3", "Greens 2", "Greens 3","Oslo", "Purple-Blue", "Red-Purple", "Red-Blue", 
             "Purple-Orange", "Purple-Yellow", "Blue-Yellow", "Green-Yellow","Red-Yellow", "Heat", 
             "Heat 2", "Terrain", "Terrain 2", "Viridis", "Plasma", "Inferno", "Rocket", "Mako", 
             "Dark Mint", "Mint", "BluGrn", "Teal", "TealGrn", "Emrld", "BluYl", "ag_GrnYl", "Peach", 
             "PinkYl", "Burg", "BurgYl", "RedOr", "OrYel", "Purp", "PurpOr", "Sunset", "Magenta", 
             "SunsetDark", "ag_Sunset", "BrwnYl", "YlOrRd", "YlOrBr", "OrRd", "Oranges", "YlGn", 
             "YlGnBu", "Reds", "RdPu", "PuRd", "Purples", "PuBuGn", "PuBu", "Greens", "BuGn", 
             "GnBu", "BuPu", "Blues", "Lajolla", "Turku", "Hawaii", "Batlow")

ui <- page_navbar(
  title = "Pokemon Statistic Comparisons",
  bg = "#2D89C8",
  inverse = TRUE,
  nav_panel(title = "Diverging Bar Plot",
  h5("Notice: Diverging bar plot and heatmap tabs above"),
            layout_sidebar(
              headerPanel(''),
              sidebar = sidebarPanel(
                width = 12,
                # Plot 1 Settings
                selectInput(inputId = 'col', label = 'Pokemon Statistic', choices = vars, selected = vars[[2]]),
                checkboxInput(inputId = "type", label = "Main Type / Secondary Type", value = TRUE),
                
                selectInput(inputId = 'col.pal', label = 'Color Palette', choices = coolers, selected = coolers[1])
              ),
              mainPanel(
                width = 12,
                plotOutput('plot', width = "100%", height = "600px")
              )
            ),
  p("Data related to the Pokemon video game focused on the Pokemon and their related statistics.")
  ),
  nav_panel(title = "Heatmap", 
  h5("Notice: Diverging bar plot and heatmap tabs above"),
            layout_sidebar(
              headerPanel(''),
              sidebar = sidebarPanel(
                width = 12,
                # Plot 2 Settings
                selectInput(inputId = 'xcol2', label = 'X Variable', choices = vars2, selected = vars2[[2]]),
                selectInput(inputId = 'ycol2', label = 'Y Variable', choices = vars2, selected = vars2[[1]]),
                
                selectInput(inputId = 'col.pal2', label = 'Color Palette', choices = coolers2, selected = coolers2[45])
              ),
              mainPanel(
                width = 12,
                plotOutput('plot2', width = "100%", height = "600px")
      )
    ),
  p("Data related to the Pokemon video game focused on the Pokemon and their related statistics.")
  )
)

# Complete App
shinyApp(ui, server)

