library(shiny)
library(igraph)
library(dplyr)
library(DT)
library(ggplot2)
library(shinythemes)

theme_set(theme_bw())

# Function to calculate uniform edge weights
## Every incoming edge of v with degree dv has weight 1/dv.
uniformWeights <- function(G) {
  # Initialize empty list to store edge weights
  Ew <- list()
  # Loop over edges in the graph
  for (e in E(G)) {
    # Get the target node of the edge
    v <- ends(G, e)[2]
    # Calculate the degree of the target node
    dv <- degree(G, v, mode = "in")
    # Assign weight to the edge
    Ew[[as.character(e)]] <- 1 / dv
  }
  return(Ew)
}

# Function to calculate random edge weights 
## Every edge has random weight. After weights assigned, we normalize weights of all incoming edges for each node so that they sum to 1.
randomWeights <- function(G) {
  Ew <- list()  # Initialize empty list to store edge weights
  # Assign random weights to edges
  for (v in V(G)) {
    in_edges <- incident(G, v, mode = "in")  # Get incoming edges for the current node
    ew <- runif(length(in_edges))  # Generate random weights for incoming edges
    total_weight <- sum(ew)  # Calculate the total weight of incoming edges
    # Normalize weights so that they sum to 1 for each node
    ew <- ew / total_weight
    # Store the weights for the incoming edges
    for (i in seq_along(in_edges)) {
      Ew[[as.character(in_edges[i])]] <- ew[i]
    }
  }
  return(Ew)
}


# Function to run linear threshold model
runLT <- function(G, S, Ew) {
  T <- unique(S)  # Targeted set with unique nodes
  lv <- sapply(V(G), function(u) runif(1))  # Threshold for nodes
  W <- rep(0, vcount(G))  # Weighted number of activated in-neighbors
  Sj <- unique(S)
  
  while (length(Sj) > 0) {
    if (length(T) >= vcount(G)) {
      break  # Break if the number of active nodes exceeds or equals the total number of nodes in G
    }
    Snew <- c()
    for (u in Sj) {
      neighbors <- neighbors(G, u, mode = "in")
      for (v in neighbors) {
        e <- as.character(get.edge.ids(G, c(v, u)))  # Define 'e' as the edge index
        if (!(v %in% T)) {
          # Calculate the total weight of the activated in-neighbors
          total_weight <- sum(Ew[[e]])
          
          # Update the weighted number of activated in-neighbors
          W[v] <- W[v] + total_weight
          
          # Check if the threshold is exceeded
          if (W[v] >= lv[v]) {
            Snew <- c(Snew, v)
            T <- c(T, v)
          }
        }
      }
    }
    Sj <- unique(Snew)  # Ensure unique nodes in the new set
  }
  return(T)  # Return all activated nodes
}


# Function to calculate the total number of active nodes at each iteration
activeNodes <- function(G, S, Ew, iterations) {
  active_df <- data.frame(iteration = integer(), 
                          total_active_nodes = integer())
  total_active_nodes <- rep(0, iterations)  # Initialize empty vector to store total active nodes
  
  for (i in 1:iterations) {
    T <- runLT(G, S, Ew)
    cat("--", i,"T:  ", T, "\n")
    total_active <- length(unique(T))  # Calculate the total active nodes in this iteration
    total_active_nodes[i] <- total_active  # Update total active nodes for current iteration
    
    # Limit total active nodes to the number of nodes in the graph
    if (total_active_nodes[i] >= vcount(G)) {
      total_active_nodes[i] <- vcount(G)  
    }
    
    # Update data frame with current iteration's total active nodes
    active_df <- rbind(active_df, data.frame(iteration = i, 
                                             total_active_nodes = total_active_nodes[i]))
    
    # Update seed set S for the next iteration
    S <- unique(c(S, T))
  }
  return(active_df)
}


#############################################################################################


# Shiny App
ui <- fluidPage(
  titlePanel("Linear Threshold Model Simulation :)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("graph_model", "Random Graph Model:",
                   choices = c("Erdős–Rényi" = "erdos_renyi", "Preferential Attachment" = "preferential"),
                   selected = "erdos_renyi"),
      conditionalPanel(
        condition = "input.graph_model == 'erdos_renyi'",
        sliderInput("n_er", "Number of Nodes (n):", min = 50, max = 500, value = 50),
        sliderInput("p_er", "Edge Probability (p):", min = 0.0001, max = 1, value = 0.05), 
        sliderInput("iterations", "Number of Iterations:", min = 1, max = 50, value = 10),
      ),
      conditionalPanel(
        condition = "input.graph_model == 'preferential'",
        sliderInput("n_pa", "Number of Nodes (n):", min = 50, max = 500, value = 50),
        sliderInput("power", "Power:", min = 1, max = 4, value = 2.3, step = 0.1),
        sliderInput("m", "Number of Edges (m):", min = 1, max = 50, value = 15), 
        sliderInput("iterations", "Number of Iterations:", min = 1, max = 50, value = 30),
      ),
      sliderInput("seed_nodes", "Number of Initial Seed Nodes:", min = 1, max = 10, value = 3),
      radioButtons("edge_weights", "Edge Weights:",
                   choices = c("Uniform" = "uniform", "Random" = "random"),
                   selected = "uniform"),
      actionButton("run_simulation", "Run Simulation")
    ),
    mainPanel(
      plotOutput("active_nodes_plot"),
      # textOutput(outputId = "active_df_table_title"),
      dataTableOutput("active_df_table")
    )
  ),
  theme = shinytheme("lumen")
)


server <- function(input, output) {
  observeEvent(input$run_simulation, {
    if (input$graph_model == "erdos_renyi") {
      random_graph <- erdos.renyi.game(input$n_er, p = input$p_er, directed = TRUE)
    } else {
      random_graph <- sample_pa(n = input$n_pa, power = input$power, m = input$m, directed = TRUE)
    }

    S <- sample(1:vcount(random_graph), input$seed_nodes)

    if (input$edge_weights == "uniform") {
      Ew <- uniformWeights(random_graph)
    } else {
      Ew <- randomWeights(random_graph)
    }

    active_df_list <- lapply(1:5, function(i) {
      activeNodes(random_graph, S, Ew, iterations = input$iterations)
    })

    # output$active_nodes_plot <- renderPlot({
    #   plot(NULL, xlim = c(1, input$iterations), ylim = c(0, max(sapply(active_df_list, max))),
    #        xlab = "Iteration", ylab = "Total Active Nodes", main = "Active Nodes over Iterations")
    #   for (i in 1:5) {
    #     lines(active_df_list[[i]]$iteration, active_df_list[[i]]$total_active_nodes, col = rainbow(5)[i])
    #   }
    #   legend("bottomright", legend = paste("df", 1:5), col = rainbow(5), lty = 1)
    # })

    active_df <- active_df_list[[1]] %>%
      left_join(active_df_list[[2]], by = "iteration") %>%
      left_join(active_df_list[[3]], by = "iteration") %>%
      left_join(active_df_list[[4]], by = "iteration") %>%
      left_join(active_df_list[[5]], by = "iteration") %>%
      rename(df1 = total_active_nodes.x,
             df2 = total_active_nodes.y,
             df3 = total_active_nodes.x.x,
             df4 = total_active_nodes.y.y,
             df5 = total_active_nodes)

    output$active_nodes_plot <- renderPlot({ 
      ggplot(active_df, aes(x = iteration)) +
        geom_line(aes(y = df1, color = "df1"), linetype = "solid") + 
        geom_line(aes(y = df2, color = "df2"), linetype = "solid") + 
        geom_line(aes(y = df3, color = "df3"), linetype = "solid") + 
        geom_line(aes(y = df4, color = "df4"), linetype = "solid") + 
        geom_line(aes(y = df5, color = "df5"), linetype = "solid") +
        scale_color_manual(values = c("#5E71C2", "#454655", "#A9AABC", "#C0535D", "#FC888F")) +
        ylab("Total Active Nodes") + 
        xlab("Iteration") + 
        scale_x_continuous(breaks = seq(min(active_df$iteration), max(active_df$iteration), by = 1)) +
        labs(color = "Data", 
             title = "Active Nodes over Iterations") +
        theme(legend.position = c(0.97, 0.02),
              legend.justification = c(1, 0),
              legend.box.background = element_rect(color = "black", linewidth = 0.5),
              legend.box.just = "top") 
    })
    
    
    # output$active_df_table_title <- renderText({
    #   "Summary of Active Nodes over Iterations"
    # })

    output$active_df_table <- renderDataTable({
      active_df
    }, rownames = FALSE, class = "hover")
  })
}


shinyApp(ui, server)