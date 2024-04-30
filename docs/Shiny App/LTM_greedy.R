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


# Function to calculate average size of activated nodes
avgLT <- function(G, S, Ew, iterations = 1) {
  avgSize <- 0
  for (i in 1:iterations) {
    T <- runLT(G, S, Ew)
    avgSize <- avgSize + length(T) / iterations
  }
  return(avgSize)
}


# Define the Greedy_LTM function
Greedy_LTM <- function(G, Ew, k, iterations) {
  start <- Sys.time()  # Record the start time
  S <- c()  # Initialize the seed set
  
  for (i in 1:k) {
    inf <- data.frame(nodes = V(G), influence = NA)  # Initialize the influence table
    
    # Calculate the influence for nodes not in S
    for (v in V(G)) {
      if (!(v %in% S)) {
        inf$influence[v] <- avgLT(G, c(S, v), Ew, iterations = 1)
      }
    }
    
    # Exclude nodes already in S
    inf_excluded <- inf[!inf$nodes %in% S, ]
    
    # Select the node with maximum influence and add it to the seed set
    u <- inf_excluded[which.max(inf_excluded$influence), ]$nodes
    cat("Selected node:", u, "with influence:", max(inf_excluded$influence), "\n")
    
    # Convert node name to numeric
    u <- as.numeric(u)
    
    # Add selected node to the seed set
    S <- c(S, u)
  }
  
  end <- Sys.time()  # Record the end time
  # Print the total time taken
  print(paste("Total time:", end - start))
  
  return(S)  # Return the seed set
}



################################################################################################

# UI
ui <- fluidPage(
  titlePanel("Greedy Algorithm of Influence Max Problem on LTM :)"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("graph_model", "Random Graph Model:",
                   choices = c("Erdős–Rényi" = "erdos_renyi", "Preferential Attachment" = "preferential"),
                   selected = "erdos_renyi"),
      conditionalPanel(
        condition = "input.graph_model == 'erdos_renyi'",
        sliderInput("n_nodes", "Number of Nodes (n):", min = 50, max = 300, value = 50),
        sliderInput("p", "Edge Probability (p):", min = 0.0001, max = 1, value = 0.05), 
        sliderInput("iterations", "Number of Iterations:", min = 1, max = 50, value = 10)
      ),
      conditionalPanel(
        condition = "input.graph_model == 'preferential'",
        sliderInput("n_nodes", "Number of Nodes (n):", min = 50, max = 300, value = 50),
        sliderInput("power", "Power:", min = 1, max = 4, value = 2.3, step = 0.1),
        sliderInput("m", "Number of Edges (m):", min = 1, max = 50, value = 15), 
        sliderInput("iterations", "Number of Iterations:", min = 1, max = 50, value = 30)
      ),
      sliderInput("seed_nodes", "Number of Initial Seed Nodes:", min = 1, max = 10, value = 3),
      radioButtons("edge_weights", "Edge Weights:",
                   choices = c("Uniform" = "uniform", "Random" = "random"),
                   selected = "uniform"),
      actionButton("run_simulation", "Run Simulation")
    ),
    mainPanel(
      plotOutput("active_nodes_plot"),
      dataTableOutput("active_df_table")
    )
  ),
  theme = shinytheme("lumen")
)

# Server
server <- function(input, output) {
  observeEvent(input$run_simulation, {
    if (input$graph_model == "erdos_renyi") {
      random_graph <- erdos.renyi.game(input$n_nodes, p = input$p, directed = TRUE)
    } else {
      random_graph <- sample_pa(n = input$n_nodes, power = input$power, m = input$m, directed = TRUE)
    }
    
    if (input$edge_weights == "uniform") {
      Ew <- uniformWeights(random_graph)
    } else {
      Ew <- randomWeights(random_graph)
    }
    
    active_df_list <- lapply(1:3, function(i) {
      S <- sample(1:vcount(random_graph), input$seed_nodes)
      activeNodes(random_graph, S, Ew, iterations = input$iterations)
    })
    
    seed_set <- Greedy_LTM(random_graph, Ew, k = input$seed_nodes, iterations = input$iterations)
    
    active_df_selectedSeed <- activeNodes(random_graph, seed_set, Ew, iterations = input$iterations)
    
    active_df <- active_df_list[[1]] %>% 
      left_join(active_df_list[[2]], by = "iteration") %>% 
      left_join(active_df_list[[3]], by = "iteration") %>% 
      left_join(active_df_selectedSeed, by = "iteration") %>% 
      rename(df1 = total_active_nodes.x,
             df2 = total_active_nodes.y,
             df3 = total_active_nodes.x.x,
             greedy = total_active_nodes.y.y)
    
    output$active_nodes_plot <- renderPlot({
      ggplot(active_df, aes(x = iteration)) +
        geom_line(aes(y = df1, color = "df1"), linetype = "solid") +
        geom_line(aes(y = df2, color = "df2"), linetype = "solid") +
        geom_line(aes(y = df3, color = "df3"), linetype = "solid") +
        geom_line(aes(y = greedy, color = "greedy"), linetype = "solid") +
        scale_color_manual(values = c("black", "black", "black", "#FC888F")) +
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
    
    output$active_df_table <- renderDataTable({
      active_df
    }, rownames = FALSE, class = "hover")
  })
}

# Run the Shiny app
shinyApp(ui, server)