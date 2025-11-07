# app.R
library(shiny)
library(ape)
library(phytools)
library(data.table)
library(plotly)
library(foreach)
library(doParallel)

# -----------------------
# 1. Function definitions
# -----------------------
analyse_length_desc_subtrees <- function(tree) {
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  internal_nodes <- (Ntip + 1):(Ntip + Nnode)
  
  out <- lapply(internal_nodes, function(node) {
    desc_tips <- phytools::getDescendants(tree, node)
    desc_tips <- desc_tips[desc_tips <= Ntip]
    
    freq <- length(desc_tips)
    subtree <- keep.tip(tree, tree$tip.label[desc_tips])
    subtree_br_len <- sum(subtree$edge.length)
    
    data.frame(node = node, freq = freq, subtree_br_len = subtree_br_len)
  })
  
  return(rbindlist(out))
}

simulate_freq_vs_diversity <- function(n_genomes, n_trees, seed = 42, num_core = 2) {
  set.seed(seed)
  
  # Set up parallel cluster
  cl <- parallel::makeCluster(num_core)
  doParallel::registerDoParallel(cl)
  
  # Export custom functions to workers
  parallel::clusterExport(cl, varlist = c("analyse_length_desc_subtrees"), envir = environment())
  
  # Run in parallel
  results <- foreach(t = seq_len(n_trees), .combine = "rbind",
                     .packages = c("ape", "phytools", "data.table")) %dopar% {
                       tree <- rcoal(n_genomes)
                       tree_nodes_dt <- analyse_length_desc_subtrees(tree)
                       tree_nodes_dt$tree_id <- t
                       # wrap inside a list so .combine=rbindlist can concatenate properly
                       tree_nodes_dt[, c("tree_id", "node", "freq", "subtree_br_len")]
                     }
  
  # Clean up
  parallel::stopCluster(cl)
  # foreach::registerDoSEQ()
  
  return(results)
}

# -----------------------
# 2. Shiny UI
# -----------------------
ui <- fluidPage(
  titlePanel("Accessory Gene Frequency vs Subtree Length"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_genomes", "Number of genomes:", 93, min = 5, max = 500),
      numericInput("n_trees", "Number of simulated trees:", 50, min = 10, max = 2000),
      numericInput("seed", "Random seed:", 42),
      numericInput("num_core", "Number of threads:", 12),
      actionButton("run_sim", "Run Simulation"),
      br(), br(),
      helpText("This app simulates random coalescent trees (using ape::rcoal) 
               and visualizes the relationship between subtree branch length 
               and accessory gene frequency.")
    ),
    
    mainPanel(
      plotlyOutput("plot", height = "600px"),
      br(),
      textOutput("summary")
    )
  )
)

# -----------------------
# 3. Server logic
# -----------------------
server <- function(input, output, session) {
  sim_data <- eventReactive(input$run_sim, {
    withProgress(message = "Simulating trees...", value = 0, {
      # run the simulation
      res <- simulate_freq_vs_diversity(
        n_genomes = input$n_genomes,
        n_trees   = input$n_trees,
        seed      = input$seed,
        num_core  = input$num_core
      )
      incProgress(1) # mark as complete
      res
    })
  })
  
  output$plot <- renderPlotly({
    req(sim_data())
    res <- sim_data()
    plot_ly(
      data = res,
      x = ~subtree_br_len,
      y = ~freq,
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 5, color = 'rgba(55,128,191,0.6)'),
      hoverinfo = 'text',
      text = ~paste('Tree:', tree_id, '<br>Node:', node,
                    '<br>Freq:', freq, '<br>Branch length:', round(subtree_br_len, 3))
    ) %>%
      layout(
        title = "Accessory Gene Frequency vs Subtree Length",
        xaxis = list(title = "Subtree total branch length"),
        yaxis = list(title = "Accessory gene frequency")
      )
  })
  
  output$summary <- renderText({
    req(sim_data())
    paste("Simulated", input$n_trees, "trees with", input$n_genomes, "genomes each.")
  })
}


# -----------------------
# 4. Run app
# -----------------------
shinyApp(ui, server)

