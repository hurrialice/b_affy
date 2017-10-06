library(plotly)
# plot locally
p <- plot_ly(midwest, x = ~percollege, color = ~state, type = "box")
p

# online plot
Sys.setenv("plotly_username"="hurrialice")
Sys.setenv("plotly_api_key"="CYOeIoOFiGxEFy8sgcwc")



p <- plot_ly(midwest, x = ~percollege, color = ~state, type = "box")
api_create(p, filename = "midwest-boxplots")



mtcars$am[which(mtcars$am == 0)] <- 'Automatic'
mtcars$am[which(mtcars$am == 1)] <- 'Manual'
mtcars$am <- as.factor(mtcars$am)

p <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E')) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'Weight'),
                      yaxis = list(title = 'Gross horsepower'),
                      zaxis = list(title = '1/4 mile time')))

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = api_create(p, filename="basic")
chart_link


library(d3heatmap)
url <- "http://datasets.flowingdata.com/ppg2008.csv"
nba_players <- read.csv(url, row.names = 1)
d3heatmap(nba_players, scale = "column")
