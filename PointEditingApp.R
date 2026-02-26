## App using shiny to allow users to edit road points interactively ## 

library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)

# make some simple example points
pts <- st_as_sf(data.frame(
  id = 1:5,
  lon = c(-123.1, -123.2, -123.3, -123.4, -123.5),
  lat = c(49.3, 49.31, 49.32, 49.33, 49.34)
), coords = c("lon", "lat"), crs = 4326)

# create ui layout
ui <- fluidPage(
  titlePanel("Editable Points App — Smooth Add/Delete with AOI"),
  fluidRow(
    column(8, leafletOutput("map", height = "650px")),   # Map panel
    column(4,                                            # Control panel
           # Mode selector: choose between adding or deleting individual points
           radioButtons("mode", "Click mode:",
                        choices = c("Add points" = "add",
                                    "Delete points" = "delete"),
                        selected = "add"),
           # Button to delete all points inside a drawn Area of Interest (AOI)
           actionButton("delete_aoi", "Delete points in AOI"),
           br(), br(),
           # Reset button to restore the original set of points
           actionButton("reset", "Reset to original"),
           # Download edited points as GeoJSON
           downloadButton("download", "Download GeoJSON"),
           br(), br(),
           strong("AOI selection:"),
           textOutput("selected_count"),  # Shows how many points fall inside AOI
           br(),
           strong("Saved / current points:"),
           verbatimTextOutput("result")   # Displays current points in console format
    )
  )
)

# create server logic
server <- function(input, output, session) {
  # Reactive values store dynamic data that changes during interaction
  rv <- reactiveVal(pts)          # Stores current points
  aoi <- reactiveVal(NULL)        # Stores AOI polygon (if drawn)
  drawing_in_progress <- reactiveVal(FALSE) # Track if user is currently drawing
  map_view <- reactiveVal(list()) # Store current map view (center & zoom)
  
  # --- Initialize Leaflet map ---
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("Esri.WorldImagery") %>% # Satellite basemap
      addCircleMarkers(
        data = rv(),       # initial points
        radius = 6,
        color = "red",
        stroke = TRUE,
        fillOpacity = 0.9,
        layerId = ~as.character(id),  # marker ID = point ID
        group = "points",
        label = ~paste0("ID: ", id)
      ) %>%
      # Drawing tools for user to create AOI (polygon or rectangle)
      addDrawToolbar(
        targetGroup = "aoi_draw",
        polylineOptions = FALSE,
        circleOptions = FALSE,
        markerOptions = FALSE,
        circleMarkerOptions = FALSE,
        polygonOptions = drawPolygonOptions(showArea = TRUE, repeatMode = FALSE),
        rectangleOptions = drawRectangleOptions(repeatMode = FALSE),
        editOptions = editToolbarOptions()
      ) %>%
      # Add layer control to toggle visibility of point and AOI groups
      addLayersControl(overlayGroups = c("points", "aoi_draw"),
                       options = layersControlOptions(collapsed = TRUE))
  })
  
  # AOI drawing state tracking
  observeEvent(input$map_draw_start, { drawing_in_progress(TRUE) })  # User started drawing
  observeEvent(input$map_draw_stop, { drawing_in_progress(FALSE) })  # Drawing finished
  
  # Track map view (zoom & center) so we can restore it after updates
  observe({
    req(input$map_zoom, input$map_center)
    map_view(list(zoom = input$map_zoom, center = input$map_center))
  })
  
  # Add a new point by clicking on the map
  observeEvent(input$map_click, {
    req(input$mode == "add")          # Only when in'Add' mode
    req(!drawing_in_progress())       # Ignore clicks during AOI drawing
    
    click <- input$map_click          # Get click coordinates
    new_id <- ifelse(nrow(rv()) == 0, 1, max(rv()$id, na.rm = TRUE) + 1)
    
    # Create new point as an sf object
    new_point <- st_sf(
      id = new_id,
      geometry = st_sfc(st_point(c(click$lng, click$lat)), crs = 4326)
    )
    rv(rbind(rv(), new_point))        # Append to reactive dataset
    
    # Update map dynamically without full re-render
    view <- map_view()
    leafletProxy("map") %>%
      addCircleMarkers(
        data = new_point,
        radius = 6,
        color = "red",
        layerId = as.character(new_id),
        group = "points",
        label = paste0("ID: ", new_id)
      ) %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Delete a single point by clicking its marker
  observeEvent(input$map_marker_click, {
    req(input$mode == "delete")       # Only when in 'Delete' mode
    click <- input$map_marker_click
    rv(rv()[rv()$id != click$id, ])   # Remove clicked point from dataset
    
    # Update map view after deletion
    view <- map_view()
    leafletProxy("map") %>%
      removeMarker(layerId = click$id) %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Handle AOI polygon or rectangle creation
  observeEvent(input$map_draw_new_feature, {
    feat <- input$map_draw_new_feature
    geo <- feat$geometry
    if (!is.null(geo) && tolower(geo$type) %in% c("polygon", "multipolygon")) {
      coords <- geo$coordinates[[1]]
      mat <- do.call(rbind, lapply(coords, function(x) c(as.numeric(x[1]), as.numeric(x[2]))))
      poly <- st_polygon(list(mat))
      aoi_poly <- st_sfc(poly, crs = 4326)
      aoi(aoi_poly)  # Store AOI geometry for later use
      
      # Draw AOI outline on the map
      view <- map_view()
      leafletProxy("map") %>%
        clearGroup("aoi") %>%  # Clear any previous AOI
        addPolygons(data = aoi_poly, group = "aoi", color = "blue", weight = 2, fill = FALSE) %>%
        setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
    }
  })
  
  # Display number of points inside AOI
  output$selected_count <- renderText({
    if (is.null(aoi())) return("No AOI drawn")    # No AOI yet
    current <- rv()
    if (nrow(current) == 0) return("0 points in AOI")
    mat <- st_intersects(current, aoi(), sparse = FALSE)
    n_sel <- if (is.matrix(mat)) sum(mat[,1]) else 0
    paste0(n_sel, " point(s) in AOI")
  })
  
  # Delete all points inside AOI when button pressed
  observeEvent(input$delete_aoi, {
    req(!is.null(aoi()))   # Require AOI exists
    current <- rv()
    mat <- st_intersects(current, aoi(), sparse = FALSE)
    sel_idx <- which(mat[,1])         # Indices of points inside AOI
    if(length(sel_idx) == 0) return() # Nothing to delete
    rv(current[-sel_idx, , drop = FALSE])  # Keep only points outside AOI
    
    # Update map by removing deleted markers
    view <- map_view()
    leafletProxy("map") %>%
      removeMarker(layerId = as.character(current$id[sel_idx])) %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Reset points to original dataset 
  observeEvent(input$reset, {
    rv(pts)      # Restore original points
    aoi(NULL)    # Clear AOI
    view <- map_view()
    leafletProxy("map") %>%
      clearGroup("aoi") %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Download current points as GeoJSON file 
  output$download <- downloadHandler(
    filename = function() paste0("edited_points_", Sys.Date(), ".geojson"),
    content = function(file) st_write(rv(), file, driver = "GeoJSON", delete_dsn = TRUE)
  )
  
  # Show current sf object (debug/info output)
  output$result <- renderPrint({ rv() })
}

# Run the app
shinyApp(ui, server)