library(shiny)
library(leaflet)
library(leaflet.extras)
library(sf)

# Initialize with empty sf object
pts <- st_sf(
  id = integer(0),
  geometry = st_sfc(crs = 4326)
)

# create ui layout
ui <- fluidPage(
  # Add custom JavaScript to clear drawn shapes
  tags$head(tags$script(HTML("
    Shiny.addCustomMessageHandler('clearDrawnShapes', function(message) {
      var map = HTMLWidgets.find('#map').getMap();
      if (map && map.eachLayer) {
        map.eachLayer(function(layer) {
          if (layer instanceof L.Path && layer.options.edit) {
            map.removeLayer(layer);
          }
        });
      }
    });
  "))),
  
  titlePanel("Editable Points App — Smooth Add/Delete with AOI"),
  fluidRow(
    column(8, leafletOutput("map", height = "650px")),
    column(4,
           strong("1. Upload GeoJSON file:"),
           fileInput("upload_geojson", NULL, accept = c(".geojson", ".json")),
           actionButton("load_upload", "Load uploaded file"),
           br(), br(),
           
           strong("2. Edit points:"),
           radioButtons("mode", "Click mode:",
                        choices = c("Add points" = "add",
                                    "Delete points" = "delete"),
                        selected = "add"),
           actionButton("delete_aoi", "Delete points in AOI"),
           br(), br(),
           actionButton("reset", "Reset to uploaded"),
           downloadButton("download", "Download GeoJSON"),
           br(), br(),
           strong("AOI selection:"),
           textOutput("selected_count"),
           br(),
           strong("Saved / current points:"),
           verbatimTextOutput("result")
    )
  )
)

# create server logic
server <- function(input, output, session) {
  rv <- reactiveVal(pts)
  seed_pts <- reactiveVal(pts)
  aoi <- reactiveVal(NULL)
  drawing_in_progress <- reactiveVal(FALSE)
  map_view <- reactiveVal(list())
  
  # --- Initialize Leaflet map (static - only renders once) ---
  output$map <- renderLeaflet({
    leaflet() %>%
      addProviderTiles("Esri.WorldImagery") %>%
      setView(lng = -123.3, lat = 49.3, zoom = 10) %>%
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
      addLayersControl(overlayGroups = c("points", "aoi_draw"),
                       options = layersControlOptions(collapsed = TRUE))
  })
  
  # Load uploaded GeoJSON file
  observeEvent(input$load_upload, {
    req(input$upload_geojson)
    
    tryCatch({
      uploaded <- st_read(input$upload_geojson$datapath, quiet = TRUE)
      
      if (!"id" %in% names(uploaded)) {
        uploaded$id <- seq_len(nrow(uploaded))
      }
      
      if (is.na(st_crs(uploaded))) {
        st_crs(uploaded) <- 4326
      } else if (st_crs(uploaded)$epsg != 4326) {
        uploaded <- st_transform(uploaded, 4326)
      }
      
      # Update reactive values
      rv(uploaded)
      seed_pts(uploaded)
      aoi(NULL)
      
      # Get bounding box
      bbox <- st_bbox(uploaded)
      
      # Use leafletProxy with clearMarkers and fitBounds in sequence
      leafletProxy("map", session) %>%
        clearMarkers() %>%
        clearGroup("points") %>%
        clearGroup("aoi") %>%
        addCircleMarkers(
          data = uploaded,
          lng = ~st_coordinates(geometry)[,1],
          lat = ~st_coordinates(geometry)[,2],
          radius = 6,
          color = "red",
          stroke = TRUE,
          fillOpacity = 0.9,
          layerId = ~as.character(id),
          group = "points",
          label = ~paste0("ID: ", id)
        ) %>%
        fitBounds(
          lng1 = bbox["xmin"],
          lat1 = bbox["ymin"],
          lng2 = bbox["xmax"],
          lat2 = bbox["ymax"],
          options = list(padding = c(50, 50), maxZoom = 18)
        )
      
      showNotification(
        paste0("Loaded ", nrow(uploaded), " points"), 
        type = "message"
      )
      
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
      print(e)
    })
  })
  
  # AOI drawing state tracking
  observeEvent(input$map_draw_start, { drawing_in_progress(TRUE) })
  observeEvent(input$map_draw_stop, { drawing_in_progress(FALSE) })
  
  # Track map view
  observe({
    req(input$map_zoom, input$map_center)
    map_view(list(zoom = input$map_zoom, center = input$map_center))
  })
  
  # Add a new point by clicking on the map
  observeEvent(input$map_click, {
    req(input$mode == "add")
    req(!drawing_in_progress())
    
    click <- input$map_click
    current <- rv()
    crs_cur <- st_crs(current)
    
    new_id <- if (nrow(current) == 0) 1 else max(current$id, na.rm = TRUE) + 1
    
    new_point <- st_sf(
      id = new_id,
      geometry = st_sfc(st_point(c(click$lng, click$lat)), crs = crs_cur)
    )
    
    # Harmonize columns
    missing_cols <- setdiff(names(current), names(new_point))
    for (nm in missing_cols) {
      new_point[[nm]] <- NA
    }
    
    extra_cols <- setdiff(names(new_point), names(current))
    if (length(extra_cols)) {
      new_point <- new_point[, setdiff(names(new_point), extra_cols), drop = FALSE]
    }
    
    new_point <- new_point[, names(current), drop = FALSE]
    rv(rbind(current, new_point))
    
    coords <- st_coordinates(new_point)
    
    view <- map_view()
    leafletProxy("map") %>%
      addCircleMarkers(
        lng = coords[1],
        lat = coords[2],
        radius = 6,
        color = "red",
        stroke = TRUE,
        fillOpacity = 0.9,
        layerId = as.character(new_id),
        group = "points",
        label = paste0("ID: ", new_id)
      ) %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Delete a single point by clicking
  observeEvent(input$map_marker_click, {
    req(input$mode == "delete")
    click <- input$map_marker_click
    rv(rv()[rv()$id != click$id, ])
    
    view <- map_view()
    leafletProxy("map") %>%
      removeMarker(layerId = click$id) %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
  })
  
  # Handle AOI polygon creation
  observeEvent(input$map_draw_new_feature, {
    feat <- input$map_draw_new_feature
    geo <- feat$geometry
    if (!is.null(geo) && tolower(geo$type) %in% c("polygon", "multipolygon")) {
      coords <- geo$coordinates[[1]]
      mat <- do.call(rbind, lapply(coords, function(x) c(as.numeric(x[1]), as.numeric(x[2]))))
      poly <- st_polygon(list(mat))
      aoi_poly <- st_sfc(poly, crs = 4326)
      aoi(aoi_poly)
      
      view <- map_view()
      leafletProxy("map") %>%
        clearGroup("aoi") %>%
        addPolygons(data = aoi_poly, group = "aoi", color = "blue", weight = 2, fill = FALSE) %>%
        setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
    }
  })
  
  # Display number of points in AOI
  output$selected_count <- renderText({
    if (is.null(aoi())) return("No AOI drawn")
    current <- rv()
    if (nrow(current) == 0) return("0 points in AOI")
    mat <- st_intersects(current, aoi(), sparse = FALSE)
    n_sel <- if (is.matrix(mat)) sum(mat[,1]) else 0
    paste0(n_sel, " point(s) in AOI")
  })
  
  # Delete points in AOI
  observeEvent(input$delete_aoi, {
    req(!is.null(aoi()))
    current <- rv()
    if (nrow(current) == 0) return()
    mat <- st_intersects(current, aoi(), sparse = FALSE)
    sel_idx <- which(mat[,1])
    if(length(sel_idx) == 0) return()
    rv(current[-sel_idx, , drop = FALSE])
    
    # Clear AOI reactive value
    aoi(NULL)
    
    view <- map_view()
    leafletProxy("map") %>%
      removeMarker(layerId = as.character(current$id[sel_idx])) %>%
      clearGroup("aoi") %>%
      setView(lng = view$center$lng, lat = view$center$lat, zoom = view$zoom)
    
    # Send custom message to clear drawn shapes (the blue box)
    session$sendCustomMessage(type = "clearDrawnShapes", message = list())
    
    showNotification("Deleted points in AOI and cleared polygon", type = "message")
  })
  
  # Reset to uploaded file
  observeEvent(input$reset, {
    uploaded <- seed_pts()
    rv(uploaded)
    aoi(NULL)
    
    if (nrow(uploaded) > 0) {
      bbox <- st_bbox(uploaded)
      
      leafletProxy("map") %>%
        clearMarkers() %>%
        clearGroup("points") %>%
        clearGroup("aoi") %>%
        addCircleMarkers(
          data = uploaded,
          lng = ~st_coordinates(geometry)[,1],
          lat = ~st_coordinates(geometry)[,2],
          radius = 6,
          color = "red",
          stroke = TRUE,
          fillOpacity = 0.9,
          layerId = ~as.character(id),
          group = "points",
          label = ~paste0("ID: ", id)
        ) %>%
        fitBounds(
          lng1 = bbox["xmin"],
          lat1 = bbox["ymin"],
          lng2 = bbox["xmax"],
          lat2 = bbox["ymax"],
          options = list(padding = c(50, 50), maxZoom = 18)
        )
    }
  })
  
  # Download current points
  output$download <- downloadHandler(
    filename = function() paste0("edited_points_", Sys.Date(), ".geojson"),
    content = function(file) {
      req(nrow(rv()) > 0)
      st_write(rv(), file, driver = "GeoJSON", delete_dsn = TRUE)
    }
  )
  
  # Show current points
  output$result <- renderPrint({ rv() })
}

# Run the app
shinyApp(ui, server)