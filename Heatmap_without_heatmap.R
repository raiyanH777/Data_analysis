# Filter relevant columns
hdata <- rt_data[, c("Patient ID", "Sex", "As Risk level", "Hematuria/Urine infection", 
                     "Cigarette/Betel Leaf habit", "Recurrent /Primary", "Corrected_Grade", 
                     "Tumor DNA_ADGRG6 mutation Status", "Professional hazard if any", 
                     "Invasiveness", "Mutation","Survival_status","Pathologic_Stage")]


# Convert necessary variables to factors for consistent ordering
hdata$`Tumor DNA_ADGRG6 mutation Status` <- as.factor(hdata$`Tumor DNA_ADGRG6 mutation Status`)
hdata$`Recurrent /Primary` <- as.factor(hdata$`Recurrent /Primary`)
hdata$Mutation <- as.factor(hdata$Mutation)
hdata$`As Risk level` <- as.factor(hdata$`As Risk level`)

# Heatmap Clustering and Annotation
dt <- t(as.matrix(hdata[, -1])) # Exclude Patient ID from clustering
order_by_sex <- order(hdata$Sex)
dt_ordered <- dt[, order_by_sex]
hdata_ordered <- hdata[order_by_sex, ]

# Define Annotations
col_annotations <- list(
  HeatmapAnnotation(
    Sex = hdata_ordered$Sex,
    col = list(Sex = c("Male" = "#2CB5C0", "Female" = "#FFB0B0")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    Pathologic_Stage = hdata_ordered$Pathologic_Stage,
    col = list(Pathologic_Stage = c("pT1" = "#3ABEF9", "pT2" = "#C65BCF")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  #HeatmapAnnotation(
   # Mutation = hdata_ordered$Mutation,
    #col = list(Mutation = c("Wild" = "#3ABEF9", "Mutant" = "#C65BCF")),
    #border = TRUE,
    #gp = gpar(col = "black")
  #),
  HeatmapAnnotation(
    Grade = hdata_ordered$Corrected_Grade,
    Grade = c("High" = "#640D5F", "Low" = "#D2FF72"),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    `Recurrent /Primary` = hdata_ordered$`Recurrent /Primary`,
    col = list(`Recurrent /Primary` = c("Primary" = "darkslategray2", "Recurrent" = "#fd306e")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    `Cigarette/Betel Leaf` = hdata_ordered$`Cigarette/Betel Leaf habit`,
    col = list(`Cigarette/Betel Leaf` = c("yes" = "white", "No" = "#3c3c64")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    `Arsenic Risk Level` = hdata_ordered$`As Risk level`,
    col = list(`Arsenic Risk Level` = c("High" = "#3ABEF9", "Moderate to High" = "#C65BCF", "Moderate" = "#640D5F", "Low" = "#47663B")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
 # HeatmapAnnotation(
  #  `Tumor DNA_ADGRG6 mutation Status` = hdata_ordered$`Tumor DNA_ADGRG6 mutation Status`,
   # col = list(`Tumor DNA_ADGRG6 mutation Status` = c("WT" = "#F9D4B7", "ADGRG6-209T" = "#B983A7", "ADGRG6-206A/209T" = "#E44144", "ADGRG6-206A" = "#41B7C4")),
    #border = TRUE,
    #gp = gpar(col = "black")
#  ),
  HeatmapAnnotation(
    `Professional Hazard` = hdata_ordered$`Professional hazard if any`,
    col = list(`Professional Hazard` = c("Yes" = "#3ABEF9", "No" = "#C65BCF")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    Invasiveness = hdata_ordered$Invasiveness,
    col = list(Invasiveness = c("Invasive" = "#3ABEF9", "Noninvasive" = "#C65BCF")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    Hematuria = hdata_ordered$`Hematuria/Urine infection`,
    col = list(Hematuria = c("No" = "#3c3c64", "Yes" = "#c39fc9")),
    border = TRUE,
    gp = gpar(col = "black")
  ),
  HeatmapAnnotation(
    Survival_status = hdata_ordered$Survival_status,
    col = list(Survival_status = c("Died" = "#3c3c64", "Alive recurred" = "#c39fc9","Alive" = "darkslategray2")),
    border = TRUE,
    gp = gpar(col = "black")
  )
)

# Generate Heatmap
ht <- Heatmap(
  dt_ordered,
  name = "Expression",
  top_annotation = do.call(c, col_annotations),
  col = c("Male" = "#FFFFFF00", "Female" = "#FFFFFF00"),
  show_heatmap_legend = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  rect_gp = gpar(type = "none"),
  column_title = "Bangladesh Bladder Cancer Cohort",
  width = unit(20, "cm"),
  column_split = hdata_ordered$Sex
)
ht
