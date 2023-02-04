#https://tomazweiss.github.io/blog/object_detection/
library(image.darknet)
library(Rcpp)
library(dplyr)
library(tidyr)
library(here)
#---------------------User inputs-------------------------------------------
#Assign working directory
setwd("/home/amplified_prog/RstudioProjects/segment_image_exp")
#Assign image path
image_path <- here("/home/amplified_prog/RstudioProjects/segment_image_exp/jpegs")
#Assign threshold for predictions
Threshold <- 0.19
#Assign model parameters
model <- image_darknet_model(
      type = "detect",
      model = "tiny-yolo-voc.cfg",
      weights = system.file(package = "image.darknet", "models", "tiny-yolo-voc.weights"),
      labels = system.file(package = "image.darknet", "include", "darknet", "data", "voc.names")
)
#----------------------End user inputs--------------------------------------------------
print(model)
all_images <- dir(path = image_path, pattern = "\\.png|\\.jpg|\\.jpeg")
all_images
dir.create('predicted_images')
#apply a function to predict on all images in path
detect_function <- function(x) {
  filename <- paste(image_path, x, sep = "/") #assign path to all files/
  prediction <- image_darknet_detect( #predict on each file
    file = filename, #feed the files
    object = model, #upload model
    threshold = Threshold #set threshold for predictions
  )
  file.rename("predictions.png", paste0("predicted_images/", x))
  return(prediction)
}
#------------------------------------------------------------------------------------
#write the cpp stdout to a .txt file
cppFunction('void redir(){FILE* F=freopen("capture.txt","w+",stdout);}')
#redirect cpp stdout back to console
cppFunction('void resetredir(){FILE* F=freopen("CON","w+",stdout);}')
redir();
#--------------------redirect cpp stdout to txt file-------------------------------------------------
#enter all images into model for predictions
std_out <- lapply(all_images, detect_function)
std_out
#---------------------cpp stdout back to console----------------------
resetredir();
#read txt file into data frame
std_out <- data.frame(txt = unlist(readLines("capture.txt")))
head(std_out)
#removes Boxes and any unnamed text
std_out_filtered <- std_out %>%
      filter(!grepl("Boxes", txt)) %>%
      filter(!grepl("pandoc", txt)) %>%
      filter(!grepl("unnamed", txt))
head(std_out_filtered)
std_out_filtered$txt
#make a logical vector column of rows with filenames=TRUE 
grepl(image_path, std_out_filtered$txt)
std_out_filtered$isfile <- grepl(image_path, std_out_filtered$txt)
head(std_out_filtered)
#remove path names and keep only file names
image_path
gsub(paste0(image_path, '/'), "", std_out_filtered$txt)
std_out_filtered$txt <- gsub(paste0(image_path, '/'), "", std_out_filtered$txt)
#replace isfile column with file names
ifelse(std_out_filtered$isfile, std_out_filtered$txt, NA)
#make file column
std_out_filtered$file <- ifelse(std_out_filtered$isfile, std_out_filtered$txt, NA)
#add a column of predicted objects
ifelse(!std_out_filtered$isfile, std_out_filtered$txt, NA)
std_out_filtered$predicted <- ifelse(!std_out_filtered$isfile, std_out_filtered$txt, NA)
data<- tidyr::fill(std_out_filtered, "file")
data
## Take out NAs and select the last two columns
data <- na.omit(data)[, 3:4]
data
# Separate the text that is held in two parts
data <- data %>% separate(file, into = c("file", "time"), sep = ":")
data
data <- data %>% separate(predicted, into = c("predicted", "prob"), sep = ":")

data <- data %>% filter(!is.na(prob))

# Keep only the prediction time
data$time <- gsub("Predicted in (.+).$", "\\1", data$time)

# Convert probabilities to numbers
data$prob <- as.numeric(sub("%", "", data$prob)) / 100

data %>% knitr::kable()

View(data)
# Optionally remove the file
file.remove("capture.txt")

break
View(std_out_filtered)
