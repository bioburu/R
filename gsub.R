#--sub out strings inbetween 
gsub(".*[beginning]([^.]+)[end].*", "\\1",genes$ID)
