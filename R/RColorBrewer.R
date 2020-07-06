### QUICK CODE FOR RCOLORBREWER

require(RColorBrewer)

# show all palettes
display.brewer.all(colorblindFriendly = TRUE)

# pick one palette
display.brewer.pal(n = 4, name = "Dark2")

# return the hex code
hex <- brewer.pal(n = 8, name = "Dark2")

# save the desired colors to your variables
treatmentcols <- c("HDG" = hex[1], "LDG" =  hex[2], "NG" = hex[3])
timecols <- c("PRE" = hex[4], "24H" = hex[5], "1WK" = hex[6], "4WK" = hex[7])

# accent colors
accentcols <- c(hex[8:12])
