# Calculate Freeman-Tukey statistic
FT_stat <- function(realised, expected){
  sum((sqrt(realised) - sqrt(expected))^2)
}