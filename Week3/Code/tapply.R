x <- 1:20 # A vector

# A factor (of the same length) defining groups
y <- factor(rep(letters[1:5], each = 4))

tapply(x, y, sum)