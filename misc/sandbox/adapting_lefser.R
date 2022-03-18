library(lefser)
data("zeller14")
se <- zeller14[, zeller14$study_condition != "adenoma"]
res <- lefser2(se, groupCol = "study_condition")

