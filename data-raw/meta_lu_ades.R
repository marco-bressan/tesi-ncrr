## code to prepare `meta_lu_ades` dataset goes here
xx <- read.csv("data-raw/smoking_cessation.csv")
smoking <- apply(xx, 1, \(x) {
  trt <- as.character(x[-(1:2)])
  ctrl <- rep(0, length(trt))
  ctrl[as.integer(x[1])] <- 1
  dd <- data.frame(study.id = as.integer(x[2]), treatment = names(xx)[-(1:2)],
                   cases = gsub(",", "", x = trt),
                   is.baseline = ctrl)
  return(dd[dd$cases != "", ])
} , simplify = FALSE) |>
  do.call(rbind, args = _)

cases.mat <- strsplit(smoking[[3]], "/", fixed = TRUE) |>
  vapply(as.numeric, numeric(2)) |>
  t()
colnames(cases.mat) <- c("rik", "nik")
smoking <- cbind(smoking[-3], cases.mat)
lvlnm <- c("No.contact..A." = "A - No contact",
           "Self.help..B." = "B - Self help",
           "Individual.counseling..C." = "C - Individual counseling",
           "Group.counseling..D." = "D - Group counseling")
smoking$treatment <- stringr::str_replace_all(smoking$treatment, lvlnm) |>
  factor(levels = lvlnm)

xx <- read.csv("data-raw/thrombolytic_drugs.csv")
thromb <- apply(xx, 1, \(x) {
  trt <- as.character(x[-(1:2)])
  ctrl <- rep(0, length(trt))
  ctrl[as.integer(x[1])] <- 1
  dd <- data.frame(study.id = as.integer(x[2]), treatment = names(xx)[-(1:2)],
                   cases = gsub(",", "", x = trt),
                   is.baseline = ctrl)
  return(dd[dd$cases != "", ])
} , simplify = FALSE) |>
  do.call(rbind, args = _)

cases.mat <- strsplit(thromb[[3]], "/", fixed = TRUE) |>
  vapply(as.numeric, numeric(2)) |>
  t()
colnames(cases.mat) <- c("rik", "nik")
thromb <- cbind(thromb[-3], cases.mat)
lvlnm <- c("SK..1." = "1 - SK",
           "AtPA..2." = "2 - AtPA",
           "t...PA..3." = "3 - t-PA",
           "SK...tPA..4." = "4 - SK+tPA",
           "Ten..5." = "5 - Ten",
           "Ret..6." = "6 - Ret",
           "UK..7." = "7 - UK",
           "ASPAC..8." = "8 - ASPAC")
thromb$treatment <- stringr::str_replace_all(thromb$treatment, lvlnm) |>
  factor(levels = lvlnm)

usethis::use_data(smoking, overwrite = TRUE)
usethis::use_data(thromb, overwrite = TRUE)
