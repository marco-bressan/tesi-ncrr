## code to prepare `meta_lu_ades` dataset goes here
xx <- read.csv("data-raw/csv/smoking_cessation.csv")
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

xx <- read.csv("data-raw/csv/thrombolytic_drugs.csv")
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

# smoke alarm data (achana et al. ex 1)

yy <- read.csv("data-raw/csv/smoke_alarm", header = FALSE,
               stringsAsFactors = FALSE)
study <- yy[, 1]
study.id <- seq_along(yy[, 1])
trts <- tolower(do.call(rbind, strsplit(yy[, 2], " versus ", fixed = TRUE)))
adj.risk <- stringr::str_extract_all(yy[, 3], r"{\(([0-9.]+)\)}") #???
raw.risk <- do.call(rbind, strsplit(gsub(r"{\( ?[0-9.]+\)}", "", yy[, 3]), "/| vs. "))

smoke.alarm <- rbind(
  data.frame(study.id = study.id, treatment = trts[, 1], is.baseline = 1,
             rik = raw.risk[, 1], nik = raw.risk[, 2]),
  data.frame(study.id = study.id, treatment = trts[, 2], is.baseline = 0,
             rik = raw.risk[, 3], nik = raw.risk[, 4])
)

smoke.alarm$treatment <- relevel(factor(trimws(smoke.alarm$treatment)), "usual care")
smoke.alarm$rik <- as.numeric(smoke.alarm$rik)
smoke.alarm$nik <- as.numeric(smoke.alarm$nik)

smoke.alarm <- smoke.alarm[order(smoke.alarm$study.id), ]


# morphine (achana)
yy <- read.csv("data-raw/csv/morphine.csv", header = TRUE,
               stringsAsFactors = FALSE)
study <- trimws(yy[, 1])
study.id <- seq_along(yy[, 1])

trts <- lapply(yy[, -1], \(x) {
  ss <- strsplit(x, "/", fixed = TRUE)
  l0 <- which(!lengths(ss))
  ss[l0] <- rep(list(c(NA, NA, NA)), length(l0))
  do.call(rbind, ss)
})

names(trts) <- NULL

morphine <- mapply(\(x, nm) {
  mode(x) <- "numeric"
  colnames(x) <- c("nik", "mik", "sik")
  data.frame(study.id = study.id, treatment = nm,
             is.baseline = as.numeric(nm == "placebo"), x)
}, trts, c("placebo", "paracetamol", "nsaid", "cox-2"), SIMPLIFY = FALSE) |>
  do.call(what = rbind)

naidx <- apply(morphine[, c("nik", "mik", "sik")], 1, \(x) all(is.na(x)))
morphine <- morphine[-which(naidx), ]
morphine <- morphine[order(morphine$study.id), ]
morphine$treatment <- relevel(factor(morphine$treatment), ref = "placebo")


usethis::use_data(smoking, overwrite = TRUE)
usethis::use_data(thromb, overwrite = TRUE)
usethis::use_data(smoke.alarm, overwrite = TRUE)
usethis::use_data(morphine, overwrite = TRUE)
