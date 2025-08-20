devtools::load_all()
data("smoking")

object <- ncrr.design(smoking)


# NMA ==========

nj <- sqrt(object$nj)
M <- do.call(cbind, lapply(object$design, \(x) rbind(x[-1], x[1]))) + 1


edges <- as.data.frame(table(M[1,], M[2,]))
edges <- subset(edges, Freq > 0)
G <- igraph::graph_from_data_frame(edges, directed = FALSE)

igraph::V(G)$size <- 50 * sqrt((as.vector(nj) / max(nj)))
igraph::E(G)$width <- edges$Freq/2
igraph::V(G)$name <- object$treatments[as.numeric(names(igraph::V(G)))]

glay <- t(sapply(pi/2 + 1:length(nj) * 2*pi/length(nj), \(x) 2 * c(cos(x), sin(x))))
plot(G, vertex.label = igraph::V(G)$name, vertex.size = igraph::V(G)$size,
     layout = glay)
plot(ncrr.design(smoke.alarm), vertex.label.cex = 0.7, vertex.label.family = "Arial",
     pos.offset = 0.6, smart.layout.seed = NA)
labbe.plot(object, xylims = c(0, 0.7))

# L'abbè plot ==========

pct <- with(object$x, rik / nik)
lc <- cumsum(ll <- lengths(object$design))
lc <- c(0, lc[-length(lc)])
pct <- mapply(\(pos, len) {
  pp <- rbind(n = object$x$nik[(pos + 1):(pos + len)],
              pct = pct[(pos + 1):(pos + len)])
  if (ncol(pp) > 2) {
    pp <- cbind(n = unname(expand.grid(pp[1, 1], pp[1, -1], KEEP.OUT.ATTRS = FALSE)),
                pct = unname(expand.grid(pp[2, 1], pp[2, -1], KEEP.OUT.ATTRS = FALSE)))
  } else {
    pp <- cbind(n = t(pp[1, ]), pct = t(pp[2, ]))
    colnames(pp) <- paste(rep(c("n", "pct"), c(2, 2)), rep(1:2, 2), sep = ".")
  }
  return(pp)
} , lc, ll, SIMPLIFY = FALSE)
pct.data <- do.call(rbind, mapply(\(x, ...) cbind(as.data.frame(x), ...),
                                  x = pct, study.id = seq_along(pct),
                                  trt = lapply(object$design, \(x) paste(x[1], x[-1], sep = "-")),
                                  SIMPLIFY = FALSE))
pct.data$trt <- factor(pct.data$trt)
pct.data <- pct.data[order(pct.data$n.2, decreasing = TRUE), ]
plot(pct.2 ~ pct.1, data = pct.data, pch = 16, cex = 1 + 3 * (n.2 / max(n.2)),
     col = 1 + as.integer(trt), xlim = range(c(pct.1, pct.2)), ylim = range(c(pct.1, pct.2)))
points(pct.2 ~ pct.1, data = pct.data,
       pch = ifelse(substr(as.character(trt), 1, 1) == "0", 1, 13),
       cex = 1 + 3 * (n.2 / max(n.2)))
abline(a = 0, b = 1, lty = 2, col = "gray30")
legend("bottomright", col = c(1 + seq_along(levels(pct.data$trt)), 1, 1),
       pch = c(rep(16, length(levels(pct.data$trt))), 1, 13),
       legend = c(gsub("-", " vs. ", levels(pct.data$trt)), "baseline", "no baseline"))
