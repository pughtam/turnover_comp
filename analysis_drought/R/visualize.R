library(xlsx)
library(ggplot2)
library(data.table)

load("data/dm_events.RData")
load("data/mortEvents.RData")

## for (p in names(mortEvent[[1]]@ts)) {
##   match = NULL
##   random.match = NULL
##   for (dme in names(dm.events)) {
##     if (length(unique(mortEvent[[dme]]@ts[[p]][,'mort.frac'])) < 2) {
##       next
##     }
##     match = append(match, mortEvent[[dme]]@wilcox[p])
##     random.match = append(random.match, mortEvent[[dme]]@random.wilcox[[p]])
##   }
  
##   model.summary[[m]][[d]] <- rbind(model.summary[[m]][[d]], c(event=sum(match < 0.05) / length(match),
##                                                               random=sum(random.match < 0.05) / length(random.match)))
## }
## row.names(model.summary[[m]][[d]]) <- names(mortEvent[[1]]@ts)

fout <- "mortality_summary.xlsx"
for (mn in names(model.summary)) {
  message(mn)
  yy = NULL
  for (dn in names(model.summary[[mn]])) {
    message(dn)
    xx=model.summary[[mn]][[dn]]
    rownames(xx) = sub("cmort", "cmort_", rownames(xx))
    rownames(xx) = sub("_$", "", rownames(xx))
    rownames(xx) = sub("cmort_", "", rownames(xx))
    rownames(xx) = sub("^_", "", rownames(xx))
    xx = as.data.frame(xx[, c("n.total", "event", "random")])
    colnames(xx) <- paste(colnames(xx), dn, sep="_")
    if (is.null(yy)) {
      yy <- xx
    } else {
      yy <- cbind(yy, xx)
    }
  }
  res <- write.xlsx(yy, fout, sheetName=mn, append=TRUE)
}

## Listing potential conflicts/errors
for (m in names(model.results)) {
  for (d in names(model.results[[m]])) {
    for (l in names(model.results[[m]][[d]])) {
      for (p in names(model.results[[m]][[d]][[l]]@ts)) {
        if (all(!is.finite(model.results[[m]][[d]][[l]]@ts[[p]][,'mort.frac']))) {
          print(paste("All NA:", m,d,l,p))
        } else if (any(model.results[[m]][[d]][[l]]@ts[[p]][,'mort.frac'] < 0, na.rm=TRUE)) {
          print(paste("< 0:   ", m,d,l,p))
        }
      }
    }
  }
}

df <- data.table()
for (m in names(model.results)) {
  for (d in names(model.results[[m]])) {
    for (l in names(model.results[[m]][[d]])) {
      for (p in names(model.results[[m]][[d]][[l]]@ts)) {
        df <- rbind(df, data.table(model=m, driver=d, event=l, process=sub("cmort|cmort_", "", p),
                                   mort.frac=as.vector(model.results[[m]][[d]][[l]]@ts[[p]][,'mort.frac']),
                                   mort.flux=as.vector(model.results[[m]][[d]][[l]]@ts[[p]][,'mort.flux']),
                                   cveg=as.vector(model.results[[m]][[d]][[l]]@ts[[p]][,'cveg']),
                                   Year=as.vector(time(model.results[[m]][[d]][[l]]@ts[[p]][,'mort.frac']))))
      }
    }
  }
}

## remove those locations that have no data at all
for (n in names(dm.events)) {
  if (all(!is.finite(df[event == n, ]$mort.frac)))
    df = df[event != n, ]
}
      


pdf(file="figures/histograms.pdf", width=8, height=10, paper="special")
for (l in unique(df$event)) {
  gg <- ggplot(df[event==l, ], aes(x=mort.frac, y=..density.., col=process, fill=process))
  gg <- gg + DGVMTools::dgvm.ggplot.theme("scatter")
  gg <- gg + geom_density(position="dodge", alpha=0.2)
  gg <- gg + labs(title=paste0(l, "(", dm.events[[l]]$period[1], "-", dm.events[[l]]$period[2],")"))
  gg <- gg + facet_grid(model~driver, scale="free")
  gg <- gg + scale_x_log10()
  gg <- gg + guides(fill = guide_legend(ncol = 3, title=NULL, override.aes=list(alpha=1)), col=FALSE)
  print(gg)
}
dev.off()


pdf(file="figures/timeseries.pdf", width=8, height=10, paper="special")
for (l in unique(df$event)) {
  gg <- ggplot(df[event==l, ], aes(x=Year, y=mort.frac, col=process))
  gg <- gg + DGVMTools::dgvm.ggplot.theme("temporal")
  gg <- gg + geom_line()
  gg <- gg + labs(title=paste0(l, "(", dm.events[[l]]$period[1], "-", dm.events[[l]]$period[2],")"))
  gg <- gg + facet_grid(model~driver, scale="free_y")
  gg <- gg + guides(col = guide_legend(ncol = 3, title=NULL))
  print(gg)
}
dev.off()


q()


df=data.frame()
for (m in names(model.results)) {
  for (d in names(model.results[[m]])) {
    for (l in names(model.results[[m]][[d]])) {
      df <- rbind(df, data.frame(model=m, driver=d, event=l,
                                 wilcox=model.results[[m]][[d]][[l]]@wilcox,
                                 random=sapply(model.results[[m]][[d]][[l]]@random.wilcox, mean, na.rm=TRUE),
                                 process=names(model.results[[m]][[d]][[l]]@wilcox)))
    }
  }
}

gg <- ggplot(df, aes(x=model, y=wilcox, fill=driver))
gg <- gg + geom_violin()
print(gg)
