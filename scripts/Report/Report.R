#+setup, include=FALSE
library(knitr)
library(ggplot2)

header <- "
\\newenvironment{knitrout}{}{} % an empty environment to be redefined in TeX
\\definecolor{fgcolor}{rgb}{0.345, 0.345, 0.345}

\\usepackage{geometry}
\\geometry{top=20mm}
\\geometry{bottom=25mm}
\\geometry{left=20mm}
\\geometry{right=20mm}

\\usepackage{titlesec}
\\titleformat
{\\section} % command
[display] % shape
{\\bf\\LARGE} % format
{} % label
{0.1ex} % sep
{ } % before-code \\thesection.
[] % after-code

\\titlespacing{\\section}{0pt}{*0.5}{*1}
"

#opts_knit$set(header=header)
opts_chunk$set(fig.path = 'figures/f-', echo=FALSE, out.width='.8\\linewidth')

#\fontsize{14}{2}\selectfont

theme_set(theme_gray(base_size = 18))

#+Report
#\section{Merge}
probs_hist <- qplot(max_merge_probs, bins=100, xlab=merge_xlabel, ylab='Number of cells', col=I("black"))
if (probs_threshold >= 0) {
    probs_hist <- probs_hist + geom_vline(aes(xintercept=probs_threshold, color='Optimal'), show.legend=TRUE) +
    scale_color_manual(name = "Threshold", values = c(Optimal = "red")) + theme(legend.position=c(0.8, 0.9))
}
print(probs_hist)

qplot(nonzero_neighbours_num, bins=10, xlab='Number of the neighboues with nonzero intersection', ylab='#Cells')