# Analysis code accompanying early detection roundup, 2019-07-01

This directory accompanies the post [Mid-2019 Early Detection Roundup, Part I](https://ihaque.org/posts/2019/07/02/early-detection-mid-2019-part-1/).

In support of the analysis in the post, I extracted figure 3 from [Oxnard et al.](http://web.archive.org/web/20190629031324/https://grail.com/wp-content/uploads/ASCO_2019_CCGA1_Outcomes_Oxnard_Poster_Final.pdf) using Inkscape and ungrouped the elements, producing drawing-ungrouped.svg.gz in this directory.

`parse_oxnard_svg.py` takes this file as a command line argument, extracts the data points in the figure, and produces the plots in the referenced blog post.

