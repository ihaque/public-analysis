import gzip
import re
import sys
import matplotlib.pyplot as plt

import numpy as np
from scipy.cluster.vq import kmeans2
from scipy.interpolate import interp1d

NOTDETECTED = 'style="fill:#3b4ba7;fill-opacity:1;fill-rule:nonzero;stroke:#3b4ba7;stroke-width:0.10230555;stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"'
DETECTED = 'style="fill:#ed1c24;fill-opacity:1;fill-rule:nonzero;stroke:#ed1c24;stroke-width:0.10230555;stroke-linecap:round;stroke-linejoin:round;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1"'


def get_groups(svg, eltname):
    ingroup = False
    groups = []
    buf = []

    for line in svg:
        if ingroup:
            buf.append(line)
            ingroup = not (line.endswith('/>') or line.endswith('</'+eltname+'>'))
        if not ingroup:
            if line.lstrip().startswith('<'+eltname):
                ingroup = True
                continue
            if buf:
                groups.append(buf)
                buf = []
    if buf:
        groups.append(buf)
        buf = []
    return groups

def find_circle_center(p):
    def floatify(c):
        x, y = c.split(',')
        return np.array([float(x), float(y)])

    directions = next(x for x in p if x.strip().startswith('d="'))
    directions = directions.split('"')[1].split(' ')
    assert directions[0] == 'm'
    start_point = np.array(floatify(directions[1]))
    assert directions[2] == 'c'
    deltas = np.vstack(list(map(floatify, directions[3:])))
    coordinates = start_point + np.cumsum(deltas, axis=0)
    return np.mean(coordinates, axis=0)

def savefig(fig, width, height, name):
    fig.set_figwidth(width)
    fig.set_figheight(height)
    fig.savefig(name, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    # SVG exported from 
    # https://grail.com/wp-content/uploads/ASCO_2019_CCGA1_Outcomes_Oxnard_Poster_Final.pdf
    # using Inkscape to open poster, crop plot, and ungroup all points
    with gzip.open(sys.argv[1], 'rt') as svgf:
        svg = [x.rstrip() for x in svgf]

    paths = get_groups(svg, 'path')
    tspans = get_groups(svg, 'tspan')
    
    points = []
    detected = []
    
    for p in paths:
        if any(DETECTED in x for x in p):
            points.append(find_circle_center(p))
            detected.append(True)
        if any(NOTDETECTED in x for x in p):
            points.append(find_circle_center(p))
            detected.append(False)

    points = np.array(points)
    detected = np.array(detected)
    
    # Get tick coordinates from tick marks
    svg_yticks = []
    svg_ticklabels = []
    for p in paths:
        directions = next(x for x in p if x.strip().startswith('d="'))
        matches = re.findall('m [-0-9.,]+ h ', directions)
        # There's only one path that has six horizontal lines; it's the ticks
        if len(matches) == 6:
            for m in matches:
                dx, dy = m.split(' ')[1].split(',')
                svg_yticks.append(float(dy))
    # Coordinates in path were relative moves; sum them to get absolute cords
    svg_yticks = np.cumsum(svg_yticks)
    
    # Now find the tick labels and match them up with the ticks
    for tspan in tspans:
        text = re.search('[^>]+>([^<]+)<.*', tspan[-1]).groups()[0]
        if not text.endswith('%'):
            continue
        ycor = next(l for l in tspan if l.strip().startswith('y='))
        ycor = float(ycor.split('"')[1])
        svg_ticklabels.append((text, ycor)) 
    
    # Remap label values to closest tick coordinate
    yticks = []
    ytickvals = []
    for val, y in svg_ticklabels:
        ytickvals.append(0.01 * float(val.rstrip('%')))
        yticks.append(svg_yticks[np.argmin(np.abs(svg_yticks - y))])
    yticks = np.array(yticks)
    ytickvals = np.array(ytickvals)
    
    # Interpolate on log scale to get TF from plot coordinates
    y2log10tf = interp1d(yticks, np.log10(ytickvals), fill_value='extrapolate')
    tf = np.power(10, y2log10tf(points[:,1]))
    
    # Find stage assignments by clustering along x
    centroids, labels = kmeans2(points[:,0], 4, iter=100,
                                missing='raise')
    sorted_cen = sorted(centroids)
    centroids = list(centroids)
    label2stage = {centroids.index(sorted_cen[0]): 1,
                   centroids.index(sorted_cen[1]): 2,
                   centroids.index(sorted_cen[2]): 3,
                   centroids.index(sorted_cen[3]): 4}
    stages = np.array([label2stage[l] for l in labels])

    with open('points.csv', 'wt') as outf:
        outf.write('Stage,Tumor Fraction,Detected by CCGA1 WGBS assay?\n')
        for i in range(len(points)):
            outf.write('%d,%.6f,%s\n' % (stages[i], tf[i], detected[i]))
    
    fig1 = plt.figure()
    # Invert y for plotting because SVG coordinate system is flipped
    plt.scatter(points[~detected,0], -points[~detected,1],
                color='b', label='Not Detected',
                alpha=0.75)
    plt.scatter(points[detected,0], -points[detected,1],
                color='r', label='Detected',
                alpha=0.75)
    plt.yticks(-yticks, ['%.3f%%' % (t * 100) for t in ytickvals])
    plt.xticks(sorted_cen,['I', 'II', 'III', 'IV'])
    plt.xlabel('Stage')
    plt.ylabel('Tumor fraction')
    plt.legend()
    savefig(fig1, 6, 6, 'asco-2019-oxnard-fig3-reproduced.png')
    
    def make_tf_detection_plot(predicate=None, weights=None):
        if predicate is None:
            predicate = np.ones(tf.shape, dtype=np.bool)
        if weights is None:
            weights = np.ones(tf.shape, dtype=np.float)
        plt.hist(
            [tf[predicate & detected], tf[predicate & ~detected]],
            bins=np.logspace(-5,0, 100),
            cumulative=True, density=True, stacked=True,
            weights=[weights[predicate & detected],
                     weights[predicate & ~detected]],
            color='rb',
            label=['detected', 'not detected'])
        plt.gca().set_xscale('log')
        plt.xlabel('Tumor Fraction')
        plt.ylabel('Cumulative prop. of samples')
        plt.legend()
    
    fig2 = plt.figure()
    for stage in range(1,5):
        plt.subplot(2,2,stage)
        make_tf_detection_plot(stages == stage)
        plt.title('Stage %d (N=%d)' %
                  (stage, sum((stages == stage).astype(np.int))))
    plt.suptitle('Cumulative distribution of tumor fraction by stage')
    savefig(fig2, 10, 10, 'oxnard-tf-detection-by-stage.png')

    # Combine stages I-III and IV
    fig3 = plt.figure()
    plt.subplot(121)
    stage123 = (stages == 1) | (stages == 2) | (stages == 3)
    make_tf_detection_plot(stage123)
    plt.title('Stages I-III (N=%d)' % (sum((stage123).astype(np.int))))
    plt.subplot(122)
    make_tf_detection_plot(stages == 4)
    plt.title('Stage IV (N=%d)' % (sum((stages == 4).astype(np.int))))
    savefig(fig3, 12, 6, 'oxnard-tf-detection-stage13-4.png')

    # Number of CRCs detected by stage in 10k screening study in Imperiale 2014
    stage_weights = np.array([np.nan, 29., 21., 10., 4.])
    sample_weights = stage_weights[stages]
    fig4 = plt.figure()
    make_tf_detection_plot(predicate=None, weights=sample_weights)
    plt.gca().set_xscale('log')
    plt.legend()
    plt.suptitle('Cumulative distribution of tumor fraction, weighted by stage incidence in screening')
    savefig(fig4, 8, 8, 'oxnard-tf-detection-incidence-weighted.png')
    plt.show()
