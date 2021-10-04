#!/usr/bin/env python

import argparse
import math
import os
import re
import sys


def parse_kicad_element(el):
    parts = el.strip("()").split(' ')
    if parts[0] in ('start', 'end'):
        return parts[0], tuple((float(pp) for pp in parts[1:]))
    return parts[0]


def rotate_point(point, center, angle_deg):
    st = (point[0]-center[0], point[1]-center[1]) # translate to (0,0)
    angle_rad = angle_deg/180.0*math.pi

    # rotate
    sr = (
            st[0]*math.cos(angle_rad) + st[1]*math.sin(angle_rad),
            -st[0]*math.sin(angle_rad) + st[1]*math.cos(angle_rad)
            )

    # translate back
    return (sr[0]+center[0], sr[1]+center[1])


def process_kicad_layer(outline, arc_subdivision=1):
    # extract data
    # example:
    #   (gr_line (start 215.75 65.5) (end 226.25 65.5) (layer Edge.Cuts) (width 0.05) (tstamp 61516B52))
    #   (gr_arc (start 208.75 46) (end 212.25 46) (angle -90) (layer Edge.Cuts) (width 0.05) (tstamp 61516B4C))
    paths = []
    for ln in outline:
        p = {
                'type': ln.strip("()").split(' ')[0].split('_')[-1]
                }

        if p['type'] not in ('line', 'arc'):
            continue

        for el in re.findall(r'\([0-9a-zA-Z\. -]*\)', ln):
            el_p = el.strip("()").split(' ')
            if el_p[0] in ('start', 'end'):
                p[el_p[0]] = tuple((float(pp) for pp in el_p[1:]))
            elif el_p[0]=='angle':
                p[el_p[0]] = float(el_p[1])

        if p['type']=='arc':
            # on arc, the 'start' is the center point and 'end' is starting point
            # transform data

            # do the rotation
            p['center'] = p['start']
            p['start'] = rotate_point(point=p['end'], center=p['start'], angle_deg=-p['angle']) # angle has wrong orientation

        paths.append(p)

    return paths


def distance(p1, p2):
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

def sort_paths(paths):
    current = paths[0]
    rest = paths[1:]

    shapes = []
    current_shape = [paths[0],]
    while len(rest):
        rem = None
        for idx, p in enumerate(rest):
            # check if it starts at the end of the previous segment
            if distance(p['start'], current_shape[-1]['end'])<0.01:
                rem = idx
                break
            if distance(p['end'], current_shape[-1]['end'])<0.01:
                s = p['start']
                p['start'] = p['end']
                p['end'] = s
                p['swapped'] = True
                rem = idx
                break
        else:
            #raise ValueError("the outline isn't a closed shape")
            shapes.append(current_shape)
            current_shape = [rest.pop(), ]
            continue

        current_shape.append(rest[idx])
        del(rest[idx])

    shapes.append(current_shape)
    return shapes


def main():
    extensions = ('.stl', '.amf', '.png', '.pdf', '.scad')
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', '--shape', help="Index of the desired shape from KiCAD file", type=int, default=0)

    parser.add_argument('-l', '--margin-left', help="Left margin (mm)", type=float, default=20)
    parser.add_argument('-r', '--margin-right', help="Right margin (mm)", type=float, default=20)
    parser.add_argument('-t', '--margin-top', help="Top margin (mm)", type=float, default=20)
    parser.add_argument('-b', '--margin-bottom', help="Bottom margin (mm)", type=float, default=20)

    parser.add_argument('-m', '--mirror', help="Mirror the PCB (to get the bottom side up)", action='store_true')


    # TODO: print frame holder
    #parser.add_argument('-m', '--mirror', help="Mirror the PCB (to get the bottom side up)", action='store_true')

    parser.add_argument('-p', '--pcb-thickness', help="Thickness of the PCB (mm)", type=float, default=1.6)

    parser.add_argument('-o', '--offset', help="Offset between the PCB/stencil and frame edge (mm)", type=float, default=0.1)

    parser.add_argument('-f', '--format', help="Output format", type=str, choices=('stl', 'amf', 'png', 'scad'), default="stl")

    parser.add_argument('infile', help="path to KiCad PCB file")
    parser.add_argument('outfile', help="path to output file (extension can be %s)"%(", ".join(extensions),))

    args = parser.parse_args()

    outline_raw = []

    if not any([args.outfile.endswith(ext) for ext in extensions]):
        print("unsupported output file extension")
        return 1

    try:
        # kicad file
        # NOTE: cutouts not supported yet
        with open(args.infile, "r") as fin:
            for line in fin:
                if "(layer Edge.Cuts)" in line:
                    outline_raw.append(line.strip())
    except IOError:
        print("input file not found")

    shapes = sort_paths(process_kicad_layer(outline_raw))

    if len(shapes)>0:
        print("Found {} closed shapes inside the file".format(len(shapes)))
    else:
        print("No shapes found in the PCB file")
        return 1


    # TODO: find the outer shape, for now just take the first one
    paths = shapes[args.shape]

    # prepare list of points for OpenSCAD polygon (i.e. expand arcs)
    pol = []
    for p in paths:
        if p['type']=='arc':
            # add points across the arc
            angle_step = 1 # 1°
            angle = 0
            while abs(angle)<abs(p['angle']):
                pol.append(rotate_point(point=p['start'], center=p['center'], angle_deg=angle))
                angle += angle_step * (1 if p.get('swapped', False) else -1)
        else:
            pol.append(p['start'])

    # find center and translate
    pth_c = (sum((p[0] for p in pol))/len(pol), sum((p[1] for p in pol))/len(pol), )
    pol = [(p[0]-pth_c[0], p[1]-pth_c[1],) for p in pol]

    if not pol:
        print("no shape found on Edge.Cuts layer")
        return 1

    if args.mirror:
        # mirror around y axis
        pol = [(-p[0], p[1]) for p in pol]

    stencil_margins= [20, 20, 20, 20]

    # find polygon boundaries
    bounds = {
            'xmin': min([p[0] for p in pol]),
            'xmax': max([p[0] for p in pol]),
            'ymin': min([p[1] for p in pol]),
            'ymax': max([p[1] for p in pol]),
            }
    pol_stencil = [
            (bounds['xmin']-args.margin_left, bounds['ymin']-args.margin_top),
            (bounds['xmin']-args.margin_left, bounds['ymax']+args.margin_bottom),
            (bounds['xmax']+args.margin_right, bounds['ymax']+args.margin_bottom),
            (bounds['xmax']+args.margin_right, bounds['ymin']-args.margin_top)
            ]

    base_margin = 5
    pol_base = [
            (bounds['xmin']-args.margin_left-base_margin, bounds['ymin']-args.margin_top-base_margin),
            (bounds['xmin']-args.margin_left-base_margin, bounds['ymax']+args.margin_bottom+base_margin),
            (bounds['xmax']+args.margin_right+base_margin, bounds['ymax']+args.margin_bottom+base_margin),
            (bounds['xmax']+args.margin_right+base_margin, bounds['ymin']-args.margin_top-base_margin)
            ]

    pcb_cutout = "linear_extrude(height=10) offset(r={offset}) polygon(points={points}, convexity=10);".format(points=str([list(p) for p in pol]), offset=args.offset)

    stencil_cutout = "translate([0, 0, {thick}]) linear_extrude(height=5) offset(r={offset}) polygon(points={points});".format(points=str([list(p) for p in pol_stencil]), offset=args.offset, thick=args.pcb_thickness)

    base = "translate([0, 0, {vert}]) linear_extrude(height={height}) polygon(points={points});".format(points=str([list(p) for p in pol_base]), height=3+args.pcb_thickness, vert=-1)

    # add hole on the longest side of the PCB for easier PCB extraction
    # FIXME: probably needs another hole or two
    maxlen = 0
    maxidx = -1
    for i in range(len(pol)):
        dd = distance(pol[i], pol[(i-1)%len(pol)])
        if dd>maxlen:
            maxidx = i
            maxlen = dd
    d = min(maxlen/2, 10)
    holes = "translate([{x}, {y}, 0]) cylinder(h=20, r={r}, center=true, $fn=100);".format(x=(pol[maxidx][0]+pol[(maxidx-1)%len(pol)][0])/2, y=(pol[maxidx][1]+pol[(maxidx-1)%len(pol)][1])/2, r=d/2)

    # add holes for the stencil removal
    for i in range(4):
        holes += "translate([{x}, {y}, 0]) cylinder(h=20, r={r}, center=true, $fn=100);".format(x=(pol_base[i][0]+pol_base[(i-1)%len(pol_base)][0])/2, y=(pol_base[i][1]+pol_base[(i-1)%len(pol_base)][1])/2, r=8+base_margin)


    code = "difference(){{ \n{base} union(){{ {pcb} {stenc} {holes} }} }}".format(base=base, pcb=pcb_cutout, stenc=stencil_cutout, holes=holes)

    if args.outfile.endswith(".scad"):
        with open(args.outfile, "w") as f:
            f.write(code)
    else:
        cmd = "openscad /dev/null -D '{code}' -o {outfile}".format(code=code, outfile=args.outfile)
        os.system(cmd)

    return 0


if __name__=='__main__':
    sys.exit(main())
