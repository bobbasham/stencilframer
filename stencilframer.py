#!/usr/bin/env python
# -*- coding: utf-8 -*-

# The MIT License (MIT)
#
# Copyright © 2021 Igor Brkic <i@hglt.ch>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the “Software”), to deal in
# the Software without restriction, including without limitation the rights to use,
# copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the
# Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import argparse
import enum
import logging
import math
import os
import re
import sys
import tempfile


class InFormat(enum.Enum):
    KICAD = 0
    GERBER = 1

class Interpolation(enum.Enum):
    LINEAR = 0
    ARC_CW = 1
    ARC_CCW = 2

class CoordFormat(enum.Enum):
    ABSOLUTE = 0
    INCREMENTAL = 1


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


def distance(p1, p2):
    # euclidean distance
    return math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)


def get_angle(s, e, c):
    arg = distance(s, e)/2 / distance(s, c)
    return math.asin(arg) * 2 / math.pi * 180


def parse_kicad_element(el):
    parts = el.strip("()").split(' ')
    if parts[0] in ('start', 'end'):
        return parts[0], tuple((float(pp) for pp in parts[1:]))
    return parts[0]


def process_kicad_layer(outline, arc_subdivision=1):
    # extract data
    # example:
    #   (gr_line (start 215.75 65.5) (end 226.25 65.5) (layer Edge.Cuts) (width 0.05) (tstamp 61516B52))
    #   (gr_arc (start 208.75 46) (end 212.25 46) (angle -90) (layer Edge.Cuts) (width 0.05) (tstamp 61516B4C))
    paths = []
    for ln in outline:
        if "(layer Edge.Cuts)" not in ln:
            # extract only the Edge.Cuts layer
            continue

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
            # FIXME: revisit arc orientation
            # on arc, the 'start' is the center point and 'end' is starting point
            # transform data

            # do the rotation
            p['center'] = p['start']
            p['start'] = rotate_point(point=p['end'], center=p['start'], angle_deg=-p['angle']) # angle has wrong orientation
            logging.debug("interpolating arc from %s to %s", p['start'], p['end'])
        else:
            logging.debug("interpolating line from %s to %s", p['start'], p['end'])

        paths.append(p)

    return paths


def process_gerber_layer(outline, arc_subdivision=1):
    apertures = {}
    integers = 4
    decimals = 6
    unit_convert = 1
    total_coord = 6
    interpolation = Interpolation.LINEAR
    point = (0,0)
    coord_format = CoordFormat.ABSOLUTE

    paths = []

    # convert raw coordinate to float in mm
    coord = lambda x: float(x)/(10**decimals)/unit_convert

    lineno = 0
    ln = outline[lineno]
    while lineno<len(outline):
        if ln.startswith('%FSLA'):
            # parse decimal format
            try:
                fmts = re.findall(r'FSLAX([0-9]+)Y([0-9]+)', ln)[0]
            except IndexError:
                raise ValueError("Invalid coordinate format specified")
            if len(fmts)!=2 or fmts[0]!=fmts[1]:
                raise ValueError("Invalid coordinate format specified")
            integers = int(fmts[0][0])
            decimals = int(fmts[0][1:])
            logging.debug("coordinate format set to %d.%d", decimals, total_coord)

        elif ln.startswith('%MOIN') or ln.startswith('G70*'):
            # units in inches - convert to mm
            unit_convert = 25.4
            logging.debug("units set to inches")

        elif ln.startswith('%MOMM') or ln.startswith('G71*'):
            # units in inches - convert to mm
            unit_convert = 1
            logging.debug("units set to mm")

        elif ln.startswith('%ADD'):
            # apertures list
            # since we're only looking for outline we can assume 0.05mm
            # aperture size and skip this for now
            pass

        elif ln.startswith('G01'):
            interpolation = Interpolation.LINEAR
            logging.debug("interpolation set to %s", interpolation)
            ln = ln[3:] # continue parsing the same line
            continue

        elif ln.startswith('G02'):
            interpolation = Interpolation.ARC_CW
            logging.debug("interpolation set to %s", interpolation)
            ln = ln[3:] # continue parsing the same line
            continue

        elif ln.startswith('G03'):
            interpolation = Interpolation.ARC_CCW
            logging.debug("interpolation set to %s", interpolation)
            ln = ln[3:] # continue parsing the same line
            continue

        elif ln.startswith('G04'):
            # just a comment, skip it
            pass

        elif ln.startswith('G75'):
            # multi quadrant mode (legacy)
            pass

        elif ln.startswith('G90'):
            # set the coordinates to absolute (legacy)
            coord_format = CoordFormat.ABSOLUTE
            logging.debug("coordinate format set to absolute")
            ln = ln[3:] # continue parsing the same line
            continue
        elif ln.startswith('G91'):
            # set the coordinates to absolute (legacy)
            coord_format = CoordFormat.INCREMENTAL
            logging.debug("coordinate format set to incremental")
            ln = ln[3:] # continue parsing the same line
            continue

        elif ln.startswith('D'):
            logging.debug("selecting aperture %s (ignored)", ln)
            # select apperture from the list - skip it (for now)
            pass

        elif ln.startswith('X') or ln.startswith('Y'):
            # move or interpolate command
            x = None
            y = None
            i = None
            j = None
            success = False

            # FIXME: take into consideration absolute/incremental format
            try:
                pts = re.findall(r'X([0-9\-]+)Y([0-9\-]+)I([0-9\-]+)J([0-9\-]+)', ln)[0]
                x = coord(pts[0])
                y = coord(pts[1])
                i = coord(pts[2])
                j = coord(pts[3])
                success = True
            except IndexError:
                pass

            if not success:
                try:
                    pts = re.findall(r'X([0-9\-]+)Y([0-9\-]+)', ln)[0]
                    x = coord(pts[0])
                    y = coord(pts[1])
                    success = True
                except IndexError:
                    pass

            if not success:
                try:
                    pts = re.findall(r'X([0-9\-]+)', ln)[0]
                    x = coord(pts)
                    y = point[1]
                    success = True
                except IndexError:
                    pass

            if not success:
                try:
                    pts = re.findall(r'Y([0-9\-]+)', ln)[0]
                    x = point[0]
                    y = coord(pts)
                    success = True
                except IndexError:
                    pass

            if not success:
                logging.warning("Failed to parse coordinates: %s", ln)

            center = None
            if i is not None:
                center = (point[0]+i, point[1]+j)

            if ln.endswith('D02*') or ln.endswith('D2*'):
                # move command
                point = (x, y)
                logging.debug("moving to %s", point)

            elif ln.endswith('D01*') or ln.endswith('D1*'):
                # interpolate command
                pend = (x, y)
                if interpolation==Interpolation.LINEAR:
                    logging.debug("interpolating line from %s to %s", point, pend)
                    paths.append({
                        'type': 'line',
                        'start': point,
                        'end': pend
                        })
                elif interpolation==Interpolation.ARC_CW:
                    logging.debug("interpolating CW arc from %s to %s with center at %s", point, pend, center)
                    paths.append({
                        'type': 'arc',
                        'start': pend,
                        'end': point,
                        'center': center,
                        'angle': get_angle(point, pend, center)
                        })
                elif interpolation==Interpolation.ARC_CCW:
                    logging.debug("interpolating CCW arc from %s to %s with center at %s", point, pend, center)
                    paths.append({
                        'type': 'arc',
                        'start': point,
                        'end': pend,
                        'center': center,
                        'angle': get_angle(point, pend, center)
                        })
                point = pend

            else:
                logging.warning("only supported commands currently are move or interpolate. Continuing...")

        elif ln.startswith('M02*'):
            # end of file
            break

        lineno += 1     # go to the next line
        ln = outline[lineno]

    return paths


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

    # stencil margins
    parser.add_argument('-l', '--margin-left', help="Left margin (mm)", type=float, default=20)
    parser.add_argument('-r', '--margin-right', help="Right margin (mm)", type=float, default=20)
    parser.add_argument('-t', '--margin-top', help="Top margin (mm)", type=float, default=20)
    parser.add_argument('-b', '--margin-bottom', help="Bottom margin (mm)", type=float, default=20)

    # pcb properties
    parser.add_argument('-m', '--mirror', help="Mirror the PCB (to get the bottom side up)", action='store_true')
    parser.add_argument('-p', '--pcb-thickness', help="Thickness of the PCB (mm)", type=float, default=1.6)
    parser.add_argument('-s', '--shape', help="Index of the desired shape from input file", type=int, default=0)

    # 3D model specifics
    parser.add_argument('-f', '--frame', help="Generate stencil holding frame instead of stencil frame", action='store_true')
    parser.add_argument('-k', '--skip-holes', help="Don't add holes for easy removal in the fixture", action='store_true')
    parser.add_argument('-o', '--offset', help="Offset between the PCB/stencil and frame edge (mm)", type=float, default=0.1)
    parser.add_argument('--stencil-offset', help="Offset between the stencil and frame edge (mm). If not specified, the --offset is used", type=float, default=None)

    parser.add_argument('-d', '--debug', help="Show debugging info", action='store_true')
    parser.add_argument('-w', '--use-temp-file', help="Use temporary file when calling OpenSCAD (used by default on Windows)", action='store_true')

    parser.add_argument('--openscad', help="Path to OpenSCAD executable", type=str, default="openscad")
    parser.add_argument('infile', help="path to KiCad PCB or gerber file (.kicad_pcb, .gbr, .gm1)")
    parser.add_argument('outfile', help="path to output file (extension can be %s)"%(", ".join(extensions),))

    args = parser.parse_args()


    logfmt = '[%(levelname)s] %(message)s'
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format=logfmt)
    else:
        logging.basicConfig(level=logging.INFO, format=logfmt)

    outline_raw = []

    if not any([args.outfile.endswith(ext) for ext in extensions]):
        logging.warning("unsupported output file extension")
        return 1

    if os.system(args.openscad+" -v")!=0:
        logging.warning("OpenSCAD executable not found on the system.")
        return 1

    if args.stencil_offset is None:
        args.stencil_offset = args.offset

    if args.infile.lower().endswith('.kicad_pcb'):
        informat = InFormat.KICAD
    elif args.infile.lower()[-4:] in ('.gbr', '.gm1'):
        informat = InFormat.GERBER
    else:
        logging.warning("invalid input file format")
        return 1

    try:
        with open(args.infile, "r") as fin:
            for line in fin:
                outline_raw.append(line.strip())
    except IOError:
        logging.warning("input file not found")
        return 1

    if informat==InFormat.KICAD:
        shapes = sort_paths(process_kicad_layer(outline_raw))
    elif informat==InFormat.GERBER:
        shapes = sort_paths(process_gerber_layer(outline_raw))

    if len(shapes)>0:
        logging.info("Found %d closed shapes inside the file", len(shapes))
    else:
        logging.warning("No shapes found in the PCB file")
        return 1


    # TODO: find the outer shape (the one with greatest surface, for now just take the first one
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
        logging.warning("no shape found on Edge.Cuts layer")
        return 1

    # since the KiCAD and OpenSCAD use different y-axis direction, the polygon is mirrored
    # by default
    if informat==InFormat.KICAD:
        pol = [(p[0], -p[1]) for p in pol]

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

    # OpenSCAD code generation start

    code = ""

    if args.frame:
        # generate just the frame to hold the stencil in place
        stencil_bounds = {
                'xmin': min([p[0] for p in pol_stencil]),
                'xmax': max([p[0] for p in pol_stencil]),
                'ymin': min([p[1] for p in pol_stencil]),
                'ymax': max([p[1] for p in pol_stencil]),
                }
        # add chamfer to frame
        chamf = min([(stencil_bounds['xmax']-stencil_bounds['xmin'])/5, (stencil_bounds['ymax']-stencil_bounds['ymin'])/5, 15])
        pol_frame_out = [
            (stencil_bounds['xmin']+chamf, stencil_bounds['ymin']),
            (stencil_bounds['xmax']-chamf, stencil_bounds['ymin']),

            (stencil_bounds['xmax'], stencil_bounds['ymin']+chamf),
            (stencil_bounds['xmax'], stencil_bounds['ymax']-chamf),

            (stencil_bounds['xmax']-chamf, stencil_bounds['ymax']),
            (stencil_bounds['xmin']+chamf, stencil_bounds['ymax']),

            (stencil_bounds['xmin'], stencil_bounds['ymax']-chamf),
            (stencil_bounds['xmin'], stencil_bounds['ymin']+chamf),
            ]

        frame_out = "linear_extrude(height=5) offset(r={offset}) polygon(points={points}, convexity=10);".format(points=str([list(p) for p in pol_frame_out]), offset=args.stencil_offset)
        frame_in = "translate([0, 0, -2]) linear_extrude(height=10) offset(r={offset}) polygon(points={points}, convexity=10);".format(points=str([list(p) for p in pol_frame_out]), offset=-5)

        code="difference(){{ {fout} {fin} }}".format(fout=frame_out, fin=frame_in)

    else:
        # generate the actual stencil frame

        # arrange cutouts for PCB and stencil
        pcb_cutout = "linear_extrude(height=10) offset(r={offset}) polygon(points={points}, convexity=10);".format(points=str([list(p) for p in pol]), offset=args.offset)

        stencil_cutout = "translate([0, 0, {thick}]) linear_extrude(height=5) offset(r={offset}) polygon(points={points});".format(points=str([list(p) for p in pol_stencil]), offset=args.stencil_offset, thick=args.pcb_thickness)

        base = "translate([0, 0, {vert}]) linear_extrude(height={height}) polygon(points={points});".format(points=str([list(p) for p in pol_base]), height=3+args.pcb_thickness, vert=-1)

        holes = ""
        if not args.skip_holes:
            # add hole on the longest side of the PCB for easier PCB extraction
            # FIXME: probably needs another hole or two and better placement
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
                holes += "translate([{x}, {y}, 0]) cylinder(h=20, r={r}, center=true, $fn=100);".format(x=(pol_base[i][0]+pol_base[(i-1)%len(pol_base)][0])/2, y=(pol_base[i][1]+pol_base[(i-1)%len(pol_base)][1])/2, r=7+base_margin)

        code = "difference(){{ \n{base} union(){{ {pcb} {stenc} {holes} }} }}".format(base=base, pcb=pcb_cutout, stenc=stencil_cutout, holes=holes)

    # produce the output file
    if args.outfile.endswith(".scad"):
        with open(args.outfile, "w") as f:
            f.write(code)
    else:
        if args.use_temp_file or os.name!='posix':
            # Windows will report "command line too long" if all code is passed as a command line argument
            # that's why we'll use temporary file
            fd, name = tempfile.mkstemp()
            os.write(fd, code.encode('utf-8'))
            os.close(fd)
            cmd = "openscad {infile} -o {outfile}".format(infile=name, outfile=args.outfile)
            os.system(cmd)
            os.unlink(name)

        else:
            cmd = "openscad /dev/null -D '{code}' -o {outfile}".format(code=code, outfile=args.outfile)
            os.system(cmd)

    return 0


if __name__=='__main__':
    sys.exit(main())
