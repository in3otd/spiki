#!/usr/bin/env python

# Licensed under GNU GPL version 2 or later

import math
import re
import subprocess

# ideas from
# https://wiki.python.org/moin/PointsAndRectangles
# and
# http://stackoverflow.com/questions/12468900/making-a-point-class-in-python


class Point(object):
    """ Point class represents and manipulates x,y coords. """

    def __init__(self, x=0, y=0):
        """ Create a new point at x, y """
        self.x = x
        self.y = y

    def length(self):
        """Calculate length of vector to point from origin."""
        return math.hypot(self.x, self.y)

    def dist(self, p):
        return (self - p).length()

    def __add__(self, p):
        """Point(x1+x2, y1+y2)"""
        return Point(self.x + p.x, self.y + p.y)

    def __sub__(self, p):
        """Point(x1-x2, y1-y2)"""
        return Point(self.x - p.x, self.y - p.y)

    def __str__(self):
        return "(%s, %s)" % (self.x, self.y)

    def __repr__(self):
        return "%s(%r, %r)" % (self.__class__.__name__, self.x, self.y)

    def copy(self):
        """Return a full copy of this point."""
        return Point(self.x, self.y)

    def slide(self, p):
        '''Move to new (x+dx,y+dy).

        Can anyone think up a better name for this function?
        slide? shift? delta? move_by?
        '''
        self.x = self.x + p.x
        self.y = self.y + p.y

    def slide_xy(self, dx, dy):
        '''Move to new (x+dx,y+dy).
        Can anyone think up a better name for this function?
        slide? shift? delta? move_by?
        '''
        self.x = self.x + dx
        self.y = self.y + dy

    def rotate(self, rad):
        """Rotate counter-clockwise by rad radians.

        Positive y goes *up,* as in traditional mathematics.

        Interestingly, you can use this in y-down computer graphics, if
        you just remember that it turns clockwise, rather than
        counter-clockwise.

        """
        s, c = [f(rad) for f in (math.sin, math.cos)]
        self.x, self.y = (c * self.x - s * self.y, s * self.x + c * self.y)

    def rotate_about(self, p, theta):
        """Rotate counter-clockwise around a point, by theta degrees.

        Positive y goes *up,* as in traditional mathematics.

        """
        self.slide_xy(-p.x, -p.y)
        self.rotate(theta)
        self.slide_xy(p.x, p.y)


class kmodule():
    """KiCad module."""

    def __init__(self, fname):
        self.fout = open(fname[0], 'w')

    def write_header(self, name='SIND', descr='spiral inductor', tags='SMD'):
        self.fout.write('(module %s (layer F.Cu)\n' % name)
        self.fout.write('  (at 0 0)\n')
        self.fout.write('  (descr "%s")\n' % descr)
        self.fout.write('  (tags "%s")\n' % tags)
        self.fout.write('  (attr smd)\n')

    def write_refs(self, x, y, ref='REF**', value='LLL'):
        self.fout.write(
            '  (fp_text reference %s (at %f %f) (layer F.SilkS)\n' %
            (ref, x, y))
        self.fout.write('    (effects (font (size 1 1) (thickness 0.15)))\n')
        self.fout.write('  )\n')
        self.fout.write(
            '  (fp_text value %s (at %f %f) (layer F.Fab)\n' %
            (value, x, y + 1))
        self.fout.write('    (effects (font (size 1 1) (thickness 0.15)))\n')
        self.fout.write('  )\n')

    def add_circ_spiral(self, vertices, layer, width):
        start = vertices[0]
        for v in vertices[1:]:
            self.add_line(start, v, layer, width)
            start = v.copy()

    def add_arc(self, centre, start, angle, layer, pen_width):
        """
        centre -- the position of the centre of the circle that is the basis of the arc
        start -- the starting point of the arc. Both the radius of the arc and the start angle are calculated from this point.
        angle -- the angular span that the arc covers in radians, from the start angle, in clock-wise direction
        layer -- the name of the layer for the line
        pen_width -- the pen width
        """
        # Kicad arc definition is a bit weird:
        #   start is the arc centre
        #   end is the arc starting point
        #   angle is the angular span, as expected
        self.fout.write(
            '  (fp_arc (start %f %f) (end %f %f) (angle %f) (layer %s) (width %f))\n' %
            (centre.x, centre.y, start.x, start.y, math.degrees(angle), layer, pen_width))

    def add_arc_spiral(self, arcs, layer, width):
        for a in arcs:
            centre, start, angle = a
            self.add_arc(centre, start, angle, layer, width)

    def add_line(self, start, end, layer, width):
        """Draw line"""
        # (fp_line (start x y ) (end x y) (layer name) (width pen )
        self.fout.write(
            '  (fp_line (start %f %f) (end %f %f) (layer %s) (width %f))\n' %
            (start.x, start.y, end.x, end.y, layer, width))

    def add_smd_pad(self, name, shape, origin, size, layer):
        self.fout.write(
            '  (pad %s smd %s (at %f %f) (size  %f %f) (layers %s.Cu %s.Paste %s.Mask))\n' %
            (str(name), shape, origin.x, origin.y, size.x, size.y, layer, layer, layer))

    def add_thru_pad(self, name, shape, origin, size, drill):
        self.fout.write(
            '  (pad %s thru_hole %s (at %f %f) (size %f %f) (drill %f) (layers *.Cu *.Mask F.SilkS))\n' %
            (str(name), shape, origin.x, origin.y, size.x, size.y, drill))

    def close(self):
        self.fout.write(')\n')
        self.fout.close()


class fh_file():
    """Fasthenry input file."""

    def __init__(self, fname):
        self.fname = fname
        self.fout = open(self.fname, 'w')
        self.p = None  # fasthenry process
        self.N = 5  # number of straight lines per arc
        self.last_node = 0
        self.arc_idx = 0
        self.spiral_idx = 1

    def write_header(self, nwinc=1, nhinc=1):
        self.fout.write('* Spiral coil\n')
        self.fout.write('\n')
        self.fout.write('.Units mm\n')
        # assume copper conductors
        self.fout.write(
            '.default sigma=5.8e4 nwinc=%i nhinc=%i\n' %
            (nwinc, nhinc))
        self.fout.write('\n')

    def add_arc(self, centre, start, angle, layer, width):
        if (layer == 'F.Cu'):
            z = 0.0
        else:
            z = -PCB_h * 1e3

        self.fout.write('* arc #%i\n' % self.arc_idx)
        self.fout.write('* %s %s %f\n' % (centre, start, angle))
        P1 = start
        for i in range(1, self.N + 1):
            P2 = start.copy()
            P2.rotate_about(centre, i * angle / float(self.N))

            if (self.last_node == 0):
                self.fout.write(
                    'N%i_%i x=%f y=%f z=%f\n' %
                    (self.spiral_idx, self.last_node, P1.x, P1.y, z))
            self.last_node = self.last_node + 1
            self.fout.write('N%i_%i x=%f y=%f z=%f\n' %
                            (self.spiral_idx, self.last_node, P2.x, P2.y, z))
            self.fout.write(
                'E%i_%i_%i N%i_%i N%i_%i w=%f h=%f\n' %
                (self.spiral_idx,
                 self.arc_idx,
                 i,
                 self.spiral_idx,
                 self.last_node -
                 1,
                 self.spiral_idx,
                 self.last_node,
                 width,
                 Cu_t *
                 1e3))
        self.arc_idx = self.arc_idx + 1

    def add_circ_spiral(
            self,
            vertices,
            layer,
            width,
            Cu_thickness,
            PCB_thickness):
        z = -(layer - 1.0) * PCB_thickness

        start = vertices[0]
        self.fout.write('N%i_%i x=%f y=%f z=%f\n' %
                        (self.spiral_idx, self.last_node, start.x, start.y, z))
        for v in vertices[1:]:
            self.last_node = self.last_node + 1
            self.fout.write('N%i_%i x=%f y=%f z=%f\n' %
                            (self.spiral_idx, self.last_node, v.x, v.y, z))
            self.fout.write(
                'E%i_%i N%i_%i N%i_%i w=%f h=%f\n' %
                (self.spiral_idx,
                 self.last_node,
                 self.spiral_idx,
                 self.last_node -
                 1,
                 self.spiral_idx,
                 self.last_node,
                 width,
                 Cu_thickness))

    def add_line(self, start, end, layer, width):
        if (layer == 'F.Cu'):
            z = 0.0
        else:
            z = -PCB_h * 1e3
        self.fout.write('* segment\n')
        if (self.last_node == 0):
            self.fout.write(
                'N%i_%i x=%f y=%f z=%f\n' %
                (self.spiral_idx, self.last_node, start.x, start.y, z))
        self.last_node = self.last_node + 1
        self.fout.write('N%i_%i x=%f y=%f z=%f\n' %
                        (self.spiral_idx, self.last_node, end.x, end.y, z))
        self.fout.write(
            'E%i_%i N%i_%i N%i_%i w=%f h=%f\n' %
            (self.spiral_idx,
             self.last_node,
             self.spiral_idx,
             self.last_node - 1,
             self.spiral_idx,
             self.last_node,
             width,
             Cu_t * 1e3))

    def add_thru_pad(self, name, shape, origin, size, drill):
        self.last_node = self.last_node + 1
        self.fout.write('* thru pad\n')
        self.fout.write('N%i x=%f y=%f z=%f\n' %
                        (self.last_node, origin.x, origin.y, 0.0))
        self.fout.write('N%i x=%f y=%f z=%f\n' %
                        (self.last_node, origin.x, origin.y, PCB_h * 1e3))

    def add_ports(self):
        self.fout.write('\n')
        self.fout.write('* ports\n')
        self.fout.write('.external N%i_%i N%i_%i\n' %
                        (self.spiral_idx, 0, self.spiral_idx, self.last_node))
        self.fout.write('\n')
        self.spiral_idx = self.spiral_idx + 1
        self.last_node = 0
        self.arc_idx = 0

    def add_frequency(self, freqs, npts=1):
        self.fout.write('* analysis frequency range\n')
        if (not isinstance(freqs, list)):
            self.fout.write(
                '.freq fmin=%e fmax=%e ndec=%i\n' %
                (freqs, freqs, npts))
        self.fout.write('\n')

    def close(self):
        self.fout.write('* The end\n')
        self.fout.write('.end\n')
        self.fout.close()

    def run(self):
        cmd = ['fasthenry', self.fname]
        logfile = 'fasthenry.log'
        with open(logfile, 'w') as flog:
            self.p = subprocess.Popen(
                cmd, shell=False, universal_newlines=True, stdout=flog)
            self.p.wait()  # should add timeout

    @staticmethod
    def readZc():
        pat = re.compile(r'^Impedance matrix for frequency = (\d*\.?\d+e?[+-]?\d+) (\d+) x (\d+)')
        freqs = []
        mats = []
        with open("Zc.mat") as fzc:
            for lin in fzc:
                match = pat.match(lin)
                if match:
                    mat = []
                    freq = float(match.group(1))
                    freqs.append(freq)
                    nrows = int(match.group(2))
                    ncols = int(match.group(3))
                    # take real and imag parts
                    pair = re.compile(r' *(\S+) +(\S+)j')
                    for ridx in range(nrows):
                        lin = fzc.readline()
                        # parse all complex data pairs
                        celms = pair.findall(lin)
                        matr = []
                        for cidx in range(ncols):
                            Zij = complex(*map(float, celms[cidx]))
                            matr.append(Zij)
                        mat.append(matr)
                    mats.append(mat)
        return (freqs, mats)


def calc_ind(n, dout, din):
    """
    formula from "Simple  Accurate Expressions for Planar Spiral Inductances", S. S. Mohan, M. del Mar Hershenson, S. P. Boyd and T. H. Lee
    """
    davg = 0.5 * (dout + din)  # average diameter
    rho = (dout - din) / (dout + din)  # fill ratio
    # coefficients for the current sheet expression
    c = [1.0, 2.46, 0.0, 0.2]  # "circular" spiral
    mu0 = 4e-7 * math.pi  # vacuum permeability
    Lgmd = 0.5 * mu0 * n**2 * davg * \
        c[0] * (math.log(c[1] / rho) + c[2] * rho + c[3] * rho**2)
    return Lgmd


def calc_mut(n, x):
    """
    Compute the coupling coefficient between two spiral inductors.

    Keyword arguments:
    n -- number of turns of each spiral inductor
    x -- distance between the spiral inductors (PCB thickness)

    Formula from "A new calculation for designing multilayer planar spiral inductors", J. Zhao
    """
    A = 0.184
    B = -0.525
    C = 1.038
    D = 1.001
    x = 1e3 * x  # in the Zhao formula the distance is in mm...
    Kc = (n**2) / ((A * x**3 + B * x**2 + C * x + D)
                   * ((1.67 * n**2 - 5.84 * n + 65) * 0.64))
    return Kc


def draw_circ_spiral(N_turns, r_in, pitch, tr_w, dir, d=0.1):
    """
     Draw a circular spiral approximation using straight segments

    Keyword arguments:
    N_turns -- number of turns
    r_in -- internal radius
    pitch -- spacing between conductors centers
    tr_w -- trace width
    N -- number of generator polygon vertices
    dir -- spiral direction/PCB side
    d -- max error w.r.t.circular arc
    """

    layer = 'F.Cu' if (dir == 1) else 'B.Cu'
    p_start = Point(r_in, 0)
    theta = 0.0
    done = False
    r = r_in
    while not done:
        #arclen = math.pi / 10.0
        # arc length to have a max error of 'd' w.r.t. a circular arc
        #   works well enough for a spiral arc, if r is not too small
        arclen = 2.0 * math.acos(1.0 - d / r)
        theta = theta + arclen
        if (theta > (2.0 * math.pi * N_turns)):
            theta = 2.0 * math.pi * N_turns
            done = True
        r = r_in + pitch * theta / (2.0 * math.pi)
        p_end = Point(r * math.cos(theta), r * math.sin(theta))
        sm.add_line(p_start, p_end, layer, tr_w)
        sf.add_line(p_start, p_end, layer, tr_w)
        p_start = p_end.copy()


def circ_spiral(N_turns, r_in, pitch, dir, d=0.1):
    """
    Draw a circular spiral approximation using straight segments

    Keyword arguments:
    N_turns -- number of turns
    r_in -- internal radius
    pitch -- spacing between conductors centers
    dir -- spiral direction
    d -- max error w.r.t.circular arc

    Returns segments vertices list
    """

    p_start = Point(r_in, 0)
    vertices = [p_start]  # first point
    theta = 0.0
    done = False
    r = r_in
    while not done:
        # arclen = math.pi / 10.0 # fixed angle step
        # arc length to have a max error of 'd' w.r.t. a circular arc
        #   works well enough for a spiral arc, if r is not too small
        arclen = 2.0 * math.acos(1.0 - d / r)
        theta = theta + arclen
        if (theta > (2.0 * math.pi * N_turns)):
            theta = 2.0 * math.pi * N_turns
            done = True
        r = r_in + pitch * theta / (2.0 * math.pi)
        p_end = Point(r * math.cos(dir * theta), r * math.sin(dir * theta))
        vertices.append(p_end)

    return vertices


def draw_arcs_spiral(N_turns, r_in, pitch, tr_w, N, dir):
    """
    Draw a spiral approximation using circular arcs.

    Keyword arguments:
    N_turns -- number of turns
    r_in -- internal radius
    pitch -- spacing between conductors centers
    tr_w -- trace width
    N -- number of generator polygon vertices
    dir -- spiral direction/PCB side

    see
    "Scan Converting Spirals", F. Taponecco and M. Alexa,
    http://wscg.zcu.cz/wscg2002/Papers_2002/F11.pdf
    and
    "Piecewise Circular Approximation of Spirals and Polar Polynomials", F. Taponecco and M. Alexa,
    http://147.228.63.9/wscg2003/papers_2003/i31.pdf
    """

    theta = 2 * math.pi / N  # arc length in radians; polygon central angle
    # polygon radius (distance from the center to a vertex)
    p = pitch / (2.0 * N * math.sin(theta / 2.0))
    b = r_in - pitch / (2.0 * N)  # initial point offset
    beta = math.pi - theta  # polygon interior angle
    # initial arc point coordinates
    end_x = p + b * math.cos(beta / 2)
    end_y = -b * math.sin(beta / 2)
    p_end = Point(end_x, end_y)
    # initial point angle
    start_angle = -math.atan2(end_y, end_x)
    start_angle = beta / 2
    delta_y = p * math.cos(theta / 2)  # apothem of the polygon

    layer = 'F.Cu' if (dir == 1) else 'B.Cu'
    p_end.rotate(start_angle)  # turn CW
    p_end.slide_xy(0, -delta_y)  # shift to align top and bottom spirals
    for n in range(1, N * N_turns + 1):
        p_start = p_end.copy()
        p_center = Point(p * math.cos(n * theta), p * math.sin(n * theta))
        p_center.rotate(start_angle)
        p_center.slide_xy(0, -delta_y)
        p_end.rotate_about(p_center, theta)  # end point of the circular arc

        if (dir == -1):
            p_center.y = -p_center.y
            p_start.y = -p_start.y

        # for debug
        # if (n <= N) : # add small pads to mark polygon vertices
        #    sm.add_smd_pad(n, 'rect', p_center, Point(0.1, 0.1))
        #sm.add_smd_pad(N+n, 'rect', p_start, Point(0.1, 0.1))

        sm.add_arc(p_center, p_start, dir * theta, layer, tr_w)
        sf.add_arc(p_center, p_start, dir * theta, layer, tr_w)


def arcs_spiral(N_turns, r_in, pitch, dir, N):
    """
    Draw a spiral approximation using circular arcs.

    Keyword arguments:
    N_turns -- number of turns
    r_in -- internal radius
    pitch -- spacing between conductors centers
    dir -- spiral direction/PCB side
    N -- number of generator polygon vertices

    Returns list of circular arcs descriptions (center, start, angle)

    see
    "Scan Converting Spirals", F. Taponecco and M. Alexa,
    http://wscg.zcu.cz/wscg2002/Papers_2002/F11.pdf
    and
    "Piecewise Circular Approximation of Spirals and Polar Polynomials", F. Taponecco and M. Alexa,
    http://147.228.63.9/wscg2003/papers_2003/i31.pdf
    """

    theta = 2 * math.pi / N  # arc length in radians; polygon central angle
    # polygon radius (distance from the center to a vertex)
    p = pitch / (2.0 * N * math.sin(theta / 2.0))
    b = r_in - pitch / (2.0 * N)  # initial point offset
    beta = math.pi - theta  # polygon interior angle
    # initial arc point coordinates
    end_x = p + b * math.cos(beta / 2)
    end_y = -b * math.sin(beta / 2)
    p_end = Point(end_x, end_y)
    # initial point angle
    start_angle = -math.atan2(end_y, end_x)
    start_angle = beta / 2
    delta_y = p * math.cos(theta / 2)  # apothem of the polygon

    layer = 'F.Cu' if (dir == 1) else 'B.Cu'
    p_end.rotate(start_angle)  # turn CW
    p_end.slide_xy(0, -delta_y)  # shift to align top and bottom spirals
    arcs = []
    # FIXME: make sure that N * N_turns is an integer
    for n in range(1, int(round(N * N_turns)) + 1):
        p_start = p_end.copy()
        p_center = Point(p * math.cos(n * theta), p * math.sin(n * theta))
        p_center.rotate(start_angle)
        p_center.slide_xy(0, -delta_y)
        p_end.rotate_about(p_center, theta)  # end point of the circular arc

        if (dir == -1):
            p_center.y = -p_center.y
            p_start.y = -p_start.y

        # for debug
        # if (n <= N) : # add small pads to mark polygon vertices
        #    sm.add_smd_pad(n, 'rect', p_center, Point(0.1, 0.1))
        #sm.add_smd_pad(N+n, 'rect', p_start, Point(0.1, 0.1))

        arcs.append([p_center, p_start, dir * theta])

        #sm.add_arc(p_center, p_start, dir * theta, layer, tr_w)
        #sf.add_arc(p_center, p_start, dir * theta, layer, tr_w)
    return arcs

########################################

if __name__ == '__main__':
    N_turns = 13  # number of turns (per side)
    tr_w = 2  # trace width
    N = 5  # number of polygon sides
    r_in = 5
    pitch = 3
    dir = 1  # CW or CCW
    PCB_h = 1.6e-3  # PCB thickness
    Cu_t = 35e-6  # copper thickness

    sm = kmodule('test.kicad_mod')
    sf = fh_file('test.inp')
    sm.write_header(name='SIND', descr='spiral inductor', tags='SMD')
    sf.write_header()
    draw_arcs_spiral(N_turns, r_in, pitch, tr_w, N, dir)
    draw_circ_spiral(N_turns, r_in, pitch, tr_w, dir)
    sf.add_ports()

    # compute inner and outer diameter for Mohan's formula
    din = 2 * r_in - tr_w + pitch / 2.0
    dout = 2 * r_in + (2 * N_turns - 0.5) * pitch + tr_w
    ind = calc_ind(N_turns, dout / 1e3, din / 1e3)
    print('din =', din)
    print('dout =', dout)
    print('ind =', ind)
    k = calc_mut(N_turns, PCB_h)
    print('mutual ind =', k * ind, k)

    draw_arcs_spiral(N_turns, r_in, pitch, tr_w, N, -dir)
    sf.add_ports()
    sm.write_refs(0, 0, ref='REF**', value='LLL')

    # center pad (debug)
    sm.add_thru_pad('lc', 'circle', Point(0, 0), Point(0.6, 0.6), 0.3)

    sf.add_frequency()
    sm.close()
    sf.close()

    sf.run()

    freqs, mats = sf.readZc()
    print(freqs)
    print(mats)
