#!/usr/bin/env python

# Licensed under GNU GPL version 2 or later

import sys
import math

useNLopt = True
try:  # check if nlopt is available
    import nlopt
except ImportError:
    useNLopt = False

if (useNLopt):
    from numpy import *  # needed by nlopt

from PyQt4 import QtGui
from PyQt4.QtCore import QThread

import dos
import design


class kSpiralCalc(QtGui.QMainWindow, design.Ui_MainWindow):

    def __init__(self):
        super(self.__class__, self).__init__()
        self.setupUi(self)

        # File
        self.actionExit.triggered.connect(QtGui.qApp.quit)
        self.actionSave_module.triggered.connect(self.writeModule)

        # some default values
        self.nTurnsLineEdit.setText('13')
        self.innerRadiusLineEdit.setText('5')
        self.pitchLineEdit.setText('3')
        self.spacingLineEdit.setText('1')
        self.traceWidthLineEdit.setText('2')
        self.cuThicknessLineEdit.setText('35')
        self.pcbThicknessLineEdit.setText('1.6')
        self.nLayersLineEdit.setText('1')

        self.freqLineEdit.textChanged.connect(self.updateSkinDepth)
        self.freqLineEdit.setText('1.0')

        self.minSpacingLineEdit.setText('0.15')

        self.drawTolLineEdit.setText('0.1')

        self.nTurnsLineEdit.textChanged.connect(self.estimateInductance)
        self.innerRadiusLineEdit.textChanged.connect(self.estimateInductance)
        self.pitchLineEdit.textChanged.connect(self.estimateInductance)
        self.spacingLineEdit.textChanged.connect(self.estimateInductance)
        self.traceWidthLineEdit.textChanged.connect(self.estimateInductance)
        self.nLayersLineEdit.textChanged.connect(self.estimateInductance)

        self.estimateInductance()

        self.pitchLineEdit.textChanged.connect(self.updateSpacing)
        # if trace width is changed keep pitch constant and update the
        # resulting spacing
        self.traceWidthLineEdit.textChanged.connect(self.updateSpacing)
        self.spacingLineEdit.textChanged.connect(self.updatePitch)

        self.runSimBtn.clicked.connect(self.runSimulation)
        self.optimizeBtn.clicked.connect(self.runOptimization)
        if (not useNLopt):
            self.optimizeBtn.setEnabled(False)  # disable optimization button

        self.statusBar().showMessage("Ready.")

    def updateSpacing(self):
        pitch = float(self.pitchLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        spacing = pitch - traceWidth
        self.spacingLineEdit.blockSignals(True)
        self.spacingLineEdit.setText(str(spacing))
        self.spacingLineEdit.blockSignals(False)

    def updatePitch(self):
        spacing = float(self.spacingLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        pitch = spacing + traceWidth
        self.pitchLineEdit.blockSignals(True)
        self.pitchLineEdit.setText(str(pitch))
        self.pitchLineEdit.blockSignals(False)

    def updateSkinDepth(self):
        freq = float(self.freqLineEdit.text()) * 1e6
        sigma = 5.8e7
        mu0 = 4.0e-7 * math.pi
        self.delta = 1.0 / math.sqrt(math.pi * freq * mu0 * sigma)  # in meters
        self.skinDepthLineEdit.setText(str(self.delta * 1e3))

    def runSimulation(self):
        self.statusBar().showMessage("Simulating...")
        # update GUI to show changes status bar message
        QtGui.QApplication.processEvents()
        self.simulate()
        self.statusBar().showMessage("Ready.")

    def simulate(self):
        nTurns = float(self.nTurnsLineEdit.text())
        innerRadius = float(self.innerRadiusLineEdit.text())
        pitch = float(self.pitchLineEdit.text())
        spacing = float(self.spacingLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        cuThickness = float(self.cuThicknessLineEdit.text())
        pcbThickness = float(self.pcbThicknessLineEdit.text())
        nLayers = int(self.nLayersLineEdit.text())
        freq = float(self.freqLineEdit.text())
        d = float(self.drawTolLineEdit.text())

        # compute smaller filament size (see fasthenry docs)
        nw = 2  # fasthenry default
        nh = 2  # fasthenry default

        nwinc = 1  # needs to be odd
        while True:
            n = (nwinc - 1) / 2  # rounding down
            nfils = (2.0 - nw**n * (1.0 + nw)) / (1 - nw)
            fw = 1e-3 * traceWidth / nfils  # smaller width filament size
            if (fw < self.delta):
                break
            nwinc = nwinc + 2

        nhinc = 1  # needs to be odd
        while True:
            n = (nhinc - 1) / 2  # rounding down
            nfils = (2.0 - nh**n * (1.0 + nh)) / (1 - nh)
            fh = 1e-6 * cuThickness / nfils  # smaller height filament
            if (fh < self.delta):
                break
            nhinc = nhinc + 2

        # print 'nwinc =', nwinc
        # print 'nhinc =', nhinc

        sf = dos.fh_file('test.inp')
        sf.write_header(nwinc, nhinc)
        #draw_arcs_spiral(N_turns, r_in, pitch, tr_w, N, dir)
        dir = 1

        vx = dos.circ_spiral(nTurns, innerRadius, pitch, dir, d)
        sf.add_circ_spiral(vx, 1, traceWidth, cuThickness * 1e-3, pcbThickness)
        sf.add_ports()
        if (nLayers == 2):
            vx = dos.circ_spiral(nTurns, innerRadius, pitch, -dir, d)
            sf.add_circ_spiral(
                vx,
                2,
                traceWidth,
                cuThickness * 1e-3,
                pcbThickness)
            sf.add_ports()

        sf.add_frequency(freq * 1e6)
        sf.close()
        sf.run()

        freqs, mats = sf.readZc()

        if (nLayers == 1):
            Xs = mats[0][0][0].imag
            Rs = mats[0][0][0].real
        else:
            Zs1 = mats[0][0][0]
            Zs2 = mats[0][1][1]
            M1 = mats[0][0][1]
            M2 = mats[0][1][0]
            # inductors in series, aiding; minus sign due to coils ports order
            Zs = Zs1 + Zs2 - (M1 + M2)
            Xs = Zs.imag
            Rs = Zs.real

        Ls = 1e6 * Xs / (2.0 * math.pi * freqs[0])  # uH
        Q = Xs / Rs

        self.simIndLineEdit.setText("%.3e" % Ls)
        self.simResLineEdit.setText("%.3e" % Rs)
        self.simQLineEdit.setText("%.1e" % Q)

        return Ls

    def runOptimization(self):
        self.statusBar().showMessage("Optimizing...")
        # update GUI to show changes status bar message
        QtGui.QApplication.processEvents()

        nTurns = float(self.nTurnsLineEdit.text())
        innerRadius = float(self.innerRadiusLineEdit.text())
        pitch = float(self.pitchLineEdit.text())
        spacing = float(self.spacingLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        cuThickness = float(self.cuThicknessLineEdit.text())
        pcbThickness = float(self.pcbThicknessLineEdit.text())

        targetInd = float(self.desiredIndLineEdit.text())

        def errfunc(x, grad):
            if grad.size > 0:
                grad = Null
            self.spacingLineEdit.setText(str(x[0]))
            QtGui.QApplication.processEvents()  # update GUI
            ind = self.simulate()
            err = math.fabs(ind - targetInd)
            return err

        opt = nlopt.opt(nlopt.LN_COBYLA, 1)
        minSpacing = float(self.minSpacingLineEdit.text())
        opt.set_lower_bounds([minSpacing])
        opt.set_min_objective(errfunc)
        opt.set_xtol_rel(1e-2)
        x = opt.optimize([spacing])
        minf = opt.last_optimum_value()
        print "optimum at ", x[0]
        print "minimum value = ", minf
        print "result code = ", opt.last_optimize_result()

        self.spacingLineEdit.setText(str(x[0]))

        self.statusBar().showMessage("Ready.")

    def writeModule(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Module', '.', 'Footprint (*.kicad_mod);;Any File (*)')
        if (not fname):
            return
        if (not fname.endsWith('.kicad_mod')):
            fname = fname + '.kicad_mod'
        
        nTurns = float(self.nTurnsLineEdit.text())
        innerRadius = float(self.innerRadiusLineEdit.text())
        pitch = float(self.pitchLineEdit.text())
        spacing = float(self.spacingLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        nLayers = int(self.nLayersLineEdit.text())
        d = float(self.drawTolLineEdit.text())

        sm = dos.kmodule(fname)
        sm.write_header(name='SIND', descr='spiral inductor', tags='SMD')
        dir = 1
        vx = dos.circ_spiral(nTurns, innerRadius, pitch, dir, d)
        pad1 = vx[0]
        sm.add_circ_spiral(vx, 'F.Cu', traceWidth)
        if (nLayers == 2):
            vx = dos.circ_spiral(nTurns, innerRadius, pitch, -dir, d)
            sm.add_circ_spiral(vx, 'B.Cu', traceWidth)
            sm.add_thru_pad('lc', 'circle', vx[0], dos.Point(0.6, 0.6), 0.3)

        pad2 = vx[-1]
        # add SMD pads at the beginning and end
        padSize = dos.Point(traceWidth/2.0, traceWidth/2.0)
        sm.add_smd_pad('1', 'rect', pad1, padSize)
        sm.add_smd_pad('2', 'rect', pad2, padSize)
        
        #draw_arcs_spiral(nTurns, innerRadius, pitch, traceWidth, N, dir)
        sm.write_refs(0, 0, ref='REF**', value='LLL')
        sm.close()

    def estimateInductance(self):
        nTurns = float(self.nTurnsLineEdit.text())
        innerRadius = float(self.innerRadiusLineEdit.text())
        pitch = float(self.pitchLineEdit.text())
        traceWidth = float(self.traceWidthLineEdit.text())
        cuThickness = float(self.cuThicknessLineEdit.text())
        pcbThickness = float(self.pcbThicknessLineEdit.text())
        nLayers = int(self.nLayersLineEdit.text())

        # compute inner and outer diameter for Mohan's formula
        din = 2 * innerRadius - traceWidth + pitch / 2.0
        dout = 2 * innerRadius + (2 * nTurns - 0.5) * pitch + traceWidth
        ind = dos.calc_ind(nTurns, dout / 1e3, din / 1e3)
        # print 'din =', din
        # print 'dout =', dout
        # print 'ind =', ind

        if (nLayers == 1):
            indtot = ind  # single-lyer inductor
        else:  # two-layer inductor
            k = dos.calc_mut(nTurns, pcbThickness * 1e-3)
            indtot = 2.0 * ind * (1.0 + k)

        self.estIndLineEdit.setText("%.3e" % (indtot * 1e6))


def main():
    app = QtGui.QApplication(sys.argv)
    form = kSpiralCalc()
    form.show()
    app.exec_()

if __name__ == '__main__':
    main()
