#! /usr/bin/env python

import matplotlib
#matplotlib.use('Agg')
from rxn import Molecules, Reaction
#from pyx import *
import matplotlib.pyplot as plt

# matplotlib customizations
import matplotlib as mpl
mpl.rc('font', size=12, weight='bold')
#mpl.rc('axes', linewidth=2)
mpl.rc('lines', linewidth=2)
mpl.rc('figure.subplot', bottom=0.05, top=0.95, left=0.075, right=0.96)

if __name__ == '__main__':
    # (surfacespecies)_(hydrogen_in_supercell)_(hydrogen_at_infinity)
    # g = gas
    # s = surface

    ev2kj = 96.4869             # eV to kJ/mole
    s = 1 #ev2kj*0.05
    #s = ev2kj#*0.05
    #s = 1
    ch3ch3_g = Molecules(r'C2H6(g)', 0.0*s)
    ch3ch3_s = Molecules(r'C2H6', -0.026*s)

    ch3ch2_1h_0h = Molecules(r'CH3CH2 + 1H', 0.228*s)
    ch3ch2_0h_1h = Molecules(r'CH3CH2 + 1H', 0.047*s)

    ch2ch2_1h_1h = Molecules(r'CH2CH2 + 2H', 0.087*s)
    ch2ch2_0h_2h = Molecules(r'CH2CH2 + 2H', -0.315*s)
    ch3ch_1h_1h = Molecules(r'CH3CH + 2H', -0.286*s)
    ch3ch_0h_2h = Molecules(r'CH3CH + 2H', -0.414*s)

    ch2ch_1h_2h = Molecules(r'CH2CH + 3H', -0.372*s)
    ch2ch_0h_3h = Molecules(r'CH2CH + 3H', -0.641*s)
    ch3c_1h_2h = Molecules(r'CH3C + 3H', -1.084*s)
    ch3c_0h_3h = Molecules(r'CH3C + 3H', -1.136*s)

    chch_1h_3h = Molecules(r'CHCH + 4H', -0.811*s)
    chch_0h_4h = Molecules(r'CHCH + 4H', -1.066*s)
    ch2c_1h_3h = Molecules(r'CH2C + 4H', -0.968*s)
    ch2c_0h_4h = Molecules(r'CH2C + 4H', -1.135*s)

    chc_1h_4h = Molecules(r'CHC + 5H', -0.842*s)
    chc_0h_5h = Molecules(r'CHC + 5H', -1.109*s)
    
    cc_1h_5h = Molecules(r'CC + 6H', -0.217*s)
    cc_0h_6h = Molecules(r'CC + 6H', -0.594*s)

    ch3_ch2_0h_1h = Molecules(r'\ce{CH3$*$ + CH2$*$ + 1H$*$}', 0.511*s)
    ch3_ch_0h_2h = Molecules(r'\ce{CH3$*$ + CH$*$ + 2H$*$}', -0.228*s)
    ch3_c_0h_3h = Molecules(r'\ce{CH3$*$ + C$*$ + 3H$*$}', 0.167*s)
    ch2_ch2_0h_2h = Molecules(r'\ce{CH2$*$ + CH2$*$ + 2H$*$}', -0.177*s)
    ch2_ch_0h_3h = Molecules(r'\ce{CH2$*$ + CH$*$ + 3H$*$}', -0.627*s)
    ch2_c_0h_4h = Molecules(r'\ce{CH2$*$ + C$*$ + 4H$*$}', -0.085*s)
    ch_ch_0h_4h = Molecules(r'\ce{CH$*$ + CH$*$ + 4H$*$}', -1.025*s)
    ch_c_0h_5h = Molecules(r'\ce{CH$*$ + C$*$ + 5H$*$}', -0.462*s)
    c_c_0h_6h = Molecules(r'\ce{C$*$ + C$*$ + 6H$*$}', 0.103*s)

    # adsorption
    r10 = ch3ch3_g_to_ch3ch3_s = Reaction(ch3ch3_g, ch3ch3_s, 0.0*s)

    # dehydrogenation reactions
    # CH3CH3* -> CH3CH2* + 1H*
    r20 = ch3ch3_s_to_ch3ch2_1h_0h = r10.continue_to(ch3ch2_1h_0h, 0.86*s)
    r30 = ch3ch2_1h_0h_to_ch3ch2_0h_1h = r20.continue_to(ch3ch2_0h_1h, 0.0*s)

    # CH3CH2* + 1H* -> CH2CH2* + 2H*
    r40 = ch3ch2_0h_1h_to_ch2ch2_1h_1h = r30.continue_to(ch2ch2_1h_1h, 0.32*s)
    r50 = ch2ch2_1h_1h_to_ch2ch2_0h_2h = r40.continue_to(ch2ch2_0h_2h, 0.0*s)

    # CH3CH2* + 1H* -> CH3CH* + 2H*
    r60 = ch3ch2_0h_1h_to_ch3ch_1h_1h = r30.continue_to(ch3ch_1h_1h, 0.29*s)
    r70 = ch3ch_1h_1h_to_ch3ch_0h_2h = r60.continue_to(ch3ch_0h_2h, 0.0*s)

    # CH3CH* + 2H* -> CH3C* + 3H*
    r80 = ch3ch_0h_2h_to_ch3c_1h_2h = r70.continue_to(ch3c_1h_2h, 0.08*s)
    r90 = ch3c_1h_2h_to_ch3c_0h_3h = r80.continue_to(ch3c_0h_3h, 0.0*s)

    # CH3CH* + 2H* -> CH2CH* + 3H*
    r100 = ch3ch_0h_2h_to_ch2ch_1h_2h = r70.continue_to(ch2ch_1h_2h, 0.49*s)
    r110 = ch2ch_1h_2h_to_ch2ch_0h_3h = r100.continue_to(ch2ch_0h_3h, 0.0*s)

    # CH2CH2* + 2H* -> CH2CH* + 3H*
    r120 = ch2ch2_0h_2h_to_ch2ch_1h_2h = r50.continue_to(ch2ch_1h_2h, 0.37*s)
    r130 = ch2ch_1h_2h_to_ch2ch_0h_3h = r120.continue_to(ch2ch_0h_3h, 0.0*s)

    # CH2CH* + 3H* -> CHCH* + 4H*
    r140 = ch2ch_0h_3h_to_chch_1h_3h = r130.continue_to(chch_1h_3h, 0.70*s)
    r150 = chch_1h_3h_to_chch_0h_4h = r140.continue_to(chch_0h_4h, 0.0*s)

    # CH2CH* + 3H* -> CH2C* + 4H*
    r160 = ch2ch_0h_3h_to_ch2c_1h_3h = r130.continue_to(ch2c_1h_3h, 0.19*s)
    r170 = ch2c_1h_3h_to_ch2c_0h_4h = r160.continue_to(ch2c_0h_4h, 0.0*s)

    # CH3C* + 3H* -> CH2C* + 4H*
    r180 = ch3c_0h_3h_to_ch2c_1h_3h = r90.continue_to(ch2c_1h_3h, 0.64*s)
    r190 = ch2c_1h_3h_to_ch2c_0h_4h = r180.continue_to(ch2c_0h_4h, 0.0*s)

    # CH2C* + 4H* -> CHC* + 5H*
    r200 = ch2c_0h_4h_to_chc_1h_4h = r190.continue_to(chc_1h_4h, 1.0*s)
    r210 = chc_1h_4h_to_chc_0h_5h = r200.continue_to(chc_0h_5h, 0.0*s)

    # CHCH* + 4H* -> CHC* + 5H*
    r220 = chch_0h_4h_to_chc_1h_4h = r150.continue_to(chc_1h_4h, 0.81*s)
    r230 = chch_1h_4h_to_chc_0h_5h = r220.continue_to(chc_0h_5h, 0.0*s)

    # CHC* + 5H* -> CC + 6H*
    r240 = chc_0h_5h_to_cc_1h_5h = r210.continue_to(cc_1h_5h, 1.59*s)
    r250 = cc_1h_5h_to_cc_0h_6h = r240.continue_to(cc_0h_6h, 0.0*s)

    for rxn in [r10, r20, r30, r40, r50, r60, r70, r80, r90, r100, r110, r120,
            r130, r140, r150, r160, r170, r180, r190, r200, r210, r220, r230,
            r240, r250]:
        #print rxn.reactants
        #print rxn.products
        print rxn

    # C-C scission reactions
    # CH3CH2* -> CH3* + CH2*
    rc10 = ch3ch2_0h_1h_to_ch3_ch2_0h_1h = r30.continue_to(ch3_ch2_0h_1h,
            1.07*s)
    # CH3CH* -> CH3* + CH*
    rc20 = ch3ch_0h_2h_to_ch3_ch_0h_2h = r70.continue_to(ch3_ch_0h_2h, 1.06*s)
    # CH3C* -> CH3* + C*
    rc30 = ch3c_0h_3h_to_ch3_c_0h_3h = r90.continue_to(ch3_c_0h_3h, 1.65*s)
    # CH2CH2* -> CH2* + CH2*
    rc40 = ch2ch2_0h_2h_to_ch2_ch2_0h_2h = r50.continue_to(ch2_ch2_0h_2h,
            1.00*s)
    # CH2CH* -> CH2* + CH*
    rc50 = ch2ch_0h_3h_to_ch2_ch_0h_3h = r110.continue_to(ch2_ch_0h_3h, 1.09*s)
    # CH2C* -> CH2* + C*
    rc60 = ch2c_0h_4h_to_ch2_c_0h_4h = r170.continue_to(ch2_c_0h_4h, 1.83*s)
    # CHCH* -> CH* + CH*
    rc70 = chch_0h_4h_to_ch_ch_0h_4h = r150.continue_to(ch_ch_0h_4h, 0.99*s)
    # CHC* -> CH* + C*
    rc80 = chc_0h_5h_to_ch_c_0h_5h = r210.continue_to(ch_c_0h_5h, 1.39*s)
    # CC* -> C* + C*
    rc90 = cc_0h_6h_to_c_c_0h_6h = r250.continue_to(c_c_0h_6h, 1.56*s)
    
    def kJ(eV):
        return eV*96.4869
    
    def update_bx(ax):
        y1, y2 = ax.get_ylim()
        bx.set_ylim(kJ(y1), kJ(y2))
        bx.figure.canvas.draw()

    fig = plt.figure(figsize=(10,7.5))
    #fig = plt.figure()
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
    ax = fig.add_subplot(111)
    bx = ax.twinx()
    sw = 1
    w1 = 3*sw*1.75
    w2 = 2*sw*1.25
    ax.callbacks.connect("ylim_changed", update_bx)
    ax.plot([0,40], [0, 0], 'k:')
    # CH3CH3(g) -> CH3CH3*
    ax, e1 = r10.plot_mpl(ax, start=0.0, state_width=sw, total_width=w2,
            annotations_reactants=True)
    # CH3CH3* -> CH3CH2* + 1H*
    ax, e2 = r20.plot_mpl(ax, start=e1, state_width=sw, total_width=w1)
    ax, e3 = r30.plot_mpl(ax, start=e2, state_width=sw, total_width=w2)
    # CH3CH2* + 1H* -> CH2CH2* + 2H*
    ax, e4 = r40.plot_mpl(ax, start=e3, state_width=sw, total_width=w1) 
    ax, e5 = r50.plot_mpl(ax, start=e4, state_width=sw, total_width=w2)
    # CH3CH2* + 1H* -> CH3CH* + 2H*
    ax, e6 = r60.plot_mpl(ax, start=e3, state_width=sw, total_width=w1) 
    ax, e7 = r70.plot_mpl(ax, start=e4, state_width=sw, total_width=w2)
    # CH3CH* + 2H* -> CH3C* + 3H*
    ax, e8 = r80.plot_mpl(ax, start=e7, state_width=sw, total_width=w1)
    ax, e9 = r90.plot_mpl(ax, start=e8, state_width=sw, total_width=w2)
    # CH3CH* + 2H* -> CH2CH* + 3H*
    ax, e10 = r100.plot_mpl(ax, start=e7, state_width=sw, total_width=w1)
    ax, e11 = r110.plot_mpl(ax, start=e10, state_width=sw, total_width=w2)
    # CH2CH2* + 2H* -> CH2CH* + 3H*
    ax, e12 = r120.plot_mpl(ax, start=e5, state_width=sw, total_width=w1)
    ax, e13 = r130.plot_mpl(ax, start=e12, state_width=sw, total_width=w2)
    # CH2CH* + 3H* -> CHCH* + 4H*
    ax, e14 = r140.plot_mpl(ax, start=e13, state_width=sw, total_width=w1)
    ax, e15 = r150.plot_mpl(ax, start=e14, state_width=sw, total_width=w2)
    # CH2CH* + 3H* -> CH2C* + 4H*
    ax, e16 = r160.plot_mpl(ax, start=e13, state_width=sw, total_width=w1)
    ax, e17 = r170.plot_mpl(ax, start=e16, state_width=sw, total_width=w2)
    # CH3C* + 3H* -> CH2C* + 4H*
    ax, e18 = r180.plot_mpl(ax, start=e9, state_width=sw, total_width=w1)
    ax, e19 = r190.plot_mpl(ax, start=e18, state_width=sw, total_width=w2)
    # CH2C* + 4H* -> CHC* + 5H*
    ax, e20 = r200.plot_mpl(ax, start=e17, state_width=sw, total_width=w1)
    ax, e21 = r210.plot_mpl(ax, start=e20, state_width=sw, total_width=w2)
    # CHCH* + 4H* -> CHC* + 5H*
    ax, e22 = r220.plot_mpl(ax, start=e15, state_width=sw, total_width=w1)
    ax, e23 = r230.plot_mpl(ax, start=e22, state_width=sw, total_width=w2)
    # CHC* + 5H* -> CC + 6H*
    ax, e24 = r240.plot_mpl(ax, start=e23, state_width=sw, total_width=w1)
    ax, e25 = r250.plot_mpl(ax, start=e24, state_width=sw, total_width=w2)

    #ax.set_ylabel(r'$\bf{E_{a}\,[kJ/mole]}$', size='x-large', weight='bold')
    #ax.set_ylabel(r'$\bf{E\,[kJ/mole]}$', size='x-large', weight='bold')
    ax.set_ylabel(r'$\bf{E\,[eV]}$', size='x-large', weight='bold')
    bx.set_ylabel(r'$\bf{E\,[kJ/mole]}$', size='x-large', weight='bold')
    ax.set_xticks([])

    #ax.set_ylim(-6, 5)
    ax.autoscale_view()
    ax.set_xlim(0, 38)
    #plt.tight_layout(h_pad=1.25)
    plt.savefig('c2h6_decomposition_mpl.png')
    plt.savefig('c2h6_decomposition_mpl.pdf')
    plt.savefig('c2h6_decomposition_mpl.eps')
    plt.savefig('c2h6_decomposition_mpl.svg')
    plt.show()



def pyx_plotting():
    # === Plotting ===
    #c = canvas.canvas()
    sw = 0.45                   # state width
    w1 = 3*sw*1.75              # with activation energy
    w2 = 2*sw*1.25              # no activation energy

    # canvas
    emin = -1.2
    emax = 1.2
    painter = graph.axis.painter.regular(
            #outerticklength=graph.axis.painter.ticklength.normal,
            basepathattrs=[style.linewidth.THick, deco.earrow.large,
                deco.barrow.large])
    c = graph.axis.pathaxis(path.line(-0.25, s*emax, -0.25, s*emin),
            graph.axis.linear(min=emax*ev2kj, max=emin*ev2kj,
            title=r'E [kJ/mole]', painter=painter))

    # zero reference
    c.stroke(path.line(-0.25, 0, 19, 0), [style.linestyle.dashed])

    # CH3CH3(g) -> CH3CH3*
    c, e1 = r10.plot(c, start=0.0, state_width=sw, total_width=w2,
            annotations_reactants=True)
    # CH3CH3* -> CH3CH2* + 1H*
    c, e2 = r20.plot(c, start=e1, state_width=sw, total_width=w1)
    c, e3 = r30.plot(c, start=e2, state_width=sw, total_width=w2)
    # CH3CH2* + 1H* -> CH2CH2* + 2H*
    c, e4 = r40.plot(c, start=e3, state_width=sw, total_width=w1) 
    c, e5 = r50.plot(c, start=e4, state_width=sw, total_width=w2)
    # CH3CH2* + 1H* -> CH3CH* + 2H*
    c, e6 = r60.plot(c, start=e3, state_width=sw, total_width=w1) 
    c, e7 = r70.plot(c, start=e4, state_width=sw, total_width=w2)
    # CH3CH* + 2H* -> CH3C* + 3H*
    c, e8 = r80.plot(c, start=e7, state_width=sw, total_width=w1)
    c, e9 = r90.plot(c, start=e8, state_width=sw, total_width=w2)
    # CH3CH* + 2H* -> CH2CH* + 3H*
    c, e10 = r100.plot(c, start=e7, state_width=sw, total_width=w1)
    c, e11 = r110.plot(c, start=e10, state_width=sw, total_width=w2)
    # CH2CH2* + 2H* -> CH2CH* + 3H*
    c, e12 = r120.plot(c, start=e5, state_width=sw, total_width=w1)
    c, e13 = r130.plot(c, start=e12, state_width=sw, total_width=w2)
    # CH2CH* + 3H* -> CHCH* + 4H*
    c, e14 = r140.plot(c, start=e13, state_width=sw, total_width=w1)
    c, e15 = r150.plot(c, start=e14, state_width=sw, total_width=w2)
    # CH2CH* + 3H* -> CH2C* + 4H*
    c, e16 = r160.plot(c, start=e13, state_width=sw, total_width=w1)
    c, e17 = r170.plot(c, start=e16, state_width=sw, total_width=w2)
    # CH3C* + 3H* -> CH2C* + 4H*
    c, e18 = r180.plot(c, start=e9, state_width=sw, total_width=w1)
    c, e19 = r190.plot(c, start=e18, state_width=sw, total_width=w2)
    # CH2C* + 4H* -> CHC* + 5H*
    c, e20 = r200.plot(c, start=e17, state_width=sw, total_width=w1)
    c, e21 = r210.plot(c, start=e20, state_width=sw, total_width=w2)
    # CHCH* + 4H* -> CHC* + 5H*
    c, e22 = r220.plot(c, start=e15, state_width=sw, total_width=w1)
    c, e23 = r230.plot(c, start=e22, state_width=sw, total_width=w2)
    # CHC* + 5H* -> CC + 6H*
    c, e24 = r240.plot(c, start=e23, state_width=sw, total_width=w1)
    c, e25 = r250.plot(c, start=e24, state_width=sw, total_width=w2)

    # C-C scission reactions
    # CH3CH2* -> CH3* + CH2*
    c, e26 = rc10.plot(c, start=e3, state_width=sw, total_width=w1)
    # CH3CH* -> CH3* + CH*
    c, e27 = rc20.plot(c, start=e7, state_width=sw, total_width=w1)
    # CH3C* -> CH3* + C*
    c, e28 = rc30.plot(c, start=e9, state_width=sw, total_width=w1)
    # CH2CH2* -> CH2* + CH2*
    c, e29 = rc40.plot(c, start=e5, state_width=sw, total_width=w1)
    # CH2CH* -> CH2* + CH*
    c, e30 = rc50.plot(c, start=e13, state_width=sw, total_width=w1)
    # CH2C* -> CH2* + C*
    c, e31 = rc60.plot(c, start=e17, state_width=sw, total_width=w1)
    # CHCH* -> CH* + CH*
    c, e32 = rc70.plot(c, start=e15, state_width=sw, total_width=w1)
    # CHC* -> CH* + C*
    c, e33 = rc80.plot(c, start=e23, state_width=sw, total_width=w1)
    # CC* -> C* + C*
    c, e34 = rc90.plot(c, start=e25, state_width=sw, total_width=w1)

    c.writePDFfile('c2h6_decomposition')

    #start = 0.0
    #for key in sorted(rxns.keys()):
    #    rxn = rxns[key]['rxn']
    #    print start, w, rxns[key]['sw']
    #    c = rxn.plot(c, start=start, state_width=w, total_width=rxns[key]['sw'])
    #    start = start + rxns[key]['sw'] - w
    #c.writePDFfile('ch4_dehydrogenation')
