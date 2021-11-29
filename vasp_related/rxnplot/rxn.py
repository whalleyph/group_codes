#! /usr/bin/env python

#from pyx import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_mypaths import line
from mol2latex import mol2latex

#text.set(mode='latex')
#text.preamble(r'\usepackage[version=3]{mhchem}')

scale_ev2kj = 96.4869             # eV to kJ/mole

class Molecules:
    def __init__(self, name, energy, ev2kj=False):
        self.name = name
        if ev2kj:
            self.energy = scale_ev2kj*energy
        else:
            self.energy = energy

    def __str__(self):
        l1 = 'Molecules: %s' % (self.name)
        l2 = 'Energy: %.2f' % (self.energy)
        l = '\n'.join([l1, l2])# + '\n'
        return l

class Reaction:
    def __init__(self, reactants, products, activation_energy = 0, ev2kj=False):
        self.reactants = reactants
        self.products = products
        if ev2kj:
            self.activation_energy = scale_ev2kj*activation_energy
        else:
            self.activation_energy = activation_energy
        self.transition_state_energy = self.activation_energy + self.reactants.energy
        self.reaction_energy = self.products.energy - self.reactants.energy

    def __str__(self):
        l1 = self.reactants.name + ' --> ' + self.products.name
        l2 = '-dE: %.2f - %.2f = %.2f' %  (self.reactants.energy,
                self.products.energy, -self.reaction_energy)
        l3 = 'Ea: %.2f' % (self.activation_energy)

        return '\n'.join([l1, l2, l3]) + '\n'

    def continue_to(self, products, activation_energy = 0):
        """
        For a reaction A + B -> C + D -> E + F, takes the reaction
        A + B -> C + D and the products E + F as input and returns the
        reaction C + D -> E + F with all the energies adjusted properly

        Yet to be tested!!
        """

        return Reaction(self.products, products, activation_energy)
        #self.reactants.name = self.products.name
        #self.reactants.energy = self.products.energy
        #self.products.name = products.name
        #self.products.energy = products.energy
        #self.activation_energy = activation_energy
        #self.transition_state_energy = self.reactants.energy + self.activation_energy
        
        #return self

    def plot_mpl(self, ax, start=0.0, state_width=0.5, total_width=5,
            annotations_reactants=False, annotations_products=True):
        linewidth = 0.1
        width = state_width

        rs = start
        ps = start + total_width - width
        re = rs + width
        pe = ps + width

        print rs, pe

        # state for reactants
        path_reactants = line(ax, (rs, self.reactants.energy),
                                  (re, self.reactants.energy), color='b', lw=2)

        # state for products
        path_products = line(ax, (ps, self.products.energy),
                                 (pe, self.products.energy), color='b', lw=2)

        if self.activation_energy:
            tss = 0.5*(rs + ps)
            tse = tss + width
            # state for transition state
            path_activation_energy = line(ax,
                    (tss, self.transition_state_energy),
                    (tse, self.transition_state_energy), color='r', lw=2)

            # line from reactants to transition state
            reactant_transition_state_connector = line(ax,
                    (re, self.reactants.energy),
                    (tss, self.transition_state_energy), ls='dotted')

            # line from products to transition state
            product_transition_state_connector = line(ax,
                    (ps, self.products.energy),
                    (tse, self.transition_state_energy), ls='dotted')
        else:
            # if no transition state, just connect the reactant and product
            # states
            reactant_product_connector = line(ax,
                    (re, self.reactants.energy),
                    (ps, self.products.energy), ls='dotted')

        # Annotations
        # offsets
        ox, oy = 0.0, 0.1
        # scaling
        #scale = 1.2
        if annotations_products:
            label = mol2latex(self.products.name, bold=True)
            print label
            ax.annotate(label,
                    xy=(pe, self.products.energy),
                    xytext=(pe+ox, self.products.energy+oy),
                    fontsize=14)
            print mol2latex(self.products.name)
        if annotations_reactants:
            label = mol2latex(self.reactants.name, bold=True)
            print label
            ax.annotate(label,
                    xy=(rs, self.reactants.energy),
                    xytext=(rs+ox, self.reactants.energy+oy),
                    fontsize=14)
            print mol2latex(self.reactants.name)

        return ax, ps


#    def plot(self, c=canvas, start=0.0, state_width=0.5, total_width=5,
#            annotations_reactants=False, annotations_products=True):
#
#        linewidth = 0.1
#        width = state_width
#
#        # rs = reactant start
#        # re = reactant end
#        # ps = product start
#        # pe = product end
#        # tss = transition state start
#        # tse = transition state end
#
#        rs = start
#        ps = start + total_width - width
#        re = rs + width
#        pe = ps + width
#
#        print rs, pe
#
#        # state for reactants
#        path_reactants = path.line(rs, self.reactants.energy,
#                re, self.reactants.energy)
#        c.stroke(path_reactants, [style.linewidth(linewidth), color.rgb.blue])
#
#        # state for products
#        path_products = path.line(ps, self.products.energy,
#                pe, self.products.energy)
#        c.stroke(path_products, [style.linewidth(linewidth), color.rgb.blue])
#
#        if self.activation_energy:
#            tss = 0.5*(rs + ps)
#            tse = tss + width
#            # state for transition state
#            path_activation_energy = path.line(
#                    tss, self.transition_state_energy,
#                    tse, self.transition_state_energy)
#            c.stroke(path_activation_energy, [style.linewidth(linewidth), color.rgb.red])
#
#            # line from reactants to transition state
#            reactant_transition_state_connector = path.line(
#                    re, self.reactants.energy,
#                    tss, self.transition_state_energy)
#            c.stroke(reactant_transition_state_connector, [style.linestyle.dashed])
#
#            # line from products to transition state
#            product_transition_state_connector = path.line(
#                    ps, self.products.energy,
#                    tse, self.transition_state_energy)
#            c.stroke(product_transition_state_connector, [style.linestyle.dashed])
#        else:
#            # if no transition state, just connect the reactant and product
#            # states
#            reactant_product_connector = path.line(
#                    re, self.reactants.energy,
#                    ps, self.products.energy)
#            c.stroke(reactant_product_connector, [style.linestyle.dashed])
#
#        # Annotations
#        # offsets
#        ox, oy = 0.0, 0.1
#        # scaling
#        scale = 1.2
#        if annotations_products:
#            c.text(ps+ox, self.products.energy+oy, self.products.name,
#                    [trafo.scale(scale, scale, 0, 0)])
#        if annotations_reactants:
#            c.text(rs+ox, self.reactants.energy+oy, self.reactants.name,
#                    [trafo.scale(scale, scale, 0, 0)])
#
#        return c, ps

class ReactionPathway:
    """ UNDER DEVELOPMENT """
    def __init__(self, rxns):
        self.nrxns = len(rxns)

    def plot(self, width=8):
        self.segments = 2*nrxns + 1
        self.segment_width = float(width)/self.segments
        self.energy_level_bar_width = 0.8*self.segment_width
        

if __name__ == '__main__':
    # (surfacespecies)_(hydrogen_in_supercell)_(hydrogen_at_infinity)
    # g = gas
    # s = surface

    s = 2.5
    ch4_g = Molecules('CH4(g)', 0.0*s)
    ch4_s = Molecules('CH4*', -0.027*s)

    ch3_1h_0h = Molecules('CH3* + H*', 0.014*s)
    ch3_0h_1h = Molecules('CH3* + H*', -0.206*s)

    ch2_1h_1h = Molecules('CH2* + 2H*', -0.122*s)
    ch2_0h_2h =   Molecules('CH2* + 2H*', -0.382*s)

    ch_1h_2h = Molecules('CH* + 3H*', -0.844*s)
    ch_0h_3h = Molecules('CH* + 3H*', -0.952*s)

    c_1h_3h = Molecules('C* + 4H*', -0.566*s)
    c_0h_4h = Molecules('C* + 4H*', -0.722*s)

    rxn10 = ch4_g_to_ch4_s = Reaction(ch4_g, ch4_s, 0.0*s)

    rxn20 = ch4_s_to_ch3_1h_0h = rxn10.continue_to(ch3_1h_0h, 0.86*s)
    rxn30 = ch3_1h_0h_to_ch3_0h_1h = rxn20.continue_to(ch3_0h_1h, 0.0*s)

    rxn40 = ch3_0h_1h_to_ch2_1h_1h = rxn30.continue_to(ch2_1h_1h, 0.65*s)
    rxn50 = ch2_1h_1h_to_ch2_0h_2h = rxn40.continue_to(ch2_0h_2h, 0.0*s)

    rxn60 = ch2_0h_2h_to_ch_1h_2h = rxn50.continue_to(ch_1h_2h, 0.16*s)
    rxn70 = ch_1h_2h_to_ch_0h_3h = rxn60.continue_to(ch_0h_3h, 0.0*s)

    rxn80 = ch_0h_3h_to_ch_1h_3h = rxn70.continue_to(c_1h_3h, 1.06*s)
    rxn90 = c_1h_3h_to_c_0h_4h = rxn80.continue_to(c_0h_4h, 0.0*s)

    for rxn in [rxn10, rxn20, rxn30, rxn40, rxn50, rxn60, rxn70, rxn80, rxn90]:
        print rxn.reactants
        print rxn.products
        print rxn


    w = 0.25
    sw1 = 2*w*1.25
    sw2 = 3*w*1.75
    
    rxns = {'rxn10': {'rxn': rxn10, 'sw': sw1},
            'rxn20': {'rxn': rxn20, 'sw': sw2},
            'rxn30': {'rxn': rxn30, 'sw': sw1},
            'rxn40': {'rxn': rxn40, 'sw': sw2},
            'rxn50': {'rxn': rxn50, 'sw': sw1},
            'rxn60': {'rxn': rxn60, 'sw': sw2},
            'rxn70': {'rxn': rxn70, 'sw': sw1},
            'rxn80': {'rxn': rxn80, 'sw': sw2},
            'rxn90': {'rxn': rxn90, 'sw': sw1}
           }

    c = canvas.canvas()

    start = 0.0
    for key in sorted(rxns.keys()):
        rxn = rxns[key]['rxn']
        print start, w, rxns[key]['sw']
        c = rxn.plot(c, start=start, state_width=w, total_width=rxns[key]['sw'])
        start = start + rxns[key]['sw'] - w
    c.writePDFfile('ch4_dehydrogenation')
