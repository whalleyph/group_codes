#! /usr/bin/latex

import re

"""script to convert molecular to formulat to latex labels suitable for use in
matplotlib"""

def mol2latex(s, bold=False, roman=False):

    # Match all numbers and replace them as subscripts
    nch = re.compile(r'([C,H])([0-9])')
    l = nch.sub(r'\1_{\2}', s.strip())

    if not bold and not roman:
        l = r'$' + l + '$'
        return l

    if bold and roman:
        l = r'$\mathbf{\mathrm{' + l + '}}$'
        return l

    if bold:
        # l = r'$\mathbf{' + l + '}$'
        l = r'$\bf{' + l + '}$'
        return l

    if roman:
        l = r'$\mathrm{' + l + '}$'
        return l

if __name__ == '__main__':
    print mol2latex('C2H6* + 2H', bold=True, roman=True)
