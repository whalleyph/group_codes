#! /bin/env python

import sys
import os.path, copy

def get_ref_eng(ref):
    refdat = '~/proj/database/refeng_eles_' + ref + '.dat'
    print refdat

    lines = open (os.path.expanduser(refdat), 'r').readlines()
    lines = [line.strip().split() for line in lines if line.strip()[0] != '#']
    refs = {}

    for line in lines:
      refs[line[0]] = float(line[1])

    return refs

def hf(out, ref):

  refs = get_ref_eng(ref)

  composition = out.composition()
  sumnions = sum(out.nions())
  e0 = out.get_final_e0()


  href = 0.0
  sumnions = 0
  for ion, nions in composition:
    if refs.has_key(ion):
      sumnions += nions
      href = href + nions*refs[ion]
    else:
      print "%s %s %s" % ('No reference available for', ion, 'returning 0.0')
      href = copy.copy(e0)
      break
  return e0, href, e0-href, (e0-href)/sumnions

if __name__ == '__main__':

  import sys
  sys.path.append('/home/cande/python/mypackages/vasp/tags/0.2')
  import outcar

  if len(sys.argv) != 3:
    print "Usage: %s <OUTCAR>" % (sys.argv[0].split(os.path.sep)[-1])
    sys.exit(1)

  o = sys.argv[1]
  ref = sys.argv[2]
  out = outcar.outcar().read(o)
  composition = out.composition()
  e0 = out.get_final_e0()
  e0, h0_href, hf_per_supercell, hf_per_atom = hf(out, ref)

  print "Composition of supercell %s   " % (composition)
  print "H0 per supercell = %.3f eV" % (e0)
  print "Hf per atom      = %.0f meV" % (hf_per_atom*1e3)
  print "Hf per supercell = %.3f eV" % (hf_per_supercell)
  print "Hf reference     = %.3f eV" % (h0_href)
