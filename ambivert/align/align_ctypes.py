#!/usr/bin/env python3
# encoding: utf-8
"""
align_ctypes.py

Created by Toby Sargeant.
Copyright (c) 2013-2015  Toby Sargeant and The University of Melbourne. All rights reserved.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

"""
from ctypes import *
import os
import sys

from distutils.ccompiler import new_compiler

__author__ = "Toby Sargeant"
__copyright__ = "Copyright 2013-2015, Toby Sargeant and The University of Melbourne"
__credits__ = ["Toby Sargeant","Matthew Wakefield",]
__license__ = "GPLv3"
__version__ = "0.5b1"
__maintainer__ = "Matthew Wakefield"
__email__ = "matthew.wakefield@unimelb.edu.au"
__status__ = "Development"

if hasattr(sys, 'pypy_version_info'):
    align_c = cdll.LoadLibrary(
      new_compiler().library_filename('align_c.' + 'pypy-{major}{minor}'.format(major=sys.pypy_version_info[0], minor=sys.pypy_version_info[1]),
                                      lib_type = 'shared',
                                      output_dir = os.path.split(__file__)[0]))
else:
    align_c = cdll.LoadLibrary(
      new_compiler().library_filename('align_c',
                                      lib_type = 'shared' if sys.platform == 'darwin' else 'dylib', #different subclasses need different flags
                                      output_dir = os.path.split(__file__)[0]))

A_GAP = 0
B_GAP = 1
MATCH = 2

class AlignFrag(Structure): pass

AlignFrag._fields_ = [
  ("next", POINTER(AlignFrag)),
  ("type", c_int),
  ("sa_start", c_int),
  ("sb_start", c_int),
  ("hsp_len", c_int)
]

class Alignment(Structure): pass

Alignment._fields_ = [
  ("align_frag", POINTER(AlignFrag)),
  ("frag_count", c_int),
  ("score", c_int)
]

align_raw = align_c.align_raw
align_raw.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, POINTER(c_int), c_int, c_int ]
align_raw.restype = POINTER(Alignment)

align = align_c.align
align.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, c_ubyte * 256, POINTER(c_int), c_int, c_int ]
align.restype = POINTER(Alignment)



local_align = align_c.local_align
local_align.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, c_ubyte * 256, POINTER(c_int), c_int, c_int ]
local_align.restype = POINTER(Alignment)

local_align_raw = align_c.local_align_raw
local_align_raw.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, POINTER(c_int), c_int, c_int ]
local_align_raw.restype = POINTER(Alignment)



global_align = align_c.global_align
global_align.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, c_ubyte * 256, POINTER(c_int), c_int, c_int ]
global_align.restype = POINTER(Alignment)

global_align_raw = align_c.global_align_raw
global_align_raw.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, POINTER(c_int), c_int, c_int ]
global_align_raw.restype = POINTER(Alignment)



glocal_align = align_c.glocal_align
glocal_align.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, c_ubyte * 256, POINTER(c_int), c_int, c_int ]
glocal_align.restype = POINTER(Alignment)

glocal_align_raw = align_c.glocal_align_raw
glocal_align_raw.argtypes = [ c_char_p, c_int, c_char_p, c_int, c_int, POINTER(c_int), c_int, c_int ]
glocal_align_raw.restype = POINTER(Alignment)



alignment_free = align_c.alignment_free
alignment_free.argtypes = [ POINTER(Alignment) ]

def format_alignment(alignment, s1, s2, sequence_map, scoring_matrix, match_character = None):
    s1 = str(s1)
    s2 = str(s2)
    frag = alignment[0].align_frag
    aln = ['','','']
    while frag:
      frag = frag[0]
      if frag.type == MATCH:
        f1 = s1[frag.sa_start:frag.sa_start + frag.hsp_len]
        f2 = s2[frag.sb_start:frag.sb_start + frag.hsp_len]
        
        aln[0] += f1
        aln[2] += f2
        for a,b in zip(f1, f2):
          if a == b and not match_character:
            aln[1] += a
          elif a == b and match_character:
            aln[1] += match_character
          elif scoring_matrix[sequence_map[1][ord(a)] * sequence_map[0] + sequence_map[1][ord(b)]] > 0:
            aln[1] += '+'
          else:
            aln[1] += ' '
      
      elif frag.type == A_GAP:
        aln[0] += '-' * frag.hsp_len
        aln[1] += ' ' * frag.hsp_len
        aln[2] += s2[frag.sb_start:frag.sb_start + frag.hsp_len]
      elif frag.type == B_GAP:
        aln[0] += s1[frag.sa_start:frag.sa_start + frag.hsp_len]
        aln[1] += ' ' * frag.hsp_len
        aln[2] += '-' * frag.hsp_len
      
      frag = frag.next
    
    return aln[0], aln[1], aln[2]


def make_map(s, unknown, case_insensitive):
  if unknown is None:
    unknown = len(s)
  elif isinstance(unknown, str):
    unknown = s.find(unknown)
    if unknown == -1:
      unknown = len(s)
  else:
    unknown = min(max(0, unknown), 255)

  MAP = (c_ubyte * 256)()

  for i in range(256): MAP[i] = unknown

  for i,c in enumerate(s):
    MAP[ord(c)] = i
    if case_insensitive: MAP[ord(c.lower())] = i

  return max(len(s), unknown), MAP

def make_DNA_scoring_matrix(match,mismatch,nmatch):
 return (c_int * (5 * 5))(*(
    match,      mismatch,   mismatch,   mismatch,   nmatch,
    mismatch,   match,      mismatch,   mismatch,   nmatch,
    mismatch,   mismatch,   match,      mismatch,   nmatch,
    mismatch,   mismatch,   mismatch,   match,      nmatch,
    nmatch,     nmatch,     nmatch,     nmatch,     nmatch
    ))
 
DNA_MAP = make_map('ACGTN', 'N', True)

DNA_SCORE = make_DNA_scoring_matrix(match=1,mismatch=-4,nmatch=0)

PROT_MAP = make_map('ARNDCQEGHILKMFPSTWYVBZX*', 'X', True)

BLOSUM45 = (c_int * (24 * 24))(*(
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X   *
   5,  -2,  -1,  -2,  -1,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -2,  -2,   0,  -1,  -1,  -1,  -5,
  -2,   7,   0,  -1,  -3,   1,   0,  -2,   0,  -3,  -2,   3,  -1,  -2,  -2,  -1,  -1,  -2,  -1,  -2,  -1,   1,  -1,  -5,
  -1,   0,   6,   2,  -2,   0,   0,   0,   1,  -2,  -3,   0,  -2,  -2,  -2,   1,   0,  -4,  -2,  -3,   5,   0,  -1,  -5,
  -2,  -1,   2,   7,  -3,   0,   2,  -1,   0,  -4,  -3,   0,  -3,  -4,  -1,   0,  -1,  -4,  -2,  -3,   6,   1,  -1,  -5,
  -1,  -3,  -2,  -3,  12,  -3,  -3,  -3,  -3,  -3,  -2,  -3,  -2,  -2,  -4,  -1,  -1,  -5,  -3,  -1,  -2,  -3,  -1,  -5,
  -1,   1,   0,   0,  -3,   6,   2,  -2,   1,  -2,  -2,   1,   0,  -4,  -1,   0,  -1,  -2,  -1,  -3,   0,   4,  -1,  -5,
  -1,   0,   0,   2,  -3,   2,   6,  -2,   0,  -3,  -2,   1,  -2,  -3,   0,   0,  -1,  -3,  -2,  -3,   1,   5,  -1,  -5,
   0,  -2,   0,  -1,  -3,  -2,  -2,   7,  -2,  -4,  -3,  -2,  -2,  -3,  -2,   0,  -2,  -2,  -3,  -3,  -1,  -2,  -1,  -5,
  -2,   0,   1,   0,  -3,   1,   0,  -2,  10,  -3,  -2,  -1,   0,  -2,  -2,  -1,  -2,  -3,   2,  -3,   0,   0,  -1,  -5,
  -1,  -3,  -2,  -4,  -3,  -2,  -3,  -4,  -3,   5,   2,  -3,   2,   0,  -2,  -2,  -1,  -2,   0,   3,  -3,  -3,  -1,  -5,
  -1,  -2,  -3,  -3,  -2,  -2,  -2,  -3,  -2,   2,   5,  -3,   2,   1,  -3,  -3,  -1,  -2,   0,   1,  -3,  -2,  -1,  -5,
  -1,   3,   0,   0,  -3,   1,   1,  -2,  -1,  -3,  -3,   5,  -1,  -3,  -1,  -1,  -1,  -2,  -1,  -2,   0,   1,  -1,  -5,
  -1,  -1,  -2,  -3,  -2,   0,  -2,  -2,   0,   2,   2,  -1,   6,   0,  -2,  -2,  -1,  -2,   0,   1,  -2,  -1,  -1,  -5,
  -2,  -2,  -2,  -4,  -2,  -4,  -3,  -3,  -2,   0,   1,  -3,   0,   8,  -3,  -2,  -1,   1,   3,   0,  -3,  -3,  -1,  -5,
  -1,  -2,  -2,  -1,  -4,  -1,   0,  -2,  -2,  -2,  -3,  -1,  -2,  -3,   9,  -1,  -1,  -3,  -3,  -3,  -2,  -1,  -1,  -5,
   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -3,  -1,  -2,  -2,  -1,   4,   2,  -4,  -2,  -1,   0,   0,  -1,  -5,
   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,   2,   5,  -3,  -1,   0,   0,  -1,  -1,  -5,
  -2,  -2,  -4,  -4,  -5,  -2,  -3,  -2,  -3,  -2,  -2,  -2,  -2,   1,  -3,  -4,  -3,  15,   3,  -3,  -4,  -2,  -1,  -5,
  -2,  -1,  -2,  -2,  -3,  -1,  -2,  -3,   2,   0,   0,  -1,   0,   3,  -3,  -2,  -1,   3,   8,  -1,  -2,  -2,  -1,  -5,
   0,  -2,  -3,  -3,  -1,  -3,  -3,  -3,  -3,   3,   1,  -2,   1,   0,  -3,  -1,   0,  -3,  -1,   5,  -3,  -3,  -1,  -5,
  -1,  -1,   5,   6,  -2,   0,   1,  -1,   0,  -3,  -3,   0,  -2,  -3,  -2,   0,   0,  -4,  -2,  -3,   5,   1,  -1,  -5,
  -1,   1,   0,   1,  -3,   4,   5,  -2,   0,  -3,  -2,   1,  -1,  -3,  -1,   0,  -1,  -2,  -2,  -3,   1,   5,  -1,  -5,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -5,
  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,   1
))

BLOSUM50 = (c_int * (24 * 24))(*(
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *
   5,  -2,  -1,  -2,  -1,  -1,  -1,   0,  -2,  -1,  -2,  -1,  -1,  -3,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,  -1,  -5,
  -2,   7,  -1,  -2,  -4,   1,   0,  -3,   0,  -4,  -3,   3,  -2,  -3,  -3,  -1,  -1,  -3,  -1,  -3,  -1,   0,  -1,  -5,
  -1,  -1,   7,   2,  -2,   0,   0,   0,   1,  -3,  -4,   0,  -2,  -4,  -2,   1,   0,  -4,  -2,  -3,   5,   0,  -1,  -5,
  -2,  -2,   2,   8,  -4,   0,   2,  -1,  -1,  -4,  -4,  -1,  -4,  -5,  -1,   0,  -1,  -5,  -3,  -4,   6,   1,  -1,  -5,
  -1,  -4,  -2,  -4,  13,  -3,  -3,  -3,  -3,  -2,  -2,  -3,  -2,  -2,  -4,  -1,  -1,  -5,  -3,  -1,  -3,  -3,  -1,  -5,
  -1,   1,   0,   0,  -3,   7,   2,  -2,   1,  -3,  -2,   2,   0,  -4,  -1,   0,  -1,  -1,  -1,  -3,   0,   4,  -1,  -5,
  -1,   0,   0,   2,  -3,   2,   6,  -3,   0,  -4,  -3,   1,  -2,  -3,  -1,  -1,  -1,  -3,  -2,  -3,   1,   5,  -1,  -5,
   0,  -3,   0,  -1,  -3,  -2,  -3,   8,  -2,  -4,  -4,  -2,  -3,  -4,  -2,   0,  -2,  -3,  -3,  -4,  -1,  -2,  -1,  -5,
  -2,   0,   1,  -1,  -3,   1,   0,  -2,  10,  -4,  -3,   0,  -1,  -1,  -2,  -1,  -2,  -3,   2,  -4,   0,   0,  -1,  -5,
  -1,  -4,  -3,  -4,  -2,  -3,  -4,  -4,  -4,   5,   2,  -3,   2,   0,  -3,  -3,  -1,  -3,  -1,   4,  -4,  -3,  -1,  -5,
  -2,  -3,  -4,  -4,  -2,  -2,  -3,  -4,  -3,   2,   5,  -3,   3,   1,  -4,  -3,  -1,  -2,  -1,   1,  -4,  -3,  -1,  -5,
  -1,   3,   0,  -1,  -3,   2,   1,  -2,   0,  -3,  -3,   6,  -2,  -4,  -1,   0,  -1,  -3,  -2,  -3,   0,   1,  -1,  -5,
  -1,  -2,  -2,  -4,  -2,   0,  -2,  -3,  -1,   2,   3,  -2,   7,   0,  -3,  -2,  -1,  -1,   0,   1,  -3,  -1,  -1,  -5,
  -3,  -3,  -4,  -5,  -2,  -4,  -3,  -4,  -1,   0,   1,  -4,   0,   8,  -4,  -3,  -2,   1,   4,  -1,  -4,  -4,  -1,  -5,
  -1,  -3,  -2,  -1,  -4,  -1,  -1,  -2,  -2,  -3,  -4,  -1,  -3,  -4,  10,  -1,  -1,  -4,  -3,  -3,  -2,  -1,  -1,  -5,
   1,  -1,   1,   0,  -1,   0,  -1,   0,  -1,  -3,  -3,   0,  -2,  -3,  -1,   5,   2,  -4,  -2,  -2,   0,   0,  -1,  -5,
   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   2,   5,  -3,  -2,   0,   0,  -1,  -1,  -5,
  -3,  -3,  -4,  -5,  -5,  -1,  -3,  -3,  -3,  -3,  -2,  -3,  -1,   1,  -4,  -4,  -3,  15,   2,  -3,  -5,  -2,  -1,  -5,
  -2,  -1,  -2,  -3,  -3,  -1,  -2,  -3,   2,  -1,  -1,  -2,   0,   4,  -3,  -2,  -2,   2,   8,  -1,  -3,  -2,  -1,  -5,
   0,  -3,  -3,  -4,  -1,  -3,  -3,  -4,  -4,   4,   1,  -3,   1,  -1,  -3,  -2,   0,  -3,  -1,   5,  -3,  -3,  -1,  -5,
  -2,  -1,   5,   6,  -3,   0,   1,  -1,   0,  -4,  -4,   0,  -3,  -4,  -2,   0,   0,  -5,  -3,  -3,   6,   1,  -1,  -5,
  -1,   0,   0,   1,  -3,   4,   5,  -2,   0,  -3,  -3,   1,  -1,  -4,  -1,   0,  -1,  -2,  -2,  -3,   1,   5,  -1,  -5,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -5,
  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,   1
))

BLOSUM62 = (c_int * (24 * 24))(*(
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *
   4,  -1,  -2,  -2,   0,  -1,  -1,   0,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,   0,  -4,
  -1,   5,   0,  -2,  -3,   1,   0,  -2,   0,  -3,  -2,   2,  -1,  -3,  -2,  -1,  -1,  -3,  -2,  -3,  -1,   0,  -1,  -4,
  -2,   0,   6,   1,  -3,   0,   0,   0,   1,  -3,  -3,   0,  -2,  -3,  -2,   1,   0,  -4,  -2,  -3,   3,   0,  -1,  -4,
  -2,  -2,   1,   6,  -3,   0,   2,  -1,  -1,  -3,  -4,  -1,  -3,  -3,  -1,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,
   0,  -3,  -3,  -3,   9,  -3,  -4,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -3,  -1,  -1,  -2,  -2,  -1,  -3,  -3,  -2,  -4,
  -1,   1,   0,   0,  -3,   5,   2,  -2,   0,  -3,  -2,   1,   0,  -3,  -1,   0,  -1,  -2,  -1,  -2,   0,   3,  -1,  -4,
  -1,   0,   0,   2,  -4,   2,   5,  -2,   0,  -3,  -3,   1,  -2,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,
   0,  -2,   0,  -1,  -3,  -2,  -2,   6,  -2,  -4,  -4,  -2,  -3,  -3,  -2,   0,  -2,  -2,  -3,  -3,  -1,  -2,  -1,  -4,
  -2,   0,   1,  -1,  -3,   0,   0,  -2,   8,  -3,  -3,  -1,  -2,  -1,  -2,  -1,  -2,  -2,   2,  -3,   0,   0,  -1,  -4,
  -1,  -3,  -3,  -3,  -1,  -3,  -3,  -4,  -3,   4,   2,  -3,   1,   0,  -3,  -2,  -1,  -3,  -1,   3,  -3,  -3,  -1,  -4,
  -1,  -2,  -3,  -4,  -1,  -2,  -3,  -4,  -3,   2,   4,  -2,   2,   0,  -3,  -2,  -1,  -2,  -1,   1,  -4,  -3,  -1,  -4,
  -1,   2,   0,  -1,  -3,   1,   1,  -2,  -1,  -3,  -2,   5,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   0,   1,  -1,  -4,
  -1,  -1,  -2,  -3,  -1,   0,  -2,  -3,  -2,   1,   2,  -1,   5,   0,  -2,  -1,  -1,  -1,  -1,   1,  -3,  -1,  -1,  -4,
  -2,  -3,  -3,  -3,  -2,  -3,  -3,  -3,  -1,   0,   0,  -3,   0,   6,  -4,  -2,  -2,   1,   3,  -1,  -3,  -3,  -1,  -4,
  -1,  -2,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -3,  -3,  -1,  -2,  -4,   7,  -1,  -1,  -4,  -3,  -2,  -2,  -1,  -2,  -4,
   1,  -1,   1,   0,  -1,   0,   0,   0,  -1,  -2,  -2,   0,  -1,  -2,  -1,   4,   1,  -3,  -2,  -2,   0,   0,   0,  -4,
   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -1,  -1,  -1,  -2,  -1,   1,   5,  -2,  -2,   0,  -1,  -1,   0,  -4,
  -3,  -3,  -4,  -4,  -2,  -2,  -3,  -2,  -2,  -3,  -2,  -3,  -1,   1,  -4,  -3,  -2,  11,   2,  -3,  -4,  -3,  -2,  -4,
  -2,  -2,  -2,  -3,  -2,  -1,  -2,  -3,   2,  -1,  -1,  -2,  -1,   3,  -3,  -2,  -2,   2,   7,  -1,  -3,  -2,  -1,  -4,
   0,  -3,  -3,  -3,  -1,  -2,  -2,  -3,  -3,   3,   1,  -2,   1,  -1,  -2,  -2,   0,  -3,  -1,   4,  -3,  -2,  -1,  -4,
  -2,  -1,   3,   4,  -3,   0,   1,  -1,   0,  -3,  -4,   0,  -3,  -3,  -2,   0,  -1,  -4,  -3,  -3,   4,   1,  -1,  -4,
  -1,   0,   0,   1,  -3,   3,   4,  -2,   0,  -3,  -3,   1,  -1,  -3,  -1,   0,  -1,  -3,  -2,  -2,   1,   4,  -1,  -4,
   0,  -1,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,   0,   0,  -2,  -1,  -1,  -1,  -1,  -1,  -4,
  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,  -4,   1
))

BLOSUM80 = (c_int * (24 * 24))(*(
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X   *
   5,  -2,  -2,  -2,  -1,  -1,  -1,   0,  -2,  -2,  -2,  -1,  -1,  -3,  -1,   1,   0,  -3,  -2,   0,  -2,  -1,  -1,  -6,
  -2,   6,  -1,  -2,  -4,   1,  -1,  -3,   0,  -3,  -3,   2,  -2,  -4,  -2,  -1,  -1,  -4,  -3,  -3,  -1,   0,  -1,  -6,
  -2,  -1,   6,   1,  -3,   0,  -1,  -1,   0,  -4,  -4,   0,  -3,  -4,  -3,   0,   0,  -4,  -3,  -4,   5,   0,  -1,  -6,
  -2,  -2,   1,   6,  -4,  -1,   1,  -2,  -2,  -4,  -5,  -1,  -4,  -4,  -2,  -1,  -1,  -6,  -4,  -4,   5,   1,  -1,  -6,
  -1,  -4,  -3,  -4,   9,  -4,  -5,  -4,  -4,  -2,  -2,  -4,  -2,  -3,  -4,  -2,  -1,  -3,  -3,  -1,  -4,  -4,  -1,  -6,
  -1,   1,   0,  -1,  -4,   6,   2,  -2,   1,  -3,  -3,   1,   0,  -4,  -2,   0,  -1,  -3,  -2,  -3,   0,   4,  -1,  -6,
  -1,  -1,  -1,   1,  -5,   2,   6,  -3,   0,  -4,  -4,   1,  -2,  -4,  -2,   0,  -1,  -4,  -3,  -3,   1,   5,  -1,  -6,
   0,  -3,  -1,  -2,  -4,  -2,  -3,   6,  -3,  -5,  -4,  -2,  -4,  -4,  -3,  -1,  -2,  -4,  -4,  -4,  -1,  -3,  -1,  -6,
  -2,   0,   0,  -2,  -4,   1,   0,  -3,   8,  -4,  -3,  -1,  -2,  -2,  -3,  -1,  -2,  -3,   2,  -4,  -1,   0,  -1,  -6,
  -2,  -3,  -4,  -4,  -2,  -3,  -4,  -5,  -4,   5,   1,  -3,   1,  -1,  -4,  -3,  -1,  -3,  -2,   3,  -4,  -4,  -1,  -6,
  -2,  -3,  -4,  -5,  -2,  -3,  -4,  -4,  -3,   1,   4,  -3,   2,   0,  -3,  -3,  -2,  -2,  -2,   1,  -4,  -3,  -1,  -6,
  -1,   2,   0,  -1,  -4,   1,   1,  -2,  -1,  -3,  -3,   5,  -2,  -4,  -1,  -1,  -1,  -4,  -3,  -3,  -1,   1,  -1,  -6,
  -1,  -2,  -3,  -4,  -2,   0,  -2,  -4,  -2,   1,   2,  -2,   6,   0,  -3,  -2,  -1,  -2,  -2,   1,  -3,  -1,  -1,  -6,
  -3,  -4,  -4,  -4,  -3,  -4,  -4,  -4,  -2,  -1,   0,  -4,   0,   6,  -4,  -3,  -2,   0,   3,  -1,  -4,  -4,  -1,  -6,
  -1,  -2,  -3,  -2,  -4,  -2,  -2,  -3,  -3,  -4,  -3,  -1,  -3,  -4,   8,  -1,  -2,  -5,  -4,  -3,  -2,  -2,  -1,  -6,
   1,  -1,   0,  -1,  -2,   0,   0,  -1,  -1,  -3,  -3,  -1,  -2,  -3,  -1,   5,   1,  -4,  -2,  -2,   0,   0,  -1,  -6,
   0,  -1,   0,  -1,  -1,  -1,  -1,  -2,  -2,  -1,  -2,  -1,  -1,  -2,  -2,   1,   5,  -4,  -2,   0,  -1,  -1,  -1,  -6,
  -3,  -4,  -4,  -6,  -3,  -3,  -4,  -4,  -3,  -3,  -2,  -4,  -2,   0,  -5,  -4,  -4,  11,   2,  -3,  -5,  -3,  -1,  -6,
  -2,  -3,  -3,  -4,  -3,  -2,  -3,  -4,   2,  -2,  -2,  -3,  -2,   3,  -4,  -2,  -2,   2,   7,  -2,  -3,  -3,  -1,  -6,
   0,  -3,  -4,  -4,  -1,  -3,  -3,  -4,  -4,   3,   1,  -3,   1,  -1,  -3,  -2,   0,  -3,  -2,   4,  -4,  -3,  -1,  -6,
  -2,  -1,   5,   5,  -4,   0,   1,  -1,  -1,  -4,  -4,  -1,  -3,  -4,  -2,   0,  -1,  -5,  -3,  -4,   5,   0,  -1,  -6,
  -1,   0,   0,   1,  -4,   4,   5,  -3,   0,  -4,  -3,   1,  -1,  -4,  -2,   0,  -1,  -3,  -3,  -3,   0,   5,  -1,  -6,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -6,
  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1
))

BLOSUM90 = (c_int * (24 * 24))(*(
#  A    R    N    D    C    Q    E    G    H    I    L    K    M    F    P    S    T    W    Y    V    B    Z    X    *
   5,  -2,  -2,  -3,  -1,  -1,  -1,   0,  -2,  -2,  -2,  -1,  -2,  -3,  -1,   1,   0,  -4,  -3,  -1,  -2,  -1,  -1,  -6,
  -2,   6,  -1,  -3,  -5,   1,  -1,  -3,   0,  -4,  -3,   2,  -2,  -4,  -3,  -1,  -2,  -4,  -3,  -3,  -2,   0,  -1,  -6,
  -2,  -1,   7,   1,  -4,   0,  -1,  -1,   0,  -4,  -4,   0,  -3,  -4,  -3,   0,   0,  -5,  -3,  -4,   5,  -1,  -1,  -6,
  -3,  -3,   1,   7,  -5,  -1,   1,  -2,  -2,  -5,  -5,  -1,  -4,  -5,  -3,  -1,  -2,  -6,  -4,  -5,   5,   1,  -1,  -6,
  -1,  -5,  -4,  -5,   9,  -4,  -6,  -4,  -5,  -2,  -2,  -4,  -2,  -3,  -4,  -2,  -2,  -4,  -4,  -2,  -4,  -5,  -1,  -6,
  -1,   1,   0,  -1,  -4,   7,   2,  -3,   1,  -4,  -3,   1,   0,  -4,  -2,  -1,  -1,  -3,  -3,  -3,  -1,   5,  -1,  -6,
  -1,  -1,  -1,   1,  -6,   2,   6,  -3,  -1,  -4,  -4,   0,  -3,  -5,  -2,  -1,  -1,  -5,  -4,  -3,   1,   5,  -1,  -6,
   0,  -3,  -1,  -2,  -4,  -3,  -3,   6,  -3,  -5,  -5,  -2,  -4,  -5,  -3,  -1,  -3,  -4,  -5,  -5,  -2,  -3,  -1,  -6,
  -2,   0,   0,  -2,  -5,   1,  -1,  -3,   8,  -4,  -4,  -1,  -3,  -2,  -3,  -2,  -2,  -3,   1,  -4,  -1,   0,  -1,  -6,
  -2,  -4,  -4,  -5,  -2,  -4,  -4,  -5,  -4,   5,   1,  -4,   1,  -1,  -4,  -3,  -1,  -4,  -2,   3,  -5,  -4,  -1,  -6,
  -2,  -3,  -4,  -5,  -2,  -3,  -4,  -5,  -4,   1,   5,  -3,   2,   0,  -4,  -3,  -2,  -3,  -2,   0,  -5,  -4,  -1,  -6,
  -1,   2,   0,  -1,  -4,   1,   0,  -2,  -1,  -4,  -3,   6,  -2,  -4,  -2,  -1,  -1,  -5,  -3,  -3,  -1,   1,  -1,  -6,
  -2,  -2,  -3,  -4,  -2,   0,  -3,  -4,  -3,   1,   2,  -2,   7,  -1,  -3,  -2,  -1,  -2,  -2,   0,  -4,  -2,  -1,  -6,
  -3,  -4,  -4,  -5,  -3,  -4,  -5,  -5,  -2,  -1,   0,  -4,  -1,   7,  -4,  -3,  -3,   0,   3,  -2,  -4,  -4,  -1,  -6,
  -1,  -3,  -3,  -3,  -4,  -2,  -2,  -3,  -3,  -4,  -4,  -2,  -3,  -4,   8,  -2,  -2,  -5,  -4,  -3,  -3,  -2,  -1,  -6,
   1,  -1,   0,  -1,  -2,  -1,  -1,  -1,  -2,  -3,  -3,  -1,  -2,  -3,  -2,   5,   1,  -4,  -3,  -2,   0,  -1,  -1,  -6,
   0,  -2,   0,  -2,  -2,  -1,  -1,  -3,  -2,  -1,  -2,  -1,  -1,  -3,  -2,   1,   6,  -4,  -2,  -1,  -1,  -1,  -1,  -6,
  -4,  -4,  -5,  -6,  -4,  -3,  -5,  -4,  -3,  -4,  -3,  -5,  -2,   0,  -5,  -4,  -4,  11,   2,  -3,  -6,  -4,  -1,  -6,
  -3,  -3,  -3,  -4,  -4,  -3,  -4,  -5,   1,  -2,  -2,  -3,  -2,   3,  -4,  -3,  -2,   2,   8,  -3,  -4,  -3,  -1,  -6,
  -1,  -3,  -4,  -5,  -2,  -3,  -3,  -5,  -4,   3,   0,  -3,   0,  -2,  -3,  -2,  -1,  -3,  -3,   5,  -4,  -3,  -1,  -6,
  -2,  -2,   5,   5,  -4,  -1,   1,  -2,  -1,  -5,  -5,  -1,  -4,  -4,  -3,   0,  -1,  -6,  -4,  -4,   5,   0,  -1,  -6,
  -1,   0,  -1,   1,  -5,   5,   5,  -3,   0,  -4,  -4,   1,  -2,  -4,  -2,  -1,  -1,  -4,  -3,  -3,   0,   5,  -1,  -6,
  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -6,
  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1
))
