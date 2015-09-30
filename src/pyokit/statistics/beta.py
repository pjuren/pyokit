#!/usr/bin/python

"""
  Date of Creation: 7th Apr 2012
  Description: Beta distribution. Some of the beta-function implementation
               adapted from https://malishoaib.wordpress.com

  Copyright (C) 2012-2014
  University of Southern California,
  Philip J. Uren

  Authors: Philip J. Uren

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

# standard python imports
import math

# pyokit imports
from pyokit.common.pyokitError import PyokitError


def beta_1(a, b):
  """use gamma function or inbuilt math.gamma() to compute vals of beta func"""
  beta = math.exp(math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
  return beta


def contfractbeta(a, b, x, ITMAX=200):
  """
  contfractbeta() evaluates the continued fraction form of the incomplete
  Beta function; incompbeta(). (Code translated from: Numerical Recipes in C.)
  """

  EPS = 3.0e-7
  bm = az = am = 1.0
  qab = a + b
  qap = a + 1.0
  qam = a - 1.0
  bz = 1.0 - qab * x / qap

  for i in range(ITMAX + 1):
    em = float(i + 1)
    tem = em + em
    d = em * (b - em) * x / ((qam + tem) * (a + tem))
    ap = az + d * am
    bp = bz + d * bm
    d = -(a + em) * (qab + em) * x / ((qap + tem) * (a + tem))
    app = ap + d * az
    bpp = bp + d * bz
    aold = az
    am = ap / bpp
    bm = bp / bpp
    az = app / bpp
    bz = 1.0
    if (abs(az - aold) < (EPS * abs(az))):
      return az
  raise PyokitError("a or b too large or given ITMAX too small for "
                    "computing incomplete beta function.")


def incompbeta(a, b, x):
  """
  Incomplete beta function; code translated from: Numerical Recipes in C.

  :param a: a > 0
  :param b: b > 0
  :param x: 0 <= x <= 1.
  """

  if (x == 0):
    return 0
  elif (x == 1):
    return 1
  else:
    lbeta = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
             a * math.log(x) + b * math.log(1 - x))
    if (x < (a + 1) / (a + b + 2)):
      return math.exp(lbeta) * contfractbeta(a, b, x) / a
    else:
      return 1 - math.exp(lbeta) * contfractbeta(b, a, 1 - x) / b
