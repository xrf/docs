#!/usr/bin/env python

def range_twoJ(twoj1, twoj2):
    '''Find the range of J formed by coupling j1 with j2.'''
    twoJMin = abs(twoj1 - twoj2)
    twoJMax = twoj1 + twoj2
    return twoJMin, twoJMax

def range_twoJ_both(twoj1, twoj2, twojj1, twojj2):
    '''Find the common range of J formed by coupling j1 with j2 and coupling
    jj1 with jj2.'''
    twoJMin, twoJMax = range_twoJ(twoj1, twoj2)
    twoJJMin, twoJJMax = range_twoJ(twojj1, twojj2)
    return max(twoJMin, twoJJMin), min(twoJMax, twoJJMax)

def range_J_both(twoj1, twoj2, twojj1, twojj2):
    '''Find the common range of J formed by coupling j1 with j2 and coupling
    jj1 with jj2, where j1 and j2 must be either both half-integers or both
    integers, and similarly for jj1 and jj2.'''
    twoJMin, twoJMax = range_twoJ_both(twoj1, twoj2, twojj1, twojj2)
    return twoJMin // 2, twoJMax // 2

def iterate_nlj(eMax, nMax, lMax):
    '''Iterator that runs through n, l, j of all single-particle states.
    Yields e, n, l, twoj per iteration.'''
    for e in range(eMax + 1):
        lMin = e % 2
        for l in range(lMin, min(e, lMax) + 1, 2):
            n = (e - l) // 2
            if n > nMax:
                continue
            twojMin, twojMax = range_twoJ(2 * l, 1)
            for twoj in range(twojMin, twojMax + 1, 2):
                yield e, n, l, twoj

def get_nljDim(eMax, nMax, lMax):
    return sum(1 for _ in iterate_nlj(eMax, nMax, lMax))

def nlj_table(eMax, nMax, lMax):
    '''Return e_nlj, n_nlj, l_nlj, twoj_nlj, which are arrays that map
    nlj to e, n, l, twoj respectively.'''
    # 'zip' acts like a transposition here
    return tuple(zip(*iterate_nlj(eMax, nMax, lMax)))

def shell_number(n, l):
    return 2 * n + l

def pack_nlj(e, twoj):
    return (e * (e + 1) + twoj) // 2

def unpack_nlj(nlj):
    from math import sqrt
    e = (int(sqrt(8 * nlj + 1)) - 1) // 2
    erem = nlj - e * (e + 1) // 2
    twoj = 1 + erem * 2
    n = (e - erem) // 2
    l = e - 2 * n
    return e, n, l, twoj

def violates_parity_2(nlj1, nlj2, nnlj1, nnlj2, l_nlj):
    return (l_nlj[nlj1] + l_nlj[nlj2] - l_nlj[nnlj1] - l_nlj[nnlj2]) % 2

def iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
    '''Iterator that runs through all 2-body matrix elements in the Darmstadt
    ME2J order.  Yields nlj1, nlj2, nnlj1, nnlj2, J, T, MT per iteration.'''
    for nlj1 in range(nljDim):
        for nlj2 in range(nlj1 + 1):
            if e_nlj[nlj1] + e_nlj[nlj2] > EMax:
                break
            for nnlj1 in range(nlj1 + 1):
                for nnlj2 in range((nlj2 if nnlj1 == nlj1 else nnlj1) + 1):
                    if e_nlj[nnlj1] + e_nlj[nnlj2] > EMax:
                        break
                    if violates_parity_2(nlj1, nlj2, nnlj1, nnlj2, l_nlj):
                        continue
                    JMin, JMax = range_J_both(twoj_nlj[nlj1], twoj_nlj[nlj2],
                                              twoj_nlj[nnlj1], twoj_nlj[nnlj2])
                    for J in range(JMin, JMax + 1):
                        for T in range(2):
                            for MT in range(-T, T + 1):
                                yield nlj1, nlj2, nnlj1, nnlj2, J, T, MT

def get_me2jDim(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
    return sum(1 for _ in iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj))

def read_numbers(stream):
    for line in stream:
        for num in line.split():
            yield float(num)

def read_me2j(filename, eMax, nMax, lMax, EMax):
    '''Read an .me2j matrix element text file.  Returns a numpy array.'''
    import numpy

    nljDim = get_nljDim(eMax, nMax, lMax)
    e_nlj, n_nlj, l_nlj, twoj_nlj = nlj_table(eMax, nMax, lMax)

    me2jDim = get_me2jDim(nljDim, EMax, e_nlj, l_nlj, twoj_nlj)
    matrixElems = numpy.empty(me2jDim)

    with open(filename) as stream:
        next(stream)                    # skip one line
        for i, num in zip(range(me2jDim), read_numbers(stream)):
            matrixElems[i] = num

    if i < me2jDim:
        raise ValueError("file has too few numbers")
    return matrixElems

def print_row(*args, **kwargs):
    import sys
    stream = kwargs.get("stream", sys.stdout)
    stream.write("\t".join(map(str, args)) + "\n")

def parse_args(argv):
    import argparse
    HUGE_VALUE = 2 ** 31 - 1
    p = argparse.ArgumentParser(
        description="note: EMax/lMax/nMax default to a huge value",
    )
    p.add_argument(
        "-o",
        dest="file",
        metavar="file",
        default="darmstadt-me2j-quantum-numbers.txt",
        help="Output filename",
    )
    p.add_argument(
        "-E",
        dest="EMax",
        metavar="EMax",
        type=int,
        default=HUGE_VALUE,
        help="Maximum combined shell number",
    )
    p.add_argument(
        "-l",
        dest="lMax",
        metavar="lMax",
        type=int,
        default=HUGE_VALUE,
        help="Maximum single-particle angular momentum quantum number",
    )
    p.add_argument(
        "-n",
        dest="nMax",
        metavar="nMax",
        type=int,
        default=HUGE_VALUE,
        help="Maximum single-particle principal quantum number",
    )
    p.add_argument(
        dest="eMax",
        metavar="eMax",
        type=int,
        help="Maximum single-particle shell number",
    )
    args = p.parse_args(argv[1:])
    return args.file, args.eMax, args.nMax, args.lMax, args.EMax

def main():
    import sys

    filename, eMax, nMax, lMax, EMax = parse_args(sys.argv)

    print("\n".join(["Input limits:",
                     "  eMax = {0}".format(eMax),
                     "  nMax = {0}".format(nMax),
                     "  lMax = {0}".format(lMax),
                     "  EMax = {0}\n".format(EMax)]))

    nljDim = get_nljDim(eMax, nMax, lMax)

    # print the single-particle states
    print_row("nlj", "e", "n", "l", "j")
    for nlj, (e, n, l, twoj) in enumerate(iterate_nlj(eMax, nMax, lMax)):
        print_row(nlj, e, n, l, twoj / 2.)
        e2, _, _, twoj2 = unpack_nlj(nlj)
        assert pack_nlj(e2, twoj2) == nlj
    print("Total of {0} single-particle state(s).\n".format(nljDim))

    # build a lookup table for converting nlj to quantum numbers
    e_nlj, n_nlj, l_nlj, twoj_nlj = nlj_table(eMax, nMax, lMax)

    # write the quantum numbers out to file
    sys.stdout.write(
        "Quantum numbers of 2-body matrix elements will be written to:\n" +
        "  {0}\nWriting... ".format(filename)
    )
    sys.stdout.flush()

    with open(filename, "w") as stream:
        print_row("nlj1", "nlj2", "nnlj1", "nnlj2",
                  "J", "T", "MT", stream=stream)
        for qnums in iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
            print_row(*qnums, stream=stream)
    print("done.")

if __name__ == "__main__":
    main()
