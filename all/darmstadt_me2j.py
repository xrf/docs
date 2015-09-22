from __future__ import print_function

def iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
    '''Iterator that runs through all 2-body matrix elements in the Darmstadt
    ME2J order.  Yields nlj1, nlj2, nnlj1, nnlj2, J, T, MT per iteration.'''

    # while looping, we must ensure that both nnlj1 >= nnlj2 and
    # (nlj1, nlj2) >= (nnlj1, nnlj2) lexicographically

    # *** bra loops ***
    for nlj1 in range(nljDim):
        for nlj2 in range(nlj1 + 1):

            # ensure that e1 + e2 <= EMax
            if e_nlj[nlj1] + e_nlj[nlj2] > EMax:
                break

            # *** ket loops ***
            for nnlj1 in range(nlj1 + 1):
                for nnlj2 in range((nlj2 if nnlj1 == nlj1 else nnlj1) + 1):

                    # ensure that ee1 + ee2 <= EMax
                    if e_nlj[nnlj1] + e_nlj[nnlj2] > EMax:
                        break

                    # ensure parity conservation
                    if (l_nlj[nlj1] + l_nlj[nlj2]
                        - l_nlj[nnlj1] - l_nlj[nnlj2]) % 2:
                        continue

                    # compute the range of J allowed by triangular condition
                    JMin = max(abs(twoj_nlj[nlj1] - twoj_nlj[nlj2]),
                               abs(twoj_nlj[nnlj1] - twoj_nlj[nnlj2])) // 2
                    JMax = min((twoj_nlj[nlj1] + twoj_nlj[nlj2]),
                               (twoj_nlj[nnlj1] + twoj_nlj[nnlj2])) // 2

                    # *** diagonal loops ***

                    # loop is automatically skipped if JMin > JMax
                    for J in range(JMin, JMax + 1):
                        for T in range(2):
                            for MT in range(-T, T + 1):

                                # antisymmetry forbidden states (swap phase
                                # for nlj1 == nlj2 must be +1, same for
                                # nnlj1 == nnlj2) are not removed in order to
                                # have constant-size T-MT-blocks

                                yield nlj1, nlj2, nnlj1, nnlj2, J, T, MT

def iterate_nlj(eMax, nMax, lMax):
    '''Iterator that runs through n, l, j of all single-particle states.
    Yields e, n, l, twoj per iteration.'''

    for e in range(eMax + 1):
        lMin = e % 2
        for l in range(lMin, min(e, lMax) + 1, 2):

            n = (e - l) // 2
            if n > nMax:
                continue

            twojMin = abs(2 * l - 1)
            twojMax = 2 * l + 1
            for twoj in range(twojMin, twojMax + 1, 2):
                yield e, n, l, twoj

def nlj_table(eMax, nMax, lMax):
    '''Return e_nlj, n_nlj, l_nlj, twoj_nlj, which are arrays that map
    nlj to e, n, l, twoj respectively.'''
    # 'zip' acts like a transposition here
    return tuple(zip(*iterate_nlj(eMax, nMax, lMax)))

def get_nljDim(eMax, nMax, lMax):
    return sum(1 for _ in iterate_nlj(eMax, nMax, lMax))

def get_me2jDim(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
    return sum(1 for _ in iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj))

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

def print_row(*args, **kwargs):
    import sys
    stream = kwargs.get("stream", sys.stdout)
    stream.write("\t".join(map(str, args)) + "\n")

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

def basic_example():
    import sys

    filename = "darmstadt-me2j-quantum-numbers.txt"
    eMax = 5
    nMax = 99999
    lMax = 99999
    EMax = 99999
    nljDim = get_nljDim(eMax, nMax, lMax)

    # print the single-particle states
    print_row("nlj", "e", "n", "l", "j")
    for nlj, (e, n, l, twoj) in enumerate(iterate_nlj(eMax, nMax, lMax)):
        print_row(nlj, e, n, l, twoj / 2.)
        assert unpack_nlj(nlj) == (e, n, l, twoj)
        assert pack_nlj(e, twoj) == nlj
    print("Total of {0} single-particle state(s).".format(nljDim))

    # build a lookup table for converting nlj to quantum numbers
    e_nlj, n_nlj, l_nlj, twoj_nlj = nlj_table(eMax, nMax, lMax)

    # write the quantum numbers out to file
    print((
        "\n" +
        "Quantum numbers of 2-body matrix elements will be written to:\n" +
        "  {0}\n".format(filename) +
        "Writing... "
    ), end="")
    sys.stdout.flush()
    with open(filename, "w") as stream:
        print_row("nlj1", "nlj2", "nnlj1", "nnlj2",
                  "J", "T", "MT", stream=stream)
        for qnums in iterate_me2j(nljDim, EMax, e_nlj, l_nlj, twoj_nlj):
            print_row(*qnums, stream=stream)
    print("done.")

if __name__ == "__main__":
    basic_example()
