import matplotlib.pyplot as plt


# Format is:
# 1000 multiplications of size <size>. Total time <time> us
def parse_out(filepath: str) -> [int]:
    results = []
    with open(filepath) as f:
        results = [ int(line.split(' ')[-2])/(1e6) for line in f ]
    return results


sizes = [ i for i in range(9, 16) ]

# Benchmarks IN ORDER
initial = parse_out("initial_O3/init.stdout")
fast_inverse = parse_out("fast_inverse/fastinv.stdout")
blocking = parse_out("blocking/init.stdout")
blocking_unrolled = parse_out("blocking_unrolled/init.stdout")
fft_unroll = parse_out("fft_unroll/init.stdout")


def relativeTo(this, to):
    return [ a/b for (a,b) in zip(this, to) ]


def absolute():
    plt.plot(sizes, initial)
    plt.plot(sizes, fast_inverse)
    plt.plot(sizes, blocking)
    plt.plot(sizes, blocking_unrolled)
    plt.plot(sizes, fft_unroll)
    plt.title("Performance of Polynomial Multiplication Implementations")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled", "fft unrolled"])
    plt.show()


def relative():
    plt.plot(sizes, relativeTo(initial, initial))
    plt.plot(sizes, relativeTo(fast_inverse, initial))
    plt.plot(sizes, relativeTo(blocking, initial))
    plt.plot(sizes, relativeTo(blocking_unrolled, initial))
    plt.plot(sizes, relativeTo(fft_unroll, initial))
    plt.title("Relative Performance of Polynomial Multiplication Implementations")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Time taken (relative to initial implementation)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled", "fft unrolled"])
    plt.show()


def init():
    plt.plot(sizes, initial)
    plt.title("Performance of Initial Multiplication Implementation")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.show()


def fast_inv():
    plt.subplot(1, 2, 1)
    plt.plot(sizes, initial)
    plt.plot(sizes, fast_inverse)
    plt.title("Absolute Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.legend(["initial", "fast inverse algorithm"])

    plt.subplot(1, 2, 2)
    plt.plot(sizes, relativeTo(initial, initial))
    plt.plot(sizes, relativeTo(fast_inverse, initial))
    plt.title("Relative Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Time taken (relative to initial implementation)")
    plt.legend(["initial", "fast inverse algorithm"])

    plt.show()


def mult_unroll():
    plt.subplot(1, 2, 1)
    plt.plot(sizes, initial)
    plt.plot(sizes, fast_inverse)
    plt.plot(sizes, blocking)
    plt.plot(sizes, blocking_unrolled)
    plt.title("Absolute Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled"])

    plt.subplot(1, 2, 2)
    plt.plot(sizes, relativeTo(initial, initial))
    plt.plot(sizes, relativeTo(fast_inverse, initial))
    plt.plot(sizes, relativeTo(blocking, initial))
    plt.plot(sizes, relativeTo(blocking_unrolled, initial))
    plt.title("Relative Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Time taken (relative to initial implementation)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled"])

    plt.show()


def fft_blocking():
    plt.subplot(1, 2, 1)
    plt.plot(sizes, initial)
    plt.plot(sizes, fast_inverse)
    plt.plot(sizes, blocking)
    plt.plot(sizes, blocking_unrolled)
    plt.plot(sizes, fft_unroll)
    plt.title("Absolute Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled", "fft blocking"])

    plt.subplot(1, 2, 2)
    plt.plot(sizes, relativeTo(initial, initial))
    plt.plot(sizes, relativeTo(fast_inverse, initial))
    plt.plot(sizes, relativeTo(blocking, initial))
    plt.plot(sizes, relativeTo(blocking_unrolled, initial))
    plt.plot(sizes, relativeTo(fft_unroll, initial))
    plt.title("Relative Performance")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Time taken (relative to initial implementation)")
    plt.legend(["initial", "fast inverse algorithm", "blocking", "blocking unrolled", "fft blocking"])

    plt.show()


fft_blocking()
