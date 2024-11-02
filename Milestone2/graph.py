import matplotlib.pyplot as plt


# Format is:
# { .... a bunch of lines .... }
# 1000 multiplications of size <size>. Total time <time> us
def parse_out(filepath: str) -> [int]:
    results = []
    with open(filepath) as f:
        for line in f:
            if line[:20] != "1000 multiplications":
                continue
            results.append(int(line.split(' ')[-2])/(1e6))
    return results


sizes = [ i for i in range(9, 16) ]

# Benchmarks IN ORDER
serial = parse_out("slurm_out/initial/initial.out")
avx_loop_bodies = parse_out("slurm_out/avx/all_vectorised.out")
avx_all_vectors = parse_out("slurm_out/avx/fft_base_cases.out")
omp = parse_out("slurm_out/omp/max_depth_2_min_8.out")


def compareSpeedup(this, base):
    return [ b/t for (t, b) in zip(this, base) ]


def init():
    plt.plot(sizes, serial)
    plt.title("Performance of Serial Multiplication Implementation")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("average time per mutliplication (milliseconds)")
    plt.show()


def avx():
    speedup_loop_bodies = compareSpeedup(avx_loop_bodies, serial)
    speedup_all_vectors = compareSpeedup(avx_all_vectors, serial)
    plt.plot(sizes, speedup_loop_bodies, marker="o")
    plt.plot(sizes, speedup_all_vectors, marker="o")
    for (x, y, time) in zip(sizes, speedup_all_vectors, avx_all_vectors):
        plt.annotate(str(round(time, 2))+"ms", (x, y-0.05))

    plt.title("Relative Performance of SIMD Implementations")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Multiple of speed of serial implementation")
    plt.legend(["only loops vectorised", "all vectorised"])
    plt.show()


def openMP():
    speedup_all_vectors = compareSpeedup(avx_all_vectors, serial)
    speedup_omp = compareSpeedup(omp, serial)
    plt.plot(sizes, speedup_all_vectors, marker="o")
    plt.plot(sizes, speedup_omp, marker="o")
    for (x, y, time) in zip(sizes, speedup_omp, omp):
        plt.annotate(str(round(time, 2))+"ms", (x, y-0.15))

    plt.title("Relative Performance of SIMD and OpenMP")
    plt.xlabel("log2 of the degree of the operand polynomials")
    plt.ylabel("Multiple of speed of serial implementation")
    plt.legend(["SIMD only", "SIMD and OpenMP"])
    plt.show()


openMP()
