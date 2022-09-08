# Hawk implementation

Hawk is a post-quantum signature scheme, which has a hash-then-sign design and is based on the lattice isomorphism problem. The scheme allows efficient signature generation and verification, even on constrained devices. The signatures and public keys are relatively compact compared to prior post-quantum signature schemes.

**Disclaimer:** *this code was written to accompany [1] and currently serves research purposes. There was no security review yet so use at own risk*.

This implementation is written in C. The implementation in avx2-optimized contains another implementation, which is optimized for processors that support both the AVX2 instruction set and floating points.
However, in both cases, the same API is implemented.

This implementation is written in C and is configurable at compile-time
through macros which are documented in `config.h`; each macro is a boolean
option and can be enabled or disabled in `config.h` and/or as a
command-line parameter to the compiler. Several implementation strategies
are available;

The AVX2 optimized version allows the following compilation flags:

- `HAWK_RECOVER_CHECK`: This flag makes sure that the decompression of a signature yields the intended first half of the signature by regenerating a signature when this check fails. However, for `logn = 9` providing NIST-1 security and `logn = 10` providing NIST-5 security, this check is highly unlikely to ever fail.
- `HAWK_FMA`: This flag enables the use for FMA ("fused multiply-add") compiler intrinsics for an extra boost to performance. Occasionally (but rarely), use of HAWK_FMA will change the keys and/or signatures generated from a given random seed, impacting reproducibility of test vectors; however, this has no bearing on the security of normal usage.

## Usage

Type `make` to compile: this will generate a binary called `bin/speed`. `bin/speed` runs performance benchmarks on Hawk for all degrees up to 1024. Note here that Hawk-512 and Hawk-1024 are part of "official" Hawk while smaller degrees are only there for experimenting.

Applications that want to use Hawk normally work on the external API, which is documented in the `hawk.h` file. This is the only file that an external application needs to use. For research purposes, the inner API is documented in `inner.h`. This API gives access to many internal functions that perform some elementary operations used in Hawk. That API also has some non-obvious requirements, such as alignment on temporary buffers, or the need to adjust FPU precision on 32-bit x86 systems.

## Authorship

This project contains parts of [Falcon](https://falcon-sign.info/). A non-complete list is as follows:

- `keygen.c`: the function `solve_NTRU` together with all the functions this function makes calls to are from Falcon, but adopted to `q = 1` and the modified standard deviation used in sampling did therefore also affect the estimation of numbers in intermediate layers of the NTRU generation.
- `fft.c`
- `fpr.c/fpr.h`: the floating point emulation of Falcon is still in use in Hawk for key generation, but signing and verifying do not make use of floating points in the reference implementation.
- `rng.c`
- `shake.c`
- `tests/speed.c`: benchmarking function has stayed the same to have a fair comparison.

The Hawk code was written by Ludo Pulles <ludo.pulles@cwi.nl>. The Falcon code was written by Thomas Pornin <thomas.pornin@nccgroup.com>.

## License

This code is provided under the MIT license.

This software includes a modified version of Falcon, a third party open source software component. See https://falcon-sign.info/ for the original version. Falcon's code is also provided under the MIT license.

## References

[1] Léo Ducas, Eamonn W. Postlethwaite, Ludo N. Pulles, Wessel van Woerden, Hawk: Module LIP makes Lattice Signatures Fast, Compact and Simple ([ePrint](https://eprint.iacr.org/2022/1155))

