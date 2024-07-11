# C implementation of RFC6330 RaptorQ Codes for Librecast

<a href="https://opensource.org"><img height="150" align="right" src="https://opensource.org/files/OSIApprovedCropped.png" alt="Open Source Initiative Approved License logo"></a>

![Librecast Logo](https://secure.gravatar.com/avatar/52295d18e59ef41aeac21f3745250288?s=200)

<a href="https://scan.coverity.com/projects/lcrq">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/25261/badge.svg"/>
</a>

IP  Multicast  is based on UDP, which is inherently unreliable. Packets may
arrive out of order, or not at all. TCP provides unicast with a reliable
messaging layer on top of this unreliable, connectionless medium.

Unicast,  however, is one-to-one only. Multicast could, in theory, use all of
the same reliability options (ACKs etc.) as TCP, at the cost of not being
scalable any more.

Fortunately there are other ways to achieve similar reliability.  RFC3208
describes Pragmatic General Multicast (PGM) based on NAKs (negative
acknowledgements). This, too, has scaling issues.

Forwards Error Correction (FEC) offers us another approach.

Thanks to parity checking in the network stack, we don't generally need to worry
about errors within packets. Every packet has a checksum, and if that doesn't
match, the packet is dropped before it reaches us. Our encryption provides
further checking of  data  received.   We  need only concern ourselves with
erasures.  ie. dropped packets.

RaptorQ is an implementation of a class of systematic erasure codes called
fountain codes.

The  data we want to send is split into blocks, and then pre-encoded into a set
of intermediate symbols.  From these intermediate symbols we can generate both
our original  source  symbols, and also additional repair symbols.

Provided  the  recipient  receives  at least a minimum value K' of these
symbols (any unique combination of source and repair) the intermediate symbols
can  be  reconstituted,  and  the original data recovered.

RaptorQ  is what is called a systematic encoding, because the set of symbols we
send includes our original data as plain text. Provided all source symbols are
received, the original  data has been transmitted with no decoding overhead.  It
is only in the case where we need to supplement the source symbols with repair
symbols that we must perform the decoding process.

## Documentation

See the man pages:

- `lcrq(7)`
- `rq_init(3)`
- `rq_free(3)`
- `rq_query(3)`
- `rq_encode(3)`
- `rq_decode(3)`
- `rq_symbol(3)`

## Installation

Oh good. You're one of those people who reads documentation. Phew! Please, do
NOT just build the sources with the usual `configure && make`. It will run like a
dog with no legs.

C convention says that I absolutely MUST NOT mess with CFLAGS and try to help
you in any way by setting sane defaults.  That means you need to do it.

How you build LCRQ depends on whether you are building it to just run on your
local machine, or if you are building binaries that need to support multiple
CPUs, such as when building packages for a binary distribution.

There should be no compiler warnings, even with `-Wall -Wextra -pedantic`.
Please raise a bug if you see any.

### Supported Compiliers

The code compiles and has been tested with gcc and clang.  The default is
whatever your default CC is, which you can set when calling `configure`. eg.
`./configure CC=clang`.

### Configure Options

#### Native Compilation (best performance)

LCRQ makes heavy use of SIMD instructions when available. For best performance,
a native build targeting your local machine is best.

If you are compiling lcrq for use only on your local machine, you will want to
set some CFLAGS appropriate to your platform. You can read more about these
options in the [GCC manual](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html).

Consider:

`-O3`
: Optimize aggressively. Highly recommended.

`-march=native`
: Enable all CPU instruction sets available on the local machine, and tuning for
your CPU. Combine this with `--enable-native` for best results.

`-mpopcnt`
: In theory this will be enabled by `-march=native`, but it doesn't hurt to be
explicit if your CPU supports it.

`--enable-native`
: Enable native optimizations in the LCRQ build. Combine supported instruction
sets into one compilation unit, rather than building separately. NB: Requires that
your CPU supports AVX2.

`-flto`
: Enable link-time optimization to allow compiler optimization across
compilation units. Generally safe to use and highly recommended.

`-funroll-loops`
: Give the compiler a bit more encouragement to unroll loops.

`-ffast-math`
: Recommended, but might be unsafe in some circumstances. Try it. If you have
problems, turn it off.

`-DNDEBUG`
: Disable assertions and other debugging code.

So, putting that all together, we get something like:

```
./configure CFLAGS="-Wall -Wextra -pedantic -O3 -march=native -mpopcnt -ffast-math -funroll-loops -flto -DNDEBUG" --enable-native
make
make test      # make sure all tests PASS
make install
```

A rough indication of encoding and decoding performance can be garnered by
running  `make speedtest`. This is useful for observing the effects of various
optimizations. It is not a proper benchmark test.

NB: if rebuilding with different compiler options, don't forget to run
`make clean` first.

#### Portable Compilation (for binary distributions)

If you are packaging LCRQ, or otherwise building it to run on a whole family of
CPUs (like x86\_64) then some of the optimizations in the previous section
aren't applicable.  You can probably safely use:

`-O3`
: Optimize aggressively. Highly recommended.

`-flto`
: Enable link-time optimization to allow compiler optimization across
compilation units. Generally safe to use and highly recommended.

`-funroll-loops`
: Give the compiler a bit more encouragement to unroll loops.

`-ffast-math`
: Recommended, but might be unsafe in some circumstances. Try it. If you have
problems, turn it off.

`-DNDEBUG`
: Disable assertions and other debugging code.

So, putting that all together, we get something like:

```
./configure CFLAGS="-Wall -Wextra -pedantic -O3 -flto -funroll-loops -ffast-math -DNDEBUG"
make
make test      # make sure all tests PASS
make install
```

Please run `make speedtest` in your build scripts and pay attention to the
output. Does it seem reasonable for the target system? Compare it to a build on
your own machine.


### \*BSD

Users of \*BSD will need to use gmake instead of make. `gmake test` requires
bash.

We test all releases on FreeBSD, NetBSD and OpenBSD.

## Other RaptorQ Implementations

TvRQ
: https://github.com/lorinder/TvRQ (C)

nanorq
: https://github.com/sleepybishop/nanorq (C)

libRaptorQ
: https://github.com/LucaFulchir/libRaptorQ (C++)

raptorq
: https://github.com/cberner/raptorq (Rust)

## Acknowledgements

LCRQ uses the fast Galois field multiplication technique described in:

 J. S. Plank and K. M. Greenan and E. L. Miller (2013)
 "Screaming Fast Galois Field Arithmetic Using Intel SIMD Instructions"
 http://web.eecs.utk.edu/~jplank/plank/papers/FAST-2013-GF.html

## References

- https://datatracker.ietf.org/doc/html/rfc6330

- http://www.ece.ubc.ca/~janm/Papers_RG/Shokrollahi_IT_June06.pdf

- https://www.cberner.com/2019/03/30/raptorq-rfc6330-rust-optimization/

- https://www.cberner.com/2020/10/12/building-fastest-raptorq-rfc6330-codec-rust/

- https://github.com/cberner/raptorq/blob/master/RFC6330_ERRATA.md

<hr />

## Funding

This project received funding through [NGI Assure](https://nlnet.nl/assure), a fund established by [NLnet](https://nlnet.nl) with financial support from the European Commission's [Next Generation Internet](https://ngi.eu) program. Learn more at the [NLnet project page](https://nlnet.nl/project/LibreCastLiveStudio).

This project received funding through [NGI0 Discovery](https://nlnet.nl/discovery), a fund established by [NLnet](https://nlnet.nl) with financial support from the European Commission's [Next Generation Internet](https://ngi.eu) program. Learn more at the [NLnet project page](https://nlnet.nl/project/LibrecastLive).

<p>
  <a href="https://nlnet.nl/project/LibreCastLiveStudio/">
      <img width="250" src="https://nlnet.nl/logo/banner.png" alt="Logo NLnet: abstract logo of four people seen from above" class="logocenter" />
  </a>
  <a href="https://ngi.eu/">
      <img width="250" align="right" src="https://nlnet.nl/image/logos/NGIAssure_tag.svg" alt="Logo NGI Assure letterlogo shaped like a tag" class="logocenter" />
  </a>
</p>
