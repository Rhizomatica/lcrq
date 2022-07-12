# C implementation of RFC6330 RaptorQ Codes for Librecast

<a href="https://opensource.org"><img height="150" align="right" src="https://opensource.org/files/OSIApprovedCropped.png" alt="Open Source Initiative Approved License logo"></a>

![Librecast Logo](https://secure.gravatar.com/avatar/52295d18e59ef41aeac21f3745250288?s=200)

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

```
./configure
make
make test # (optional)
sudo make install
```

The code compiles using either gcc or clang.  There is a `make clang` target.
The default is whatever your default CC is.

To best performance, you will want to set some CFLAGS appropriate to your
platform.  Something like:

```
./configure CFLAGS="-O3 -march=native -mpopcnt -pipe -ffast-math -funroll-loops -flto -DNDEBUG"
```
can make a huge difference.

To build with SIMD enabled (highly recommended), append `--enable-simd` to your
configure command. You must enable at least `-mssse3` to use this.
`-march=native` is probably a better bet.

Putting that all together:

```
./configure CFLAGS="-O3 -march=native -mpopcnt -pipe -ffast-math -funroll-loops -flto -DNDEBUG" --enable-simd
make
make install
```

A rough indication of encoding and decoding performance can be garnered by
running  `make speedtest`. This is useful for observing the effects of various
optimizations. It is not a proper benchmark test.

NB: if rebuilding with different compiler options, don't forget to run
`make realclean` first.

### *BSD

Users of *BSD will need to use gmake instead of make. `gmake test` requires
bash.

## Other RaptorQ Implementations

TvRQ
: https://github.com/lorinder/TvRQ (C)

nanorq
: https://github.com/sleepybishop/nanorq (C)

libRaptorQ
: https://github.com/LucaFulchir/libRaptorQ (C++)

raptorq
: https://github.com/cberner/raptorq (Rust)

## References

- https://datatracker.ietf.org/doc/html/rfc6330

- http://www.ece.ubc.ca/~janm/Papers_RG/Shokrollahi_IT_June06.pdf

- https://www.cberner.com/2019/03/30/raptorq-rfc6330-rust-optimization/

- https://www.cberner.com/2020/10/12/building-fastest-raptorq-rfc6330-codec-rust/

- https://github.com/cberner/raptorq/blob/master/RFC6330_ERRATA.md

<hr />

<p class="bigbreak">
This project was funded through the <a href="https://nlnet.nl/discovery"> NGI0 Discovery </a> Fund, a fund established by NLnet with financial support from the European
Commission's <a href="https://ngi.eu">Next Generation Internet</a> programme, under the aegis of DG Communications Networks, Content and Technology under grant agreement No 825322. *Applications are still open, you can <a href="https://nlnet.nl/propose">apply today</a>*
</p>

<p>
  <a href="https://nlnet.nl/project/LibrecastLive/">
      <img width="250" src="https://nlnet.nl/logo/banner.png" alt="Logo NLnet: abstract logo of four people seen from above" class="logocenter" />
  </a>
  <a href="https://ngi.eu/">
      <img width="250" align="right" src="https://nlnet.nl/image/logos/NGI0_tag.png" alt="Logo NGI Zero: letterlogo shaped like a tag" class="logocenter" />
  </a>
</p>
