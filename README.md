# C implementation of RFC6330 RaptorQ Codes for Librecast

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
make install
```

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
