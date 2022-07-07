# C implementation of RFC6330 RaptorQ Codes for Librecast

## Documentation

See `man lcrq(7)`

## Installation

```
./configure
make
make test # (optional)
make install
```

To build with SIMD enabled, use:
`CFLAGS="-DINTEL_SSE3 " make`

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
